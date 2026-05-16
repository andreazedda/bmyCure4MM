from __future__ import annotations

import json
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Any

from django.conf import settings
from django.db.models import Count

from clinic.models import Assessment, Patient, PatientTherapy

from .models import (
    AdverseEvent,
    CausalAssumptionSet,
    CounterfactualRun,
    LongitudinalLabResult,
    ObservationResidual,
    PatientTwinState,
    SimulationRunMetadata,
    TherapyInterruption,
)
from .privacy import check_path_payload, is_sensitive_path, scan_text_for_sensitive_markers


def build_check(status: str, label: str, detail: str, *, object_ids=None, next_action: str = "", raw_details=None) -> dict[str, Any]:
    return {
        "status": status,
        "label": label,
        "detail": detail,
        "object_ids": list(object_ids or []),
        "next_action": next_action,
        "raw_details": raw_details or {},
    }


def run_developer_checks(patient: Patient | None = None) -> dict[str, list[dict[str, Any]]]:
    return {
        "data": run_data_consistency_checks(patient),
        "model": run_model_consistency_checks(patient),
        "causal": run_causal_checks(patient),
        "scientific": run_scientific_reference_checks(),
        "privacy": run_privacy_security_checks(),
        "artifact": run_artifact_checks(patient),
    }


def run_data_consistency_checks(patient: Patient | None = None) -> list[dict[str, Any]]:
    patients = Patient.objects.filter(pk=patient.pk) if patient else Patient.objects.all()
    patient_ids = list(patients.values_list("id", flat=True))
    checks: list[dict[str, Any]] = []

    orphan_assessments = Assessment.objects.filter(patient__isnull=True).values_list("id", flat=True)
    checks.append(build_check("pass" if not orphan_assessments else "fail", "orphan assessments", "Assessments should belong to a patient.", object_ids=orphan_assessments, next_action="Attach or delete orphan assessments."))

    therapies_without_dose = PatientTherapy.objects.filter(patient_id__in=patient_ids).filter(doses={}).values_list("id", flat=True)
    checks.append(build_check("pass" if not therapies_without_dose else "warn", "therapies without dose schedule", "Therapies need structured dose schedules for interpretable exposure comparisons.", object_ids=therapies_without_dose, next_action="Backfill PatientTherapy.doses, cycle_length_days, and days_on."))

    interruptions_without_match = TherapyInterruption.objects.filter(patient_id__in=patient_ids, patient_therapy__isnull=True).values_list("id", flat=True)
    checks.append(build_check("pass" if not interruptions_without_match else "warn", "interruptions without matching therapy", "Interruptions should link to a therapy when possible.", object_ids=interruptions_without_match, next_action="Resolve patient_therapy links where source records allow."))

    adverse_without_date = AdverseEvent.objects.filter(patient_id__in=patient_ids, date__isnull=True).values_list("id", flat=True)
    checks.append(build_check("pass" if not adverse_without_date else "fail", "adverse events without date", "Adverse events must be dated to support timeline overlays.", object_ids=adverse_without_date, next_action="Insert event dates or mark as unavailable outside the structured event model."))

    labs_without_unit = LongitudinalLabResult.objects.filter(patient_id__in=patient_ids, unit="").values_list("id", flat=True)
    checks.append(build_check("pass" if not labs_without_unit else "warn", "labs without unit", "Structured lab values should include units.", object_ids=labs_without_unit, next_action="Backfill units for lab rows."))

    duplicate_labs = find_duplicate_lab_records(patient)
    checks.append(build_check("pass" if not duplicate_labs else "fail", "duplicate lab records", "There should be one lab record per patient/date/analyte.", object_ids=[item["ids"] for item in duplicate_labs], next_action="Merge duplicates or preserve only the curated record.", raw_details={"duplicates": duplicate_labs}))

    missing_biomarkers = []
    for item in patients:
        marker_count = LongitudinalLabResult.objects.filter(
            patient=item,
            analyte__in=[LongitudinalLabResult.ANALYTE_M_PROTEIN, LongitudinalLabResult.ANALYTE_FLC_RATIO, LongitudinalLabResult.ANALYTE_HB],
        ).exclude(value__isnull=True).count()
        if marker_count == 0:
            missing_biomarkers.append(item.id)
    checks.append(build_check("pass" if not missing_biomarkers else "warn", "missing key biomarkers by patient", "Patients need M-protein, FLC ratio, or hemoglobin for model fit panels.", object_ids=missing_biomarkers, next_action="Insert longitudinal lab records or assessments with modeled markers."))
    return checks


def find_duplicate_lab_records(patient: Patient | None = None) -> list[dict[str, Any]]:
    queryset = LongitudinalLabResult.objects.all()
    if patient is not None:
        queryset = queryset.filter(patient=patient)
    duplicates = []
    groups = queryset.values("patient_id", "date", "analyte").annotate(count=Count("id")).filter(count__gt=1)
    for group in groups:
        ids = list(queryset.filter(patient_id=group["patient_id"], date=group["date"], analyte=group["analyte"]).values_list("id", flat=True))
        duplicates.append({**group, "ids": ids})
    return duplicates


def run_model_consistency_checks(patient: Patient | None = None) -> list[dict[str, Any]]:
    patients = [patient] if patient is not None else list(Patient.objects.all())
    checks: list[dict[str, Any]] = []
    missing_current = []
    multi_current = []
    no_pre_post = []
    no_improvement = []
    underdetermined = []
    missing_trajectory = []
    missing_predicted = []
    missing_utility = []

    for item in patients:
        current_count = PatientTwinState.objects.filter(patient=item, is_current=True).count()
        if current_count == 0:
            missing_current.append(item.id)
        if current_count > 1:
            multi_current.append(item.id)
        current_state = PatientTwinState.objects.filter(patient=item, is_current=True).first()
        if current_state:
            pre = ObservationResidual.objects.filter(patient=item, stage=ObservationResidual.STAGE_PRE_CALIBRATION).count()
            post = ObservationResidual.objects.filter(patient=item, stage=ObservationResidual.STAGE_POST_CALIBRATION).count()
            if pre == 0 or post == 0:
                no_pre_post.append(item.id)
            diagnostics = current_state.parameter_uncertainty or {}
            if diagnostics.get("rmse_after") is not None and diagnostics.get("rmse_before") is not None:
                if float(diagnostics["rmse_after"]) >= float(diagnostics["rmse_before"]):
                    no_improvement.append(item.id)
            if int(diagnostics.get("n_observations") or 0) < int(diagnostics.get("n_parameters") or 0) + 1:
                underdetermined.append(item.id)
        for run in CounterfactualRun.objects.filter(patient=item, status=CounterfactualRun.STATUS_COMPLETED):
            if not run.trajectory_artifact:
                missing_trajectory.append(run.id)
            if not (run.simulation_summary or {}).get("predicted_biomarkers"):
                missing_predicted.append(run.id)
            if (run.comparison_metrics or {}).get("research_utility") is None:
                missing_utility.append(run.id)

    checks.append(build_check("pass" if not missing_current else "fail", "current twin state exists", "Each active research patient should have a current twin state.", object_ids=missing_current, next_action="Initialize or calibrate a twin state."))
    checks.append(build_check("pass" if not multi_current else "fail", "one current twin state", "A patient should have exactly one current twin state.", object_ids=multi_current, next_action="Repair current-state flags."))
    checks.append(build_check("pass" if not no_pre_post else "warn", "pre/post residuals exist", "Calibration should persist both pre- and post-calibration residuals.", object_ids=no_pre_post, next_action="Rerun calibration."))
    checks.append(build_check("pass" if not no_improvement else "warn", "RMSE improves after calibration", "Calibration should improve or explicitly mark the state unreliable.", object_ids=no_improvement, next_action="Inspect calibration diagnostics."))
    checks.append(build_check("pass" if not underdetermined else "warn", "observations exceed parameters", "Calibration has better support when observations >= parameters + 1.", object_ids=underdetermined, next_action="Add observations or reduce fitted parameters."))
    checks.append(build_check("pass" if not missing_trajectory else "warn", "trajectory artifacts exist", "Completed runs should have trajectory artifacts for comparison charts.", object_ids=missing_trajectory, next_action="Re-run affected scenarios."))
    checks.append(build_check("pass" if not missing_predicted else "warn", "predicted biomarkers exist", "Completed runs should expose predicted biomarker summaries.", object_ids=missing_predicted, next_action="Re-run affected scenarios with current bridge."))
    checks.append(build_check("pass" if not missing_utility else "warn", "utility computed", "Completed runs should include heuristic research utility.", object_ids=missing_utility, next_action="Re-run or backfill comparison metrics."))
    collapse = detect_schedule_collapse(patient)
    checks.append(build_check("pass" if not collapse else "warn", "schedule-resolution collapse", "Some scenarios produce indistinguishable trajectories under current exposure resolution.", object_ids=[item["run_ids"] for item in collapse], next_action="Increase schedule resolution in the exposure bridge.", raw_details={"pairs": collapse}))
    return checks


def detect_schedule_collapse(patient: Patient | None = None, *, tolerance: float = 1.0e-8) -> list[dict[str, Any]]:
    runs = list(_latest_completed_runs(patient).values())
    fingerprints: dict[str, list[CounterfactualRun]] = defaultdict(list)
    for run in runs:
        trajectory = _load_trajectory(run)
        if not trajectory:
            continue
        alternative = trajectory.get("alternative_trajectory") or {}
        tumor = alternative.get("tumor_cells") or []
        healthy = alternative.get("healthy_cells") or []
        if not tumor or not healthy:
            continue
        key = json.dumps({
            "tumor_end": round(float(tumor[-1]), 6),
            "healthy_end": round(float(healthy[-1]), 6),
            "tumor_mid": round(float(tumor[len(tumor)//2]), 6),
            "healthy_mid": round(float(healthy[len(healthy)//2]), 6),
        }, sort_keys=True)
        fingerprints[key].append(run)
    collapsed = []
    for group in fingerprints.values():
        if len(group) > 1:
            collapsed.append({"run_ids": [run.id for run in group], "labels": [_scenario_label(run) for run in group], "tolerance": tolerance})
    return collapsed


def run_causal_checks(patient: Patient | None = None) -> list[dict[str, Any]]:
    queryset = CounterfactualRun.objects.filter(status=CounterfactualRun.STATUS_COMPLETED)
    assumption_sets = CausalAssumptionSet.objects.all()
    if patient is not None:
        queryset = queryset.filter(patient=patient)
        assumption_sets = assumption_sets.filter(patient=patient)
    missing_classification = [run.id for run in queryset if not (run.simulation_summary or {}).get("classification")]
    causal_claims = [run.id for run in queryset if (run.simulation_summary or {}).get("causal_effect_estimated") is True]
    checks = [
        build_check("pass" if not missing_classification else "warn", "counterfactual causal classification", "Every completed run should state mechanistic-vs-causal status.", object_ids=missing_classification, next_action="Re-run affected counterfactuals."),
        build_check("pass" if assumption_sets.exists() else "warn", "causal assumption set exists", "Research cases should document assumptions and identification status.", object_ids=list(assumption_sets.values_list("id", flat=True)), next_action="Create or review CausalAssumptionSet."),
        build_check("pass" if not causal_claims else "fail", "no causal effect overclaim", "Runs must not claim causal effects without an external causal design.", object_ids=causal_claims, next_action="Remove unsupported causal-effect flags."),
    ]
    return checks


def run_scientific_reference_checks() -> list[dict[str, Any]]:
    references = load_model_references()
    checks = []
    for reference in references:
        status = "pass" if reference.get("status") == "implemented" else "warn"
        checks.append(build_check(status, f"reference: {reference.get('component')}", reference.get("title", ""), next_action="Attach reviewed citation before claiming external validity.", raw_details=reference))
    return checks


def run_privacy_security_checks() -> list[dict[str, Any]]:
    checks = []
    repo_root = Path(settings.BASE_DIR)
    ignored = _git_check_ignore(repo_root, ["local_private", "media", "db.sqlite3"])
    checks.append(build_check("pass" if len(ignored) >= 3 else "fail", "ignored sensitive paths", "local_private, media, and db.sqlite3 should be ignored.", raw_details={"ignored": ignored}, next_action="Update .gitignore."))
    staged = _git_diff_name_only(repo_root, "--cached")
    unsafe_staged = [path for path in staged if is_sensitive_path(path)]
    checks.append(build_check("pass" if not unsafe_staged else "fail", "no sensitive staged paths", "Sensitive local paths must not be staged.", object_ids=unsafe_staged, next_action="Unstage local/private/media/db/PDF/image files."))
    text_findings = []
    for path in staged:
        full_path = repo_root / path
        if full_path.exists() and full_path.is_file() and full_path.stat().st_size < 500_000:
            try:
                text_findings.extend({**finding, "path": path} for finding in scan_text_for_sensitive_markers(full_path.read_text(encoding="utf-8", errors="ignore")))
            except OSError:
                continue
    checks.append(build_check("pass" if not text_findings else "fail", "staged text PHI scan", "Staged text should not contain obvious patient identifiers.", raw_details={"findings": text_findings}, next_action="Remove PHI from staged text."))
    return checks


def run_artifact_checks(patient: Patient | None = None) -> list[dict[str, Any]]:
    queryset = CounterfactualRun.objects.filter(status=CounterfactualRun.STATUS_COMPLETED)
    metadata = SimulationRunMetadata.objects.all()
    if patient is not None:
        queryset = queryset.filter(patient=patient)
        metadata = metadata.filter(counterfactual_run__patient=patient)
    missing_report = [run.id for run in queryset if not run.report_artifact]
    missing_trajectory = [run.id for run in queryset if not run.trajectory_artifact]
    missing_metadata = [run.id for run in queryset if not run.metadata_records.exists()]
    missing_hash = [record.id for record in metadata if not record.output_hash]
    return [
        build_check("pass" if not missing_report else "warn", "research reports exist", "Completed runs should have report artifacts.", object_ids=missing_report, next_action="Regenerate reports."),
        build_check("pass" if not missing_trajectory else "warn", "trajectory artifacts exist", "Completed runs should have trajectory artifacts.", object_ids=missing_trajectory, next_action="Regenerate trajectories."),
        build_check("pass" if not missing_metadata else "warn", "provenance metadata exists", "Completed runs should have provenance metadata.", object_ids=missing_metadata, next_action="Re-run affected scenarios."),
        build_check("pass" if not missing_hash else "warn", "output hashes exist", "Metadata should include output hashes where implemented.", object_ids=missing_hash, next_action="Backfill or regenerate metadata."),
    ]


def load_model_references() -> list[dict[str, Any]]:
    path = Path(__file__).resolve().parent / "model_references.json"
    return json.loads(path.read_text(encoding="utf-8"))


def _latest_completed_runs(patient: Patient | None = None) -> dict[str, CounterfactualRun]:
    queryset = CounterfactualRun.objects.filter(status=CounterfactualRun.STATUS_COMPLETED).order_by("-id")
    if patient is not None:
        queryset = queryset.filter(patient=patient)
    latest: dict[str, CounterfactualRun] = {}
    for run in queryset:
        latest.setdefault(_scenario_label(run), run)
    return latest


def _scenario_label(run: CounterfactualRun) -> str:
    definition = run.intervention_definition or {}
    return str(definition.get("label") or definition.get("intervention", {}).get("label") or f"run_{run.id}")


def _load_trajectory(run: CounterfactualRun) -> dict[str, Any] | None:
    if not run.trajectory_artifact:
        return None
    media_url = settings.MEDIA_URL.rstrip("/")
    relative = run.trajectory_artifact[len(media_url) + 1 :] if run.trajectory_artifact.startswith(media_url + "/") else run.trajectory_artifact.lstrip("/")
    path = Path(settings.MEDIA_ROOT) / relative
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None


def _git_diff_name_only(repo_root: Path, *args: str) -> list[str]:
    try:
        completed = subprocess.run(["git", "diff", *args, "--name-only"], cwd=repo_root, capture_output=True, text=True, check=False)
    except Exception:
        return []
    return [line.strip() for line in completed.stdout.splitlines() if line.strip()]


def _git_check_ignore(repo_root: Path, paths: list[str]) -> list[str]:
    ignored = []
    for path in paths:
        try:
            completed = subprocess.run(["git", "check-ignore", path], cwd=repo_root, capture_output=True, text=True, check=False)
        except Exception:
            continue
        if completed.returncode == 0:
            ignored.append(path)
    return ignored