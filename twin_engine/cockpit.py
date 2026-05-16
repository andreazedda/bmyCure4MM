from __future__ import annotations

import json
from collections import defaultdict
from pathlib import Path
from typing import Any

from django.conf import settings
from django.urls import reverse

from clinic.models import Assessment, Patient, PatientTherapy

from .developer_checks import detect_schedule_collapse, load_model_references, run_developer_checks
from .models import (
    AdverseEvent,
    CausalAssumptionSet,
    CounterfactualRun,
    LongitudinalLabResult,
    ObservationResidual,
    SimulationRunMetadata,
    TherapyInterruption,
)
from .state_model import get_current_twin_state
from .toxicity_model import compute_toxicity_constraints


ASSESSMENT_FIELDS = [
    ("m_protein_g_dl", "M-protein"),
    ("flc_ratio", "FLC ratio"),
    ("hemoglobin_g_dl", "Hemoglobin"),
    ("creatinine_mg_dl", "Creatinine"),
    ("beta2m_mg_l", "Beta-2 microglobulin"),
    ("ldH_u_l", "LDH"),
    ("albumin_g_dl", "Albumin"),
]
INITIALIZATION_KEY_FIELDS = ("m_protein_g_dl", "flc_ratio", "hemoglobin_g_dl")

ANALYTE_GROUPS = {
    "disease": {
        "title": "Disease markers",
        "analytes": [
            LongitudinalLabResult.ANALYTE_M_PROTEIN,
            LongitudinalLabResult.ANALYTE_FLC_RATIO,
            LongitudinalLabResult.ANALYTE_KAPPA_FLC,
            LongitudinalLabResult.ANALYTE_LAMBDA_FLC,
        ],
    },
    "hematology": {
        "title": "Hematology",
        "analytes": [
            LongitudinalLabResult.ANALYTE_HB,
            LongitudinalLabResult.ANALYTE_WBC,
            LongitudinalLabResult.ANALYTE_NEU,
            LongitudinalLabResult.ANALYTE_PLT,
        ],
    },
    "toxicity": {
        "title": "Toxicity and chemistry",
        "analytes": [
            LongitudinalLabResult.ANALYTE_AST,
            LongitudinalLabResult.ANALYTE_ALT,
            LongitudinalLabResult.ANALYTE_CREATININE,
            LongitudinalLabResult.ANALYTE_LDH,
            LongitudinalLabResult.ANALYTE_BETA2M,
        ],
    },
}

ASSESSMENT_TO_ANALYTE = {
    "m_protein_g_dl": LongitudinalLabResult.ANALYTE_M_PROTEIN,
    "flc_ratio": LongitudinalLabResult.ANALYTE_FLC_RATIO,
    "hemoglobin_g_dl": LongitudinalLabResult.ANALYTE_HB,
    "creatinine_mg_dl": LongitudinalLabResult.ANALYTE_CREATININE,
    "beta2m_mg_l": LongitudinalLabResult.ANALYTE_BETA2M,
    "ldH_u_l": LongitudinalLabResult.ANALYTE_LDH,
    "albumin_g_dl": LongitudinalLabResult.ANALYTE_ALBUMIN,
}

ANALYTE_LABELS = dict(LongitudinalLabResult.ANALYTE_CHOICES)


def build_research_cockpit_context(patient: Patient, *, include_developer_checks: bool = False) -> dict[str, Any]:
    current_state = get_current_twin_state(patient)
    recommendation = build_assessment_recommendations(patient)
    latest_runs_by_label = list_latest_completed_runs_by_label(patient)
    collapse_warnings = detect_schedule_collapse(patient)
    scenario_rows = build_scenario_rows(patient, latest_runs_by_label, collapse_warnings)
    toxicity_constraints = compute_toxicity_constraints(patient)
    metadata_records = _metadata_records_for_patient(patient, current_state)
    causal_sets = list(patient.causal_assumption_sets.order_by("-created_at")[:10])
    context = {
        "patient": patient,
        "current_twin_state": current_state,
        "data_availability": build_data_availability(patient),
        "lab_chart_data": build_lab_chart_data(patient),
        "event_overlay_data": build_event_overlay_data(patient),
        "assessment_recommendation": recommendation,
        "assessment_options": recommendation["assessment_rows"],
        "calibration_panel": build_calibration_panel(patient, current_state),
        "scenario_rows": scenario_rows,
        "trajectory_chart_data": build_trajectory_chart_data(scenario_rows),
        "toxicity_panel": build_toxicity_panel(patient, toxicity_constraints),
        "causal_panel": build_causal_panel(causal_sets),
        "scientific_references": load_model_references(),
        "provenance_records": metadata_records,
        "developer_check_summary": summarize_checks(run_developer_checks(patient)) if include_developer_checks else None,
        "next_actions": build_next_actions(patient, current_state, recommendation, scenario_rows),
        "raw_developer_payload": {
            "patient_id": patient.id,
            "current_twin_state_id": getattr(current_state, "id", None),
            "run_ids": [row["run"].id for row in scenario_rows],
        },
    }
    return context


def build_assessment_recommendations(patient: Patient) -> dict[str, Any]:
    rows = []
    for assessment in patient.assessments.order_by("date"):
        present = [field for field, _label in ASSESSMENT_FIELDS if _has_value(getattr(assessment, field, None))]
        missing = [label for field, label in ASSESSMENT_FIELDS if field not in present]
        key_present = [field for field in INITIALIZATION_KEY_FIELDS if field in present]
        rows.append(
            {
                "assessment": assessment,
                "date": assessment.date.isoformat(),
                "score": len(present),
                "key_score": len(key_present),
                "available_labels": [label for field, label in ASSESSMENT_FIELDS if field in present],
                "missing_labels": missing,
                "is_usable": len(key_present) >= 2,
                "reason": _assessment_reason(present, missing),
            }
        )
    recommended = None
    usable = [row for row in rows if row["is_usable"]]
    if usable:
        recommended = sorted(
            usable,
            key=lambda row: (-row["key_score"], -row["score"], row["assessment"].date, row["assessment"].id),
        )[0]
    return {
        "recommended": recommended,
        "assessment_rows": rows,
        "has_usable_assessment": recommended is not None,
    }


def build_data_availability(patient: Patient) -> list[dict[str, Any]]:
    analyte_counts = LongitudinalLabResult.objects.filter(patient=patient).values("analyte").order_by("analyte")
    analyte_summary = {item["analyte"]: LongitudinalLabResult.objects.filter(patient=patient, analyte=item["analyte"]).count() for item in analyte_counts.distinct()}
    items = [
        _availability_item("Assessments", patient.assessments.count(), "Clinical assessment rows available for initialization and residual plots."),
        _availability_item("Structured labs", patient.longitudinal_lab_results.count(), "LongitudinalLabResult rows used for cockpit charts."),
        _availability_item("Therapies", patient.therapies.count(), "PatientTherapy rows with regimen and schedule context."),
        _availability_item("Therapy interruptions", patient.therapy_interruptions.count(), "Interruption rows used in toxicity and event overlays."),
        _availability_item("Adverse events", patient.adverse_events.count(), "Structured toxicity/adverse event rows."),
        _availability_item("Twin states", patient.twin_states.count(), "Persisted research twin states."),
        _availability_item("Residual records", patient.observation_residuals.count(), "Observed-vs-predicted calibration residuals."),
        _availability_item("Completed what-if runs", patient.counterfactual_runs.filter(status=CounterfactualRun.STATUS_COMPLETED).count(), "Completed mechanistic counterfactual simulations."),
    ]
    items.append(
        {
            "label": "Analytes present",
            "count": len(analyte_summary),
            "status": "available" if analyte_summary else "missing",
            "detail": ", ".join(sorted(analyte_summary)) or "No structured analytes available.",
            "raw": analyte_summary,
        }
    )
    return items


def build_lab_chart_data(patient: Patient) -> dict[str, Any]:
    lab_points = defaultdict(list)
    units = {}
    for result in patient.longitudinal_lab_results.order_by("date", "analyte", "id"):
        if result.value is None:
            continue
        lab_points[result.analyte].append(
            {
                "x": result.date.isoformat(),
                "y": float(result.value),
                "source": result.get_source_quality_display(),
                "id": result.id,
            }
        )
        if result.unit:
            units[result.analyte] = result.unit

    for assessment in patient.assessments.order_by("date"):
        for field, analyte in ASSESSMENT_TO_ANALYTE.items():
            value = getattr(assessment, field, None)
            if not _has_value(value):
                continue
            if not any(point["x"] == assessment.date.isoformat() for point in lab_points[analyte]):
                lab_points[analyte].append(
                    {
                        "x": assessment.date.isoformat(),
                        "y": float(value),
                        "source": "Assessment fallback",
                        "id": f"assessment-{assessment.id}",
                    }
                )

    groups = {}
    for group_key, group in ANALYTE_GROUPS.items():
        series = []
        missing = []
        for analyte in group["analytes"]:
            points = sorted(lab_points.get(analyte, []), key=lambda item: item["x"])
            if points:
                series.append({"analyte": analyte, "label": ANALYTE_LABELS.get(analyte, analyte), "unit": units.get(analyte, ""), "points": points})
            else:
                missing.append(ANALYTE_LABELS.get(analyte, analyte))
        groups[group_key] = {"title": group["title"], "series": series, "missing": missing}
    return groups


def build_event_overlay_data(patient: Patient) -> dict[str, Any]:
    therapies = [
        {
            "id": therapy.id,
            "label": therapy.regimen.name,
            "start": therapy.start_date.isoformat(),
            "end": therapy.end_date.isoformat() if therapy.end_date else None,
            "source_quality": therapy.get_source_quality_display(),
        }
        for therapy in patient.therapies.select_related("regimen").order_by("start_date")
    ]
    interruptions = [
        {
            "id": interruption.id,
            "label": interruption.get_reason_display(),
            "drug": interruption.drug,
            "start": interruption.start_date.isoformat(),
            "end": interruption.end_date.isoformat() if interruption.end_date else None,
        }
        for interruption in patient.therapy_interruptions.order_by("start_date")
    ]
    adverse_events = [
        {
            "id": event.id,
            "label": event.get_event_type_display(),
            "date": event.date.isoformat(),
            "grade": event.grade,
            "suspected_drug": event.suspected_drug,
        }
        for event in patient.adverse_events.order_by("date")
    ]
    return {"therapies": therapies, "interruptions": interruptions, "adverse_events": adverse_events}


def build_calibration_panel(patient: Patient, current_state) -> dict[str, Any]:
    if current_state is None:
        return {
            "status": "not_initialized",
            "headline": "No current twin state",
            "detail": "Initialize a research twin before calibration quality can be assessed.",
            "residual_rows": [],
            "diagnostics": {},
        }
    diagnostics = current_state.parameter_uncertainty or {}
    pre_rows = ObservationResidual.objects.filter(twin_state=current_state, stage=ObservationResidual.STAGE_PRE_CALIBRATION).select_related("assessment").order_by("assessment__date", "created_at")
    post_rows = ObservationResidual.objects.filter(twin_state=current_state, stage=ObservationResidual.STAGE_POST_CALIBRATION).select_related("assessment").order_by("assessment__date", "created_at")
    residual_rows = []
    for residual in list(pre_rows) + list(post_rows):
        residual_rows.append(
            {
                "id": residual.id,
                "stage": residual.get_stage_display(),
                "date": residual.assessment.date.isoformat() if residual.assessment else "",
                "rmse": _round_or_none(residual.rmse),
                "mae": _round_or_none(residual.mae),
                "observed": residual.observed_values,
                "predicted": residual.predicted_values,
            }
        )
    rmse_before = diagnostics.get("rmse_before")
    rmse_after = diagnostics.get("rmse_after")
    improvement = None
    if rmse_before is not None and rmse_after is not None:
        improvement = float(rmse_before) - float(rmse_after)
    status = diagnostics.get("calibration_status") or current_state.get_method_display()
    headline = "Calibrated with residual diagnostics" if residual_rows else "Current state has no residual diagnostics"
    return {
        "status": status,
        "headline": headline,
        "detail": "Calibration is a model-fit diagnostic, not clinical validation.",
        "diagnostics": diagnostics,
        "rmse_before": _round_or_none(rmse_before),
        "rmse_after": _round_or_none(rmse_after),
        "rmse_delta": _round_or_none(improvement),
        "residual_rows": residual_rows,
        "residual_count": len(residual_rows),
    }


def list_latest_completed_runs_by_label(patient: Patient) -> dict[str, CounterfactualRun]:
    latest: dict[str, CounterfactualRun] = {}
    queryset = patient.counterfactual_runs.select_related("base_twin_state", "alternative_regimen").filter(status=CounterfactualRun.STATUS_COMPLETED).order_by("-id")
    for run in queryset:
        latest.setdefault(scenario_label(run), run)
    return latest


def build_scenario_rows(patient: Patient, latest_runs_by_label: dict[str, CounterfactualRun], collapse_warnings: list[dict[str, Any]]) -> list[dict[str, Any]]:
    collapsed_run_ids = {run_id for warning in collapse_warnings for run_id in warning.get("run_ids", [])}
    rows = []
    for label, run in latest_runs_by_label.items():
        metrics = run.comparison_metrics or {}
        summary = run.simulation_summary or {}
        predicted = summary.get("predicted_biomarkers") or {}
        baseline_predicted = summary.get("baseline_predicted_biomarkers") or {}
        trajectory_payload = load_json_artifact(run.trajectory_artifact)
        report_payload = load_json_artifact(run.report_artifact)
        issues = []
        if run.id in collapsed_run_ids:
            issues.append("Trajectory matches another scenario at current schedule resolution.")
        if not trajectory_payload:
            issues.append("Trajectory artifact unavailable.")
        if not predicted:
            issues.append("Predicted biomarker summary unavailable.")
        rows.append(
            {
                "run": run,
                "label": label,
                "utility": _round_or_none(metrics.get("research_utility")),
                "toxicity_penalty": _round_or_none(metrics.get("toxicity_constraint_penalty")),
                "tumor_reduction_delta": _metric_delta(metrics, "tumor_reduction"),
                "healthy_loss_delta": _metric_delta(metrics, "healthy_loss"),
                "durability_delta": _metric_delta(metrics, "durability_index"),
                "predicted_biomarkers": predicted,
                "predicted_biomarker_summary": _format_predicted_biomarkers(predicted),
                "baseline_predicted_biomarkers": baseline_predicted,
                "classification": summary.get("classification") or {},
                "issues": issues,
                "interpretation": interpret_scenario(metrics, predicted, issues),
                "trajectory_payload": trajectory_payload,
                "report_payload": report_payload,
                "report_url": reverse("twin_engine:counterfactual_report", args=[patient.id, run.id]),
            }
        )
    rows.sort(key=lambda item: (item["utility"] is None, -(item["utility"] or -999999), -item["run"].id))
    return rows


def build_trajectory_chart_data(scenario_rows: list[dict[str, Any]]) -> dict[str, Any]:
    series = []
    baseline_added = False
    for row in scenario_rows[:6]:
        payload = row.get("trajectory_payload") or {}
        baseline = payload.get("baseline_trajectory") or {}
        alternative = payload.get("alternative_trajectory") or {}
        if baseline and not baseline_added:
            series.append(_trajectory_series("Recorded therapy baseline", baseline, "baseline"))
            baseline_added = True
        if alternative:
            series.append(_trajectory_series(row["label"], alternative, "alternative"))
    return {"series": [item for item in series if item.get("x")], "has_data": bool(series)}


def build_toxicity_panel(patient: Patient, toxicity_constraints: dict[str, Any]) -> dict[str, Any]:
    labs = build_lab_chart_data(patient).get("toxicity", {})
    return {
        "summary": toxicity_constraints,
        "lab_series": labs.get("series", []),
        "interruptions": list(patient.therapy_interruptions.order_by("-start_date")[:10]),
        "adverse_events": list(patient.adverse_events.order_by("-date")[:10]),
        "disclaimer": "Descriptive toxicity constraints only; no predictive AST/ALT/neutrophil dynamics are implemented.",
    }


def build_causal_panel(causal_sets: list[CausalAssumptionSet]) -> dict[str, Any]:
    if not causal_sets:
        return {
            "status": "missing",
            "headline": "No causal assumption set recorded",
            "detail": "What-if runs remain mechanistic simulations and are not causal effect estimates.",
            "sets": [],
        }
    return {
        "status": "mechanistic_only",
        "headline": "Mechanistic what-if status documented",
        "detail": "The current assumption set labels this workflow as mechanistic only unless a separate causal design is added.",
        "sets": causal_sets,
    }


def summarize_checks(check_groups: dict[str, list[dict[str, Any]]]) -> dict[str, Any]:
    summary = {"pass": 0, "warn": 0, "fail": 0}
    total = 0
    for checks in check_groups.values():
        for check in checks:
            status = check.get("status", "warn")
            summary[status] = summary.get(status, 0) + 1
            total += 1
    return {"groups": check_groups, "summary": summary, "total": total}


def build_next_actions(patient: Patient, current_state, recommendation: dict[str, Any], scenario_rows: list[dict[str, Any]]) -> list[dict[str, str]]:
    actions = []
    if current_state is None:
        actions.append({"label": "Initialize twin", "detail": "Create the first research twin state from the recommended assessment.", "href": reverse("twin_engine:initialize_twin_from_assessment", args=[patient.id])})
    else:
        actions.append({"label": "Review calibration", "detail": "Check residuals and recalibrate after new longitudinal observations.", "href": "#calibration-quality"})
    if not scenario_rows:
        actions.append({"label": "Run what-if scenarios", "detail": "Generate completed counterfactual runs before comparing alternatives.", "href": "#what-if-scenarios"})
    else:
        actions.append({"label": "Compare scenario trajectories", "detail": "Use the trajectory panel to inspect model behavior and schedule-resolution limits.", "href": "#trajectory-comparison"})
    if not recommendation.get("has_usable_assessment"):
        actions.append({"label": "Add modeled markers", "detail": "M-protein, FLC ratio, and hemoglobin improve initialization and calibration interpretability.", "href": "#data-availability"})
    actions.append({"label": "Open developer console", "detail": "Run internal checks before sharing artifacts or pushing code.", "href": reverse("twin_engine:developer_console")})
    return actions


def write_local_feedback(user, payload: dict[str, Any]) -> Path:
    path = Path(settings.BASE_DIR) / "local_private" / "feedback" / "research_cockpit_feedback.jsonl"
    path.parent.mkdir(parents=True, exist_ok=True)
    record = {
        "user_id": getattr(user, "id", None),
        "payload": payload,
    }
    with path.open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(record, sort_keys=True, default=str) + "\n")
    return path


def load_json_artifact(artifact_url: str) -> dict[str, Any] | None:
    if not artifact_url:
        return None
    media_url = settings.MEDIA_URL.rstrip("/")
    relative = artifact_url[len(media_url) + 1 :] if artifact_url.startswith(media_url + "/") else artifact_url.lstrip("/")
    path = Path(settings.MEDIA_ROOT) / relative
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except json.JSONDecodeError:
        return None


def scenario_label(run: CounterfactualRun) -> str:
    definition = run.intervention_definition or {}
    return str(definition.get("label") or definition.get("intervention", {}).get("label") or f"Run {run.id}")


def interpret_scenario(metrics: dict[str, Any], predicted: dict[str, Any], issues: list[str]) -> str:
    utility = metrics.get("research_utility")
    penalty = metrics.get("toxicity_constraint_penalty")
    fragments = ["Research simulation only; heuristic ranking is for model exploration."]
    if utility is not None:
        fragments.append(f"Utility score {float(utility):.3f} combines modeled disease control, healthy-cell preservation, durability, and toxicity penalty.")
    if penalty:
        fragments.append(f"Toxicity history subtracts {float(penalty):.3f} from the exploratory utility score.")
    if predicted:
        fragments.append("Predicted biomarker endpoints are available for inspection.")
    if issues:
        fragments.append("Model-execution caveat: " + " ".join(issues))
    return " ".join(fragments)


def _availability_item(label: str, count: int, detail: str) -> dict[str, Any]:
    return {"label": label, "count": count, "status": "available" if count else "missing", "detail": detail}


def _assessment_reason(present: list[str], missing: list[str]) -> str:
    if all(field in present for field in INITIALIZATION_KEY_FIELDS):
        return "Contains M-protein, FLC ratio, and hemoglobin for initialization."
    if len([field for field in INITIALIZATION_KEY_FIELDS if field in present]) >= 2:
        return "Contains at least two modeled initialization markers."
    return "Missing multiple modeled initialization markers: " + ", ".join(missing[:4])


def _metadata_records_for_patient(patient: Patient, current_state) -> list[SimulationRunMetadata]:
    queryset = SimulationRunMetadata.objects.filter(counterfactual_run__patient=patient)
    if current_state is not None:
        queryset = queryset | SimulationRunMetadata.objects.filter(twin_state=current_state)
    return list(queryset.distinct().order_by("-created_at")[:12])


def _trajectory_series(label: str, trajectory: dict[str, Any], kind: str) -> dict[str, Any]:
    time = trajectory.get("time") or trajectory.get("days") or []
    tumor = trajectory.get("tumor_cells") or []
    healthy = trajectory.get("healthy_cells") or []
    return {
        "label": label,
        "kind": kind,
        "x": [float(value) for value in time[: len(tumor)]],
        "tumor": [float(value) for value in tumor],
        "healthy": [float(value) for value in healthy],
    }


def _metric_delta(metrics: dict[str, Any], key: str) -> float | None:
    value = (metrics or {}).get(key) or {}
    return _round_or_none(value.get("delta"))


def _format_predicted_biomarkers(predicted: dict[str, Any]) -> str:
    if not predicted:
        return "Unavailable"
    labels = {
        "m_protein_g_dl": "M-protein",
        "flc_ratio": "FLC ratio",
        "hemoglobin_g_dl": "Hemoglobin",
    }
    parts = []
    for key, label in labels.items():
        if predicted.get(key) is not None:
            parts.append(f"{label}: {_round_or_none(predicted.get(key))}")
    return "; ".join(parts) if parts else "Available in artifact"


def _round_or_none(value) -> float | None:
    if value is None:
        return None
    try:
        return round(float(value), 4)
    except (TypeError, ValueError):
        return None


def _has_value(value) -> bool:
    return value is not None and value != ""