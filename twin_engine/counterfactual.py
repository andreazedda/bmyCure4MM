from __future__ import annotations

import logging
from datetime import timedelta
from typing import Any

from django.db.models import Q

from clinic.models import Regimen

from .causal import distinguish_mechanistic_counterfactual_vs_causal_estimand
from .models import CounterfactualRun
from .provenance import CURRENT_MODEL_VERSION, record_simulation_metadata
from .report_builder import build_counterfactual_report_payload, write_json_artifact
from .simulation_bridge import run_patient_simulation
from .therapy_schedule import SUPPORTED_DRUG_ALIASES, build_therapy_schedule
from .toxicity_model import compute_toxicity_constraints
from .validators import validate_research_run_inputs


logger = logging.getLogger("twin_engine.research")


def run_counterfactual(patient, base_twin_state, intervention_definition, horizon_days, user=None):
    execution_definition = _normalize_intervention_definition(
        intervention_definition,
        base_twin_state.state_date,
        horizon_days,
    )
    validate_research_run_inputs(patient, base_twin_state, execution_definition, horizon_days)

    alternative_regimen = None
    alternative_regimen_id = execution_definition.get("alternative_regimen_id")
    if alternative_regimen_id:
        alternative_regimen = Regimen.objects.filter(pk=alternative_regimen_id).first()

    actual_therapy = _get_active_therapy(patient, base_twin_state.state_date)

    counterfactual_run = CounterfactualRun.objects.create(
        patient=patient,
        base_twin_state=base_twin_state,
        actual_therapy=actual_therapy,
        alternative_regimen=alternative_regimen,
        alternative_parameters=intervention_definition.get("parameter_overrides", {}) or {},
        intervention_definition=intervention_definition,
        status=CounterfactualRun.STATUS_RUNNING,
        created_by=user,
    )

    try:
        end_date = base_twin_state.state_date + timedelta(days=int(horizon_days))
        baseline_schedule = build_therapy_schedule(patient, base_twin_state.state_date, end_date)
        baseline_result = run_patient_simulation(
            base_twin_state,
            therapy_schedule=baseline_schedule,
            horizon_days=int(horizon_days),
        )

        alternative_overrides = dict(execution_definition.get("parameter_overrides", {}) or {})
        if execution_definition.get("drug_doses"):
            alternative_overrides["drug_doses"] = execution_definition["drug_doses"]
        if execution_definition.get("random_seed") is not None:
            alternative_overrides["seed"] = int(execution_definition["random_seed"])

        alternative_result = run_patient_simulation(
            base_twin_state,
            therapy_schedule=execution_definition.get("schedule"),
            horizon_days=int(horizon_days),
            overrides=alternative_overrides,
        )

        toxicity_constraints = compute_toxicity_constraints(patient)
        comparison_metrics = compare_counterfactual_to_baseline(
            baseline_result["summary"],
            alternative_result["summary"],
            baseline_solver_inputs=baseline_result["solver_inputs"]["raw_parameters"],
            alternative_solver_inputs=alternative_result["solver_inputs"]["raw_parameters"],
            toxicity_constraints=toxicity_constraints,
        )
        classification = distinguish_mechanistic_counterfactual_vs_causal_estimand(
            graph_definition=execution_definition.get("causal_graph"),
            intervention=intervention_definition,
            outcome={"metrics": ["tumor_reduction", "healthy_loss", "time_to_recurrence"]},
            adjustment_set=execution_definition.get("adjustment_set"),
            identification_status=execution_definition.get("identification_status"),
        )
        warning_block = [
            "Research simulation only.",
            classification.get("warning"),
            toxicity_constraints.get("reason"),
        ]

        counterfactual_run.simulation_summary = {
            "label": "research simulation",
            "counterfactual_type": "mechanistic model counterfactual",
            "causal_interpretation": "unvalidated exploratory hypothesis",
            "classification": classification,
            "baseline": baseline_result["summary"],
            "alternative": alternative_result["summary"],
            "baseline_predicted_biomarkers": baseline_result["summary"].get("predicted_biomarkers"),
            "predicted_biomarkers": alternative_result["summary"].get("predicted_biomarkers"),
            "toxicity_constraints": toxicity_constraints,
            "warning_block": warning_block,
        }
        counterfactual_run.comparison_metrics = comparison_metrics

        trajectory_payload = {
            "label": "research simulation",
            "counterfactual_type": "mechanistic model counterfactual",
            "base_twin_state_id": base_twin_state.id,
            "baseline_trajectory": baseline_result["trajectory"],
            "alternative_trajectory": alternative_result["trajectory"],
            "comparison_metrics": comparison_metrics,
        }
        trajectory_url, _ = write_json_artifact(
            "counterfactual_trajectory",
            trajectory_payload,
            patient_id=patient.id,
            run_id=counterfactual_run.id,
            folder_name="research_trajectories",
        )

        report_metadata = {
            "model_version": CURRENT_MODEL_VERSION,
            "solver_name": "MathematicalModel",
            "input_hash_source": {
                "base_twin_state_id": base_twin_state.id,
                "intervention_definition": _sanitize_payload_for_artifact(intervention_definition),
            },
        }
        report_payload = build_counterfactual_report_payload(
            counterfactual_run,
            baseline_summary=baseline_result["summary"],
            alternative_summary=alternative_result["summary"],
            comparison_metrics=comparison_metrics,
            metadata=report_metadata,
            toxicity_constraints=toxicity_constraints,
            warnings=warning_block,
        )
        report_url, _ = write_json_artifact(
            "counterfactual_report",
            report_payload,
            patient_id=patient.id,
            run_id=counterfactual_run.id,
            folder_name="research_reports",
        )

        counterfactual_run.trajectory_artifact = trajectory_url
        counterfactual_run.report_artifact = report_url
        counterfactual_run.status = CounterfactualRun.STATUS_COMPLETED
        counterfactual_run.error_message = ""
        counterfactual_run.save(
            update_fields=[
                "alternative_regimen",
                "alternative_parameters",
                "simulation_summary",
                "comparison_metrics",
                "trajectory_artifact",
                "report_artifact",
                "status",
                "error_message",
            ]
        )

        record_simulation_metadata(
            counterfactual_run=counterfactual_run,
            twin_state=base_twin_state,
            model_version=CURRENT_MODEL_VERSION,
            solver_name="MathematicalModel",
            input_payload={
                "base_twin_state_id": base_twin_state.id,
                "intervention_definition": _sanitize_payload_for_artifact(intervention_definition),
                "execution_definition": _sanitize_payload_for_artifact(execution_definition),
            },
            solver_parameters=alternative_result["solver_inputs"]["raw_parameters"],
            output_payload=report_payload,
            random_seed=execution_definition.get("random_seed"),
        )

        logger.info(
            "counterfactual_completed patient_id=%s state_id=%s run_id=%s",
            patient.id,
            base_twin_state.id,
            counterfactual_run.id,
        )
        return counterfactual_run
    except Exception as exc:
        counterfactual_run.status = CounterfactualRun.STATUS_FAILED
        counterfactual_run.error_message = str(exc)
        counterfactual_run.save(update_fields=["status", "error_message"])
        logger.exception(
            "counterfactual_failed patient_id=%s state_id=%s run_id=%s",
            patient.id,
            base_twin_state.id,
            counterfactual_run.id,
        )
        raise


def compare_counterfactual_to_baseline(
    baseline_summary,
    alternative_summary,
    *,
    baseline_solver_inputs: dict[str, Any] | None = None,
    alternative_solver_inputs: dict[str, Any] | None = None,
    toxicity_constraints: dict[str, Any] | None = None,
) -> dict[str, Any]:
    comparison = {}
    for key in ("tumor_reduction", "healthy_loss", "time_to_recurrence", "durability_index"):
        baseline_value = baseline_summary.get(key)
        alternative_value = alternative_summary.get(key)
        if baseline_value is None or alternative_value is None:
            comparison[key] = {
                "baseline": baseline_value,
                "alternative": alternative_value,
                "delta": None,
            }
            continue
        comparison[key] = {
            "baseline": float(baseline_value),
            "alternative": float(alternative_value),
            "delta": float(alternative_value) - float(baseline_value),
        }

    biomarker_keys = ("m_protein_g_dl", "flc_ratio", "hemoglobin_g_dl")
    baseline_predicted = dict(baseline_summary.get("predicted_biomarkers") or {})
    alternative_predicted = dict(alternative_summary.get("predicted_biomarkers") or {})
    comparison["predicted_biomarkers"] = {}
    for key in biomarker_keys:
        baseline_value = baseline_predicted.get(key)
        alternative_value = alternative_predicted.get(key)
        if baseline_value is None or alternative_value is None:
            comparison["predicted_biomarkers"][key] = {
                "baseline": baseline_value,
                "alternative": alternative_value,
                "delta": None,
            }
            continue
        comparison["predicted_biomarkers"][key] = {
            "baseline": float(baseline_value),
            "alternative": float(alternative_value),
            "delta": float(alternative_value) - float(baseline_value),
        }

    penalty = _determine_toxicity_constraint_penalty(
        baseline_solver_inputs or {},
        alternative_solver_inputs or {},
        toxicity_constraints or {},
    )
    disease_control_score = float(alternative_summary.get("tumor_reduction") or 0.0)
    healthy_preservation_score = 1.0 - float(alternative_summary.get("healthy_loss") or 0.0)
    durability_score = float(alternative_summary.get("durability_index") or 0.0)
    comparison["research_utility"] = disease_control_score + healthy_preservation_score + durability_score - penalty
    comparison["toxicity_constraint_penalty"] = penalty
    comparison["utility_formula"] = (
        "disease_control_score + healthy_preservation_score + durability_score - "
        "heuristic_exposure_penalty_based_on_documented_toxicity_history"
    )
    comparison["utility_validity"] = "heuristic research ranking only; not clinical recommendation"
    return comparison


def _get_active_therapy(patient, state_date):
    return patient.therapies.filter(
        start_date__lte=state_date,
    ).filter(
        Q(end_date__isnull=True) | Q(end_date__gte=state_date)
    ).order_by("-start_date").first()


def _normalize_intervention_definition(intervention_definition, base_state_date, horizon_days: int) -> dict[str, Any]:
    if not isinstance(intervention_definition, dict):
        return intervention_definition

    if "intervention" not in intervention_definition:
        return dict(intervention_definition)

    intervention_payload = intervention_definition.get("intervention") or {}
    normalized = {
        "start_day": int(intervention_payload.get("start_day", 0) or 0),
        "duration_days": int(intervention_payload.get("duration_days", horizon_days) or horizon_days),
        "schedule": _build_schedule_payload(
            intervention_definition,
            intervention_payload,
            base_state_date,
            horizon_days,
        ),
        "parameter_overrides": dict(intervention_definition.get("parameter_overrides", {}) or {}),
        "random_seed": intervention_definition.get("random_seed"),
    }

    classification = intervention_definition.get("classification")
    if classification == "mechanistic_simulation_only":
        normalized["identification_status"] = "mechanistic_only"

    if intervention_definition.get("adjustment_set") is not None:
        normalized["adjustment_set"] = intervention_definition.get("adjustment_set")
    if intervention_definition.get("causal_graph") is not None:
        normalized["causal_graph"] = intervention_definition.get("causal_graph")
    if intervention_definition.get("alternative_regimen_id") is not None:
        normalized["alternative_regimen_id"] = intervention_definition.get("alternative_regimen_id")

    return normalized


def _build_schedule_payload(intervention_definition, intervention_payload, base_state_date, horizon_days: int) -> dict[str, Any]:
    phases = intervention_payload.get("phases") or []
    if phases:
        items = [
            _phase_to_schedule_item(
                intervention_definition=intervention_definition,
                intervention_payload=intervention_payload,
                phase=phase,
                base_state_date=base_state_date,
                default_label=intervention_definition.get("label") or "research_intervention",
            )
            for phase in phases
        ]
    else:
        items = [
            _phase_to_schedule_item(
                intervention_definition=intervention_definition,
                intervention_payload=intervention_payload,
                phase={
                    "start_day": intervention_payload.get("start_day", 0),
                    "duration_days": intervention_payload.get("duration_days", horizon_days),
                    "dose_mg": intervention_payload.get("dose_mg"),
                    "schedule": intervention_payload.get("schedule") or {"type": "daily"},
                    "drug": intervention_payload.get("drug"),
                },
                base_state_date=base_state_date,
                default_label=intervention_definition.get("label") or "research_intervention",
            )
        ]

    max_end_date = max(item["end_date"] for item in items)
    return {
        "patient_id": None,
        "start_date": min(item["start_date"] for item in items),
        "end_date": max_end_date,
        "items": items,
        "missing_doses": [],
        "validation": {"is_valid": True, "missing_doses": []},
    }


def _phase_to_schedule_item(*, intervention_definition, intervention_payload, phase, base_state_date, default_label: str) -> dict[str, Any]:
    drug_name = phase.get("drug") or intervention_payload.get("drug")
    canonical_drug = SUPPORTED_DRUG_ALIASES.get(str(drug_name or "").lower(), str(drug_name or "").lower())
    schedule_spec = dict(phase.get("schedule") or intervention_payload.get("schedule") or {"type": "daily"})
    start_day = int(phase.get("start_day", 0) or 0)
    duration_days = int(phase.get("duration_days", intervention_payload.get("duration_days", 1)) or 1)
    dose_mg = phase.get("dose_mg")

    cycle_length_days, days_on = _schedule_spec_to_cycle(schedule_spec)
    start_date = base_state_date + timedelta(days=start_day)
    end_date = start_date + timedelta(days=max(duration_days - 1, 0))

    return {
        "therapy_id": None,
        "regimen_id": None,
        "regimen_name": default_label,
        "start_date": start_date.isoformat(),
        "end_date": end_date.isoformat(),
        "cycle_length_days": cycle_length_days,
        "days_on": days_on,
        "components": [
            {
                "label": canonical_drug,
                "drug": canonical_drug,
                "supported_by_solver": True,
            }
        ],
        "doses": {
            canonical_drug: {
                "dose": float(dose_mg or 0.0),
                "unit": "mg",
                "schedule": schedule_spec,
            }
        },
        "source_quality": "curated_research",
        "provenance": {
            "source": "uploaded_clinical_documentation",
            "extraction_status": "manual_from_documents",
            "case_label": "MM_RESEARCH_CASE_001",
            "scenario_label": intervention_definition.get("label", default_label),
        },
    }


def _schedule_spec_to_cycle(schedule_spec: dict[str, Any]) -> tuple[int, list[int]]:
    schedule_type = str(schedule_spec.get("type") or "daily").lower()

    if schedule_type == "daily":
        return 1, [1]
    if schedule_type == "cycle":
        cycle_length_days = int(schedule_spec.get("cycle_length_days") or 28)
        raw_days_on = schedule_spec.get("days_on")
        if isinstance(raw_days_on, int):
            return cycle_length_days, list(range(1, raw_days_on + 1))
        if isinstance(raw_days_on, list):
            return cycle_length_days, [int(day) for day in raw_days_on]
        return cycle_length_days, list(range(1, cycle_length_days + 1))
    if schedule_type == "interval":
        every_days = max(int(schedule_spec.get("every_days") or 1), 1)
        return every_days, [1]
    if schedule_type == "hold":
        return 1, [1]

    raise ValueError(f"Unsupported intervention schedule type: {schedule_type}")


def _sanitize_payload_for_artifact(payload: Any) -> Any:
    banned_keys = {"notes", "note", "schedule_notes", "first_name", "last_name", "mrn", "birth_date"}

    if isinstance(payload, dict):
        sanitized = {}
        for key, value in payload.items():
            if str(key).lower() in banned_keys:
                continue
            sanitized[key] = _sanitize_payload_for_artifact(value)
        return sanitized
    if isinstance(payload, list):
        return [_sanitize_payload_for_artifact(item) for item in payload]
    return payload


def _determine_toxicity_constraint_penalty(
    baseline_solver_inputs: dict[str, Any],
    alternative_solver_inputs: dict[str, Any],
    toxicity_constraints: dict[str, Any],
) -> float:
    if not toxicity_constraints.get("lenalidomide_toxicity_limited"):
        return 0.0

    baseline_exposure = float(baseline_solver_inputs.get("lenalidomide_dose") or 0.0)
    alternative_exposure = float(alternative_solver_inputs.get("lenalidomide_dose") or 0.0)
    tolerance = max(abs(baseline_exposure), 1.0e-9) * 0.05
    if alternative_exposure > baseline_exposure + tolerance:
        return 1.0
    if abs(alternative_exposure - baseline_exposure) <= tolerance:
        return 0.5
    return 0.0
