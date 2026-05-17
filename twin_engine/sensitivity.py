from __future__ import annotations

import copy
from typing import Any

from .counterfactual import _normalize_intervention_definition
from .simulation_bridge import build_solver_inputs_from_twin_state
from .therapy_schedule import build_therapy_schedule
from .toxicity_model import compute_toxicity_constraints
from .uncertainty import _flatten_numeric_parameters, _restore_flat_parameters, _run_uncertainty_sample


COMPLETED_STATUS = "completed"
FAILED_STATUS = "failed"


def run_counterfactual_sensitivity(
    patient,
    twin_state,
    intervention_definition,
    horizon_days,
    *,
    perturbation_fractions: tuple[float, ...] = (0.10, 0.20),
) -> dict[str, Any]:
    horizon = int(horizon_days)
    execution_definition = _normalize_intervention_definition(intervention_definition, twin_state.state_date, horizon)
    end_date = twin_state.state_date + __import__("datetime").timedelta(days=max(horizon - 1, 0))
    baseline_schedule = build_therapy_schedule(patient, twin_state.state_date, end_date)
    baseline_solver_inputs = build_solver_inputs_from_twin_state(twin_state, baseline_schedule, horizon)
    toxicity_constraints = compute_toxicity_constraints(patient)

    baseline_metrics = _run_uncertainty_sample(
        patient=patient,
        twin_state=twin_state,
        execution_definition=execution_definition,
        horizon_days=horizon,
        baseline_solver_inputs=baseline_solver_inputs,
        toxicity_constraints=toxicity_constraints,
        sampled_parameters=copy.deepcopy(dict(twin_state.parameters or {})),
        sampled_toxicity_coefficients=None,
    )

    drivers = [
        *_model_parameter_drivers(twin_state),
        *_schedule_knob_drivers(intervention_definition),
    ]
    rows = []
    for driver in drivers:
        perturbations = []
        for fraction in perturbation_fractions:
            for direction in (-1.0, 1.0):
                perturbed_value = _perturbed_value(driver["baseline_value"], fraction, direction, driver["kind"])
                perturbed_metrics = _run_driver_perturbation(
                    patient=patient,
                    twin_state=twin_state,
                    intervention_definition=intervention_definition,
                    execution_definition=execution_definition,
                    horizon_days=horizon,
                    baseline_solver_inputs=baseline_solver_inputs,
                    toxicity_constraints=toxicity_constraints,
                    driver=driver,
                    perturbed_value=perturbed_value,
                )
                delta_utility_v2 = None
                if perturbed_metrics.get("research_utility_v2") is not None and baseline_metrics.get("research_utility_v2") is not None:
                    delta_utility_v2 = float(perturbed_metrics["research_utility_v2"]) - float(baseline_metrics["research_utility_v2"])
                perturbations.append(
                    {
                        "fraction": float(fraction),
                        "direction": "increase" if direction > 0 else "decrease",
                        "perturbed_value": perturbed_value,
                        "research_utility_v2": perturbed_metrics.get("research_utility_v2"),
                        "delta_utility_v2": delta_utility_v2,
                        "tumor_reduction": perturbed_metrics.get("tumor_reduction"),
                        "healthy_loss": perturbed_metrics.get("healthy_loss"),
                        "durability_index": perturbed_metrics.get("durability_index"),
                    }
                )

        deltas = [abs(item["delta_utility_v2"]) for item in perturbations if item["delta_utility_v2"] is not None]
        positive_change = next((item for item in perturbations if item["direction"] == "increase"), None)
        elasticity = None
        if positive_change and positive_change.get("delta_utility_v2") is not None and baseline_metrics.get("research_utility_v2") not in {None, 0.0}:
            elasticity = (
                float(positive_change["delta_utility_v2"]) / float(baseline_metrics["research_utility_v2"])
            ) / max(float(positive_change["fraction"]), 1.0e-9)
        rows.append(
            {
                "parameter": driver["name"],
                "kind": driver["kind"],
                "baseline_value": driver["baseline_value"],
                "baseline_utility_v2": baseline_metrics.get("research_utility_v2"),
                "max_abs_utility_v2_delta": max(deltas) if deltas else None,
                "elasticity": elasticity,
                "direction": _driver_direction(positive_change),
                "perturbations": perturbations,
                "sensitivity_classification": _classify_sensitivity(max(deltas) if deltas else None),
            }
        )

    rows.sort(key=lambda item: (item["max_abs_utility_v2_delta"] is None, -(item["max_abs_utility_v2_delta"] or 0.0), item["parameter"]))
    for rank, row in enumerate(rows, start=1):
        row["rank"] = rank

    return {
        "status": COMPLETED_STATUS,
        "baseline_metrics": baseline_metrics,
        "rows": rows,
        "top_drivers": rows[:5],
        "unstable_parameters": [row["parameter"] for row in rows if row["sensitivity_classification"] == "high"],
        "limitations": [
            "Sensitivity analysis perturbs one parameter or schedule knob at a time around the current mechanistic configuration.",
            "Large sensitivity indicates fragile exploratory ranking, not clinical effect modification.",
        ],
    }


def _run_driver_perturbation(
    *,
    patient,
    twin_state,
    intervention_definition,
    execution_definition,
    horizon_days: int,
    baseline_solver_inputs: dict[str, Any],
    toxicity_constraints: dict[str, Any],
    driver: dict[str, Any],
    perturbed_value,
) -> dict[str, Any]:
    parameters = copy.deepcopy(dict(twin_state.parameters or {}))
    effective_execution_definition = copy.deepcopy(execution_definition)

    if driver["kind"] == "model_parameter":
        flat = _flatten_numeric_parameters(parameters)
        flat[driver["path"]] = float(perturbed_value)
        parameters = _restore_flat_parameters(flat)
    else:
        modified_intervention = copy.deepcopy(intervention_definition)
        _set_path_value(modified_intervention, driver["path"], perturbed_value)
        effective_execution_definition = _normalize_intervention_definition(modified_intervention, twin_state.state_date, horizon_days)

    return _run_uncertainty_sample(
        patient=patient,
        twin_state=twin_state,
        execution_definition=effective_execution_definition,
        horizon_days=horizon_days,
        baseline_solver_inputs=baseline_solver_inputs,
        toxicity_constraints=toxicity_constraints,
        sampled_parameters=parameters,
        sampled_toxicity_coefficients=None,
    )


def _model_parameter_drivers(twin_state) -> list[dict[str, Any]]:
    flat = _flatten_numeric_parameters(dict(twin_state.parameters or {}))
    drivers = []
    for name, value in sorted(flat.items()):
        if name.startswith("assumptions."):
            continue
        drivers.append({
            "name": name,
            "path": name,
            "kind": "model_parameter",
            "baseline_value": float(value),
        })
    return drivers


def _schedule_knob_drivers(intervention_definition: dict[str, Any]) -> list[dict[str, Any]]:
    intervention = dict(intervention_definition.get("intervention") or {})
    drivers = []
    if intervention.get("dose_mg") is not None:
        drivers.append({
            "name": "intervention.dose_mg",
            "path": ["intervention", "dose_mg"],
            "kind": "schedule_knob",
            "baseline_value": float(intervention["dose_mg"]),
        })
    if intervention.get("duration_days") is not None:
        drivers.append({
            "name": "intervention.duration_days",
            "path": ["intervention", "duration_days"],
            "kind": "schedule_knob",
            "baseline_value": float(intervention["duration_days"]),
        })
    for index, phase in enumerate(intervention.get("phases") or []):
        if phase.get("dose_mg") is not None:
            drivers.append({
                "name": f"intervention.phases[{index}].dose_mg",
                "path": ["intervention", "phases", index, "dose_mg"],
                "kind": "schedule_knob",
                "baseline_value": float(phase["dose_mg"]),
            })
    return drivers


def _perturbed_value(baseline_value: float, fraction: float, direction: float, kind: str):
    candidate = float(baseline_value) * (1.0 + (float(fraction) * float(direction)))
    if kind == "schedule_knob":
        return max(int(round(candidate)), 1)
    return max(candidate, 1.0e-9) if baseline_value >= 0 else candidate


def _driver_direction(positive_change) -> str:
    if not positive_change or positive_change.get("delta_utility_v2") is None:
        return "mixed_or_unavailable"
    return "increase_helps" if float(positive_change["delta_utility_v2"]) >= 0.0 else "increase_harms"


def _classify_sensitivity(delta: float | None) -> str:
    if delta is None:
        return "unavailable"
    if float(delta) <= 0.05:
        return "low"
    if float(delta) <= 0.20:
        return "moderate"
    return "high"


def _set_path_value(payload: dict[str, Any], path: list[Any], value) -> None:
    cursor = payload
    for part in path[:-1]:
        if isinstance(part, int):
            cursor = cursor[part]
        else:
            cursor = cursor[part]
    cursor[path[-1]] = value