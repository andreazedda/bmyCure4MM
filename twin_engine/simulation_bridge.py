from __future__ import annotations

from typing import Any

import numpy as np

from simulator.models import DEFAULT_PD_PARAMS, DEFAULT_PK_PARAMS, SimulationAttempt
from simulator.models_simulation import MathematicalModel
from simulator.pharmaco import registry as pharmaco_registry

from .observation_model import build_predicted_biomarker_projection
from .therapy_schedule import convert_patient_therapies_to_drug_doses


def build_solver_inputs_from_twin_state(
    twin_state,
    therapy_schedule,
    horizon_days: int,
    overrides: dict[str, Any] | None = None,
) -> dict[str, Any]:
    overrides = dict(overrides or {})
    parameters = dict(twin_state.parameters or {})
    state_vector = dict(twin_state.state_vector or {})

    raw_parameters = {
        "baseline_tumor_cells": state_vector.get("tumor_cells", 1.0e9),
        "baseline_healthy_cells": state_vector.get("healthy_cells", 5.0e11),
        "time_horizon": horizon_days,
        "tumor_growth_rate": parameters.get("tumor_growth_rate", 0.023),
        "healthy_growth_rate": parameters.get("healthy_growth_rate", 0.015),
        "carrying_capacity_tumor": parameters.get("carrying_capacity_tumor"),
        "carrying_capacity_healthy": parameters.get("carrying_capacity_healthy"),
        "immune_compromise_index": parameters.get("immune_compromise_index", 1.0),
        "interaction_strength": parameters.get("interaction_strength", 0.05),
        "lenalidomide_dose": 0.0,
        "bortezomib_dose": 0.0,
        "daratumumab_dose": 0.0,
        "carfilzomib_dose": 0.0,
    }

    missing_doses: list[dict[str, Any]] = []
    if therapy_schedule:
        schedule_payload = convert_patient_therapies_to_drug_doses(therapy_schedule, strict=False)
        missing_doses = list(schedule_payload.get("missing_doses", []))
        for drug, dose in schedule_payload.get("drug_doses", {}).items():
            raw_parameters[f"{drug}_dose"] = float(dose)

    for drug, dose in (overrides.pop("drug_doses", {}) or {}).items():
        raw_parameters[f"{drug}_dose"] = float(dose)

    raw_parameters.update(overrides)

    solver_inputs = SimulationAttempt._resolve_solver_inputs(raw_parameters)
    solver_inputs["raw_parameters"] = raw_parameters
    solver_inputs["missing_doses"] = missing_doses
    solver_inputs["therapy_schedule"] = therapy_schedule
    return solver_inputs


def run_patient_simulation(
    twin_state,
    therapy_schedule=None,
    horizon_days: int = 90,
    overrides: dict[str, Any] | None = None,
) -> dict[str, Any]:
    solver_inputs = build_solver_inputs_from_twin_state(
        twin_state,
        therapy_schedule,
        horizon_days,
        overrides=overrides,
    )

    pk_params, pd_params, dose_functions = pharmaco_registry.resolve(
        solver_inputs["drug_doses"],
        solver_inputs["time_horizon"],
        DEFAULT_PK_PARAMS,
        DEFAULT_PD_PARAMS,
    )
    interaction_matrix = np.eye(len(solver_inputs["drug_doses"]), dtype=float) * float(
        solver_inputs["interaction_strength"]
    )

    model = MathematicalModel(
        baseline_tumor_cells=float(solver_inputs["baseline_tumor_cells"]),
        baseline_healthy_cells=float(solver_inputs["baseline_healthy_cells"]),
        drug_doses=solver_inputs["drug_doses"],
        pk_params=pk_params,
        pd_params=pd_params,
        growth_rates=solver_inputs["growth_rates"],
        interaction_matrix=interaction_matrix,
        time_span=(0.0, float(solver_inputs["time_horizon"])),
        carrying_capacity_tumor=solver_inputs["carrying_capacity_tumor"],
        carrying_capacity_healthy=solver_inputs["carrying_capacity_healthy"],
        immune_compromise_index=float(solver_inputs["immune_compromise_index"]),
        dose_functions=dose_functions,
    )
    trajectory = model.simulate()
    summary = summarize_trajectory(
        trajectory,
        time_horizon_days=float(solver_inputs["time_horizon"]),
        tumor_baseline=float(solver_inputs["baseline_tumor_cells"]),
        healthy_baseline=float(solver_inputs["baseline_healthy_cells"]),
    )
    baseline_assessment = getattr(twin_state, "assessment", None)
    if baseline_assessment is None and hasattr(twin_state, "source_assessments"):
        baseline_assessment = twin_state.source_assessments.order_by("date").first()
    predicted_projection = build_predicted_biomarker_projection(
        trajectory,
        baseline_assessment,
        prediction_day=float(solver_inputs["time_horizon"]),
        calibration_parameters=dict(twin_state.parameters or {}) | {
            key: value
            for key, value in (overrides or {}).items()
            if key not in {"drug_doses", "seed"}
        },
        model_version=getattr(twin_state, "model_version", None),
        config_hash=getattr(twin_state, "config_hash", None),
        milestone_days={
            label: float(item["day"])
            for label, item in summary.get("milestones", {}).items()
        },
    )
    summary["predicted_biomarkers"] = predicted_projection["end"]
    summary["predicted_biomarkers_milestones"] = predicted_projection["milestones"]
    return {
        "trajectory_frame": trajectory,
        "trajectory": trajectory.to_dict(orient="list"),
        "summary": summary,
        "solver_inputs": solver_inputs,
    }


def summarize_trajectory(trajectory, *, time_horizon_days: float, tumor_baseline: float, healthy_baseline: float) -> dict[str, Any]:
    tumor_series = trajectory["tumor_cells"]
    healthy_series = trajectory["healthy_cells"]
    time_series = trajectory["time"]

    tumor_start = float(tumor_series.iloc[0])
    tumor_end = float(tumor_series.iloc[-1])
    healthy_start = float(healthy_series.iloc[0])
    healthy_end = float(healthy_series.iloc[-1])

    tumor_reduction = float(1.0 - tumor_end / max(tumor_start, 1.0e-9))
    healthy_loss = float(1.0 - healthy_end / max(healthy_start, 1.0e-9))

    time_values = np.asarray(time_series, dtype=float)
    tumor_values = np.asarray(tumor_series, dtype=float)
    healthy_values = np.asarray(healthy_series, dtype=float)

    nadir_index = int(np.argmin(tumor_values))
    nadir_time = float(time_values[nadir_index])
    nadir_tumor = float(tumor_values[nadir_index])

    milestones: dict[str, dict[str, float]] = {}
    for day in (30.0, 60.0, 90.0, float(time_horizon_days)):
        if day > float(time_horizon_days):
            continue
        position = int(np.argmin(np.abs(time_values - float(day))))
        tumor_at_day = float(tumor_values[position])
        healthy_at_day = float(healthy_values[position])
        label = "end" if day == float(time_horizon_days) else f"day_{int(day)}"
        milestones[label] = {
            "day": float(time_values[position]),
            "tumor_cells": tumor_at_day,
            "healthy_cells": healthy_at_day,
            "tumor_reduction": float(1.0 - tumor_at_day / max(tumor_start, 1.0e-9)),
            "healthy_loss": float(1.0 - healthy_at_day / max(healthy_start, 1.0e-9)),
        }

    recurrence_threshold = 0.5 * float(tumor_start)
    time_to_recurrence = None
    for index in range(nadir_index + 1, len(time_values)):
        if float(tumor_values[index]) > recurrence_threshold:
            time_to_recurrence = float(time_values[index])
            break

    auc = {}
    for column in trajectory.columns:
        if column.endswith("_concentration"):
            auc[column.replace("_concentration", "")] = _integrate_trapezoid(
                np.asarray(trajectory[column], dtype=float),
                time_values,
            )

    return {
        "label": "research simulation",
        "counterfactual_type": "mechanistic model counterfactual",
        "causal_interpretation": "unvalidated exploratory hypothesis",
        "tumor_reduction": tumor_reduction,
        "healthy_loss": healthy_loss,
        "time_to_recurrence": time_to_recurrence,
        "durability_index": _integrate_trapezoid(
            (tumor_values < float(tumor_baseline)).astype(float),
            time_values,
        ) / max(float(time_horizon_days), 1.0e-9),
        "nadir": {
            "day": nadir_time,
            "tumor_cells": nadir_tumor,
            "tumor_reduction": float(1.0 - nadir_tumor / max(tumor_start, 1.0e-9)),
        },
        "milestones": milestones,
        "auc": auc,
        "baseline_tumor_cells": float(tumor_baseline),
        "baseline_healthy_cells": float(healthy_baseline),
    }


def _integrate_trapezoid(y_values, x_values) -> float:
    trapezoid_fn = getattr(np, "trapezoid", None)
    if callable(trapezoid_fn):
        return float(trapezoid_fn(y_values, x_values))
    return float(np.trapz(y_values, x_values))
