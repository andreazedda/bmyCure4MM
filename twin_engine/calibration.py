from __future__ import annotations

import copy
import math
from typing import Any

from django.core.exceptions import ValidationError
from django.db import transaction

from .models import ObservationResidual, PatientTwinState
from .observation_model import compute_residuals, observed_from_assessment, predict_biomarkers
from .provenance import CURRENT_MODEL_VERSION, collect_twin_config_hash
from .simulation_bridge import run_patient_simulation
from .therapy_schedule import build_therapy_schedule


EPSILON = 1.0e-9
CALIBRATION_STATUS_USABLE = "usable"
CALIBRATION_STATUS_FAILED_UNRELIABLE = "failed_or_unreliable"
CALIBRATION_STATUS_FAILED_NO_IMPROVEMENT = "failed_no_improvement"

DEFAULT_PARAMETER_BOUNDS = {
    "tumor_growth_rate": (0.001, 0.10),
    "immune_compromise_index": (0.1, 5.0),
    "carrying_capacity_tumor": (1.0e7, 1.0e11),
    "observation__alpha_M": (1.0e-12, 1.0e-8),
    "observation__beta_M": (0.0, 2.0),
}


def calibrate_patient_parameters(patient, start_state, assessments, therapies, bounds: dict[str, tuple[float, float]] | None = None):
    assessments = sorted(list(assessments), key=lambda item: item.date)
    if len(assessments) < 2:
        raise ValidationError("Calibration requires at least two assessments.")

    parameter_bounds = dict(DEFAULT_PARAMETER_BOUNDS)
    if bounds:
        parameter_bounds.update(bounds)
    parameter_names = list(parameter_bounds)
    baseline_parameters = copy.deepcopy(dict(start_state.parameters or {}))
    initial_theta = [_get_parameter_value_from_parameters(baseline_parameters, name) for name in parameter_names]

    baseline_diagnostics = _build_residual_diagnostics(
        patient=patient,
        start_state=start_state,
        assessments=assessments,
        therapies=therapies,
        overrides=baseline_parameters,
    )

    with transaction.atomic():
        _upsert_residual_records(
            patient=patient,
            twin_state=start_state,
            residual_entries=baseline_diagnostics["residuals"],
            stage=ObservationResidual.STAGE_PRE_CALIBRATION,
        )

    best_attempt = _run_optimizer_sequence(
        parameter_names=parameter_names,
        parameter_bounds=parameter_bounds,
        initial_theta=initial_theta,
        patient=patient,
        start_state=start_state,
        assessments=assessments,
        therapies=therapies,
        baseline_diagnostics=baseline_diagnostics,
    )

    if best_attempt is None or float(best_attempt["diagnostics"]["objective"]) >= float(baseline_diagnostics["objective"]):
        chosen_parameters = baseline_parameters
        chosen_diagnostics = baseline_diagnostics
        optimizer_name = best_attempt["optimizer"] if best_attempt else "none"
        optimizer_success = False
        optimizer_message = best_attempt["message"] if best_attempt else "No optimizer produced an improved solution."
        optimizer_status = CALIBRATION_STATUS_FAILED_NO_IMPROVEMENT
    else:
        chosen_parameters = best_attempt["parameters"]
        chosen_diagnostics = best_attempt["diagnostics"]
        optimizer_name = best_attempt["optimizer"]
        optimizer_success = bool(best_attempt["success"])
        optimizer_message = str(best_attempt["message"])
        optimizer_status = (
            CALIBRATION_STATUS_USABLE
            if _is_usable_calibration(
                optimizer_success=optimizer_success,
                baseline_diagnostics=baseline_diagnostics,
                calibrated_diagnostics=chosen_diagnostics,
                n_parameters=len(parameter_names),
            )
            else CALIBRATION_STATUS_FAILED_UNRELIABLE
        )

    diagnostics_payload = _build_parameter_uncertainty_payload(
        optimizer=optimizer_name,
        optimizer_success=optimizer_success,
        optimizer_message=optimizer_message,
        baseline_diagnostics=baseline_diagnostics,
        calibrated_diagnostics=chosen_diagnostics,
        n_parameters=len(parameter_names),
        calibration_status=optimizer_status,
    )

    with transaction.atomic():
        calibrated_state = PatientTwinState.objects.create(
            patient=patient,
            assessment=assessments[-1],
            state_date=assessments[-1].date,
            is_current=False,
            state_vector=start_state.state_vector,
            parameters=chosen_parameters,
            parameter_uncertainty=diagnostics_payload,
            risk_score=start_state.risk_score,
            method=PatientTwinState.METHOD_RESIDUAL_MINIMIZATION,
            model_version=CURRENT_MODEL_VERSION,
            config_hash=collect_twin_config_hash(),
            created_by=start_state.created_by,
        )
        calibrated_state.source_assessments.add(*assessments)
        from .state_model import set_current_state

        set_current_state(calibrated_state)

        residual_records = _upsert_residual_records(
            patient=patient,
            twin_state=calibrated_state,
            residual_entries=chosen_diagnostics["residuals"],
            stage=ObservationResidual.STAGE_POST_CALIBRATION,
        )

    return {
        "state": calibrated_state,
        "residuals": residual_records,
        "baseline_diagnostics": baseline_diagnostics,
        "diagnostics": chosen_diagnostics,
        "optimizer": {
            "optimizer": optimizer_name,
            "success": optimizer_success,
            "message": optimizer_message,
            "status": optimizer_status,
            "objective_before": baseline_diagnostics["objective"],
            "objective_after": chosen_diagnostics["objective"],
            "rmse_before": baseline_diagnostics["rmse"],
            "rmse_after": chosen_diagnostics["rmse"],
            "mae_before": baseline_diagnostics["mae"],
            "mae_after": chosen_diagnostics["mae"],
        },
    }


def objective_function(theta, parameter_names, patient, start_state, assessments, therapies, parameter_bounds=None) -> float:
    parameter_overrides = _parameters_from_theta(
        start_state.parameters or {},
        parameter_names,
        theta,
        parameter_bounds=parameter_bounds,
    )

    diagnostics = _build_residual_diagnostics(
        patient=patient,
        start_state=start_state,
        assessments=assessments,
        therapies=therapies,
        overrides=parameter_overrides,
    )
    penalty = _bounds_penalty(theta, parameter_names, parameter_bounds or {})
    return float(diagnostics["objective"]) + penalty


def _build_residual_diagnostics(patient, start_state, assessments, therapies, overrides: dict[str, Any]) -> dict[str, Any]:
    start_date = start_state.state_date
    end_date = assessments[-1].date
    horizon_days = max((end_date - start_date).days, 1)

    therapy_schedule = build_therapy_schedule(patient, start_date, end_date)
    simulation_result = run_patient_simulation(
        start_state,
        therapy_schedule=therapy_schedule,
        horizon_days=horizon_days,
        overrides=overrides,
    )

    objective = 0.0
    weighted_square_sum = 0.0
    weighted_abs_sum = 0.0
    total_weight = 0.0
    residual_vector: list[float] = []
    weighted_residual_vector: list[float] = []
    data_points_by_analyte: dict[str, int] = {}
    residual_entries: list[dict[str, Any]] = []
    for assessment in assessments:
        offset_days = max((assessment.date - start_date).days, 0)
        predicted = predict_biomarkers(
            simulation_result["trajectory_frame"],
            start_state.assessment or assessments[0],
            offset_days,
            calibration_parameters=overrides,
        )
        observed = observed_from_assessment(assessment)
        residual_payload = compute_residuals(predicted, observed)
        for name, residual in residual_payload["residuals"].items():
            weight = float(residual_payload["weights"].get(name, 1.0))
            objective += weight * (float(residual) ** 2)
            weighted_square_sum += weight * (float(residual) ** 2)
            weighted_abs_sum += weight * abs(float(residual))
            total_weight += weight
            residual_vector.append(float(residual))
            weighted_residual_vector.append((weight ** 0.5) * float(residual))
            data_points_by_analyte[name] = data_points_by_analyte.get(name, 0) + 1
        residual_entries.append(
            {
                "assessment": assessment,
                "residual_payload": residual_payload,
            }
        )

    return {
        "objective": objective,
        "simulation_result": simulation_result,
        "residuals": residual_entries,
        "rmse": math.sqrt(weighted_square_sum / total_weight) if total_weight else None,
        "mae": (weighted_abs_sum / total_weight) if total_weight else None,
        "n_observations": sum(data_points_by_analyte.values()),
        "data_points_by_analyte": data_points_by_analyte,
        "residual_vector": residual_vector,
        "weighted_residual_vector": weighted_residual_vector,
    }


def _get_parameter_value(state, name: str) -> float:
    parameters = dict(state.parameters or {})
    if name in {"observation__alpha_M", "observation__beta_M"}:
        observation = parameters.get("observation", {})
        observation_key = name.split("__", 1)[1]
        default_value = 1.0e-9 if observation_key == "alpha_M" else 0.0
        return float(observation.get(observation_key, default_value))
    value = parameters.get(name)
    if value is None:
        value = DEFAULT_PARAMETER_BOUNDS[name][0]
    return float(value)


def _set_parameter_value(parameters: dict[str, Any], name: str, value: float) -> None:
    if name in {"observation__alpha_M", "observation__beta_M"}:
        observation = dict(parameters.get("observation", {}))
        observation[name.split("__", 1)[1]] = float(value)
        parameters["observation"] = observation
        return
    parameters[name] = float(value)


def _run_optimizer_sequence(
    *,
    parameter_names,
    parameter_bounds,
    initial_theta,
    patient,
    start_state,
    assessments,
    therapies,
    baseline_diagnostics,
):
    try:
        from scipy.optimize import differential_evolution, least_squares, minimize
    except Exception as exc:  # pragma: no cover
        raise NotImplementedError("scipy.optimize is required for calibration.") from exc

    lower_bounds = [parameter_bounds[name][0] for name in parameter_names]
    upper_bounds = [parameter_bounds[name][1] for name in parameter_names]

    attempts: list[dict[str, Any]] = []

    try:
        result = least_squares(
            _residual_vector_objective,
            x0=initial_theta,
            bounds=(lower_bounds, upper_bounds),
            args=(parameter_names, patient, start_state, assessments, therapies, parameter_bounds),
        )
        attempt = _build_attempt_payload(
            optimizer="least_squares",
            theta=result.x,
            success=bool(result.success),
            message=str(result.message),
            parameter_names=parameter_names,
            parameter_bounds=parameter_bounds,
            patient=patient,
            start_state=start_state,
            assessments=assessments,
            therapies=therapies,
        )
        attempts.append(attempt)
        if _attempt_improves_over_baseline(attempt, baseline_diagnostics):
            return attempt
    except Exception as exc:
        attempts.append({"optimizer": "least_squares", "success": False, "message": str(exc), "diagnostics": baseline_diagnostics})

    try:
        result = minimize(
            objective_function,
            x0=initial_theta,
            method="Nelder-Mead",
            args=(parameter_names, patient, start_state, assessments, therapies, parameter_bounds),
        )
        attempt = _build_attempt_payload(
            optimizer="Nelder-Mead",
            theta=result.x,
            success=bool(result.success),
            message=str(result.message),
            parameter_names=parameter_names,
            parameter_bounds=parameter_bounds,
            patient=patient,
            start_state=start_state,
            assessments=assessments,
            therapies=therapies,
        )
        attempts.append(attempt)
        if _attempt_improves_over_baseline(attempt, baseline_diagnostics):
            return attempt
    except Exception as exc:
        attempts.append({"optimizer": "Nelder-Mead", "success": False, "message": str(exc), "diagnostics": baseline_diagnostics})

    if baseline_diagnostics["n_observations"] >= (len(parameter_names) * 2):
        try:
            result = differential_evolution(
                objective_function,
                bounds=[parameter_bounds[name] for name in parameter_names],
                args=(parameter_names, patient, start_state, assessments, therapies, parameter_bounds),
                maxiter=20,
                polish=False,
                seed=17,
            )
            attempt = _build_attempt_payload(
                optimizer="differential_evolution",
                theta=result.x,
                success=bool(result.success),
                message=str(result.message),
                parameter_names=parameter_names,
                parameter_bounds=parameter_bounds,
                patient=patient,
                start_state=start_state,
                assessments=assessments,
                therapies=therapies,
            )
            attempts.append(attempt)
            if _attempt_improves_over_baseline(attempt, baseline_diagnostics):
                return attempt
        except Exception as exc:
            attempts.append({"optimizer": "differential_evolution", "success": False, "message": str(exc), "diagnostics": baseline_diagnostics})

    viable_attempts = [item for item in attempts if "parameters" in item]
    if not viable_attempts:
        return None
    return min(viable_attempts, key=lambda item: float(item["diagnostics"]["objective"]))


def _build_attempt_payload(
    *,
    optimizer,
    theta,
    success,
    message,
    parameter_names,
    parameter_bounds,
    patient,
    start_state,
    assessments,
    therapies,
):
    parameters = _parameters_from_theta(start_state.parameters or {}, parameter_names, theta, parameter_bounds=parameter_bounds)
    diagnostics = _build_residual_diagnostics(
        patient=patient,
        start_state=start_state,
        assessments=assessments,
        therapies=therapies,
        overrides=parameters,
    )
    return {
        "optimizer": optimizer,
        "success": success,
        "message": message,
        "parameters": parameters,
        "diagnostics": diagnostics,
    }


def _residual_vector_objective(theta, parameter_names, patient, start_state, assessments, therapies, parameter_bounds):
    parameters = _parameters_from_theta(start_state.parameters or {}, parameter_names, theta, parameter_bounds=parameter_bounds)
    diagnostics = _build_residual_diagnostics(
        patient=patient,
        start_state=start_state,
        assessments=assessments,
        therapies=therapies,
        overrides=parameters,
    )
    return diagnostics["weighted_residual_vector"] or [0.0]


def _parameters_from_theta(base_parameters, parameter_names, theta, parameter_bounds=None) -> dict[str, Any]:
    parameters = copy.deepcopy(dict(base_parameters or {}))
    for name, value in zip(parameter_names, theta):
        bounded_value = _clip_to_bounds(name, float(value), parameter_bounds or {})
        _set_parameter_value(parameters, name, bounded_value)
    return parameters


def _clip_to_bounds(name: str, value: float, parameter_bounds: dict[str, tuple[float, float]]) -> float:
    lower, upper = parameter_bounds.get(name, (value, value))
    return min(max(float(value), float(lower)), float(upper))


def _bounds_penalty(theta, parameter_names, parameter_bounds: dict[str, tuple[float, float]]) -> float:
    penalty = 0.0
    for name, value in zip(parameter_names, theta):
        lower, upper = parameter_bounds.get(name, (value, value))
        if float(value) < float(lower):
            penalty += (float(lower) - float(value)) ** 2 * 1.0e6
        if float(value) > float(upper):
            penalty += (float(value) - float(upper)) ** 2 * 1.0e6
    return penalty


def _is_usable_calibration(*, optimizer_success: bool, baseline_diagnostics, calibrated_diagnostics, n_parameters: int) -> bool:
    objective_before = float(baseline_diagnostics["objective"])
    objective_after = float(calibrated_diagnostics["objective"])
    rmse_before = baseline_diagnostics.get("rmse")
    rmse_after = calibrated_diagnostics.get("rmse")
    enough_observations = int(calibrated_diagnostics.get("n_observations") or 0) >= int(n_parameters) + 1
    objective_improved = objective_after < objective_before
    rmse_improved = (
        rmse_before is not None
        and rmse_after is not None
        and float(rmse_after) < float(rmse_before)
    )
    return bool((optimizer_success or objective_improved) and rmse_improved and enough_observations)


def _attempt_improves_over_baseline(attempt, baseline_diagnostics) -> bool:
    diagnostics = attempt.get("diagnostics") or {}
    return float(diagnostics.get("objective") or 0.0) < float(baseline_diagnostics.get("objective") or 0.0)


def _build_parameter_uncertainty_payload(
    *,
    optimizer,
    optimizer_success,
    optimizer_message,
    baseline_diagnostics,
    calibrated_diagnostics,
    n_parameters: int,
    calibration_status: str,
) -> dict[str, Any]:
    rmse_before = baseline_diagnostics.get("rmse")
    rmse_after = calibrated_diagnostics.get("rmse")
    mae_before = baseline_diagnostics.get("mae")
    mae_after = calibrated_diagnostics.get("mae")
    improvement_ratio = None
    if rmse_before is not None and rmse_after is not None:
        improvement_ratio = (float(rmse_before) - float(rmse_after)) / max(float(rmse_before), EPSILON)
    return {
        "optimizer": optimizer,
        "optimizer_success": bool(optimizer_success),
        "success": calibration_status == CALIBRATION_STATUS_USABLE,
        "status": calibration_status,
        "calibration_status": calibration_status,
        "message": str(optimizer_message),
        "objective_before": float(baseline_diagnostics["objective"]),
        "objective_after": float(calibrated_diagnostics["objective"]),
        "rmse_before": rmse_before,
        "rmse_after": rmse_after,
        "mae_before": mae_before,
        "mae_after": mae_after,
        "n_observations": int(calibrated_diagnostics.get("n_observations") or 0),
        "n_parameters": int(n_parameters),
        "data_points_by_analyte": dict(calibrated_diagnostics.get("data_points_by_analyte") or {}),
        "improvement_ratio": improvement_ratio,
    }


def _upsert_residual_records(*, patient, twin_state, residual_entries, stage: str):
    records = []
    for item in residual_entries:
        payload = item["residual_payload"]
        assessment = item["assessment"]
        record = ObservationResidual.objects.filter(
            patient=patient,
            twin_state=twin_state,
            assessment=assessment,
            stage=stage,
        ).first()
        if record is None:
            record = ObservationResidual.objects.create(
                patient=patient,
                twin_state=twin_state,
                assessment=assessment,
                predicted_values=payload["predicted"],
                observed_values=payload["observed"],
                residuals=payload["residuals"],
                normalized_residuals=payload["normalized_residuals"],
                rmse=payload["rmse"],
                mae=payload["mae"],
                biomarker_weights=payload["weights"],
                stage=stage,
            )
        else:
            record.predicted_values = payload["predicted"]
            record.observed_values = payload["observed"]
            record.residuals = payload["residuals"]
            record.normalized_residuals = payload["normalized_residuals"]
            record.rmse = payload["rmse"]
            record.mae = payload["mae"]
            record.biomarker_weights = payload["weights"]
            record.stage = stage
            record.save(
                update_fields=[
                    "predicted_values",
                    "observed_values",
                    "residuals",
                    "normalized_residuals",
                    "rmse",
                    "mae",
                    "biomarker_weights",
                    "stage",
                ]
            )
        records.append(record)
    return records


def _get_parameter_value_from_parameters(parameters: dict[str, Any], name: str) -> float:
    if name in {"observation__alpha_M", "observation__beta_M"}:
        observation = dict(parameters.get("observation", {}))
        observation_key = name.split("__", 1)[1]
        default_value = 1.0e-9 if observation_key == "alpha_M" else 0.0
        return float(observation.get(observation_key, default_value))
    value = parameters.get(name)
    if value is None:
        value = DEFAULT_PARAMETER_BOUNDS[name][0]
    return float(value)
