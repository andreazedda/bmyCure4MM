from __future__ import annotations

import copy
from collections import defaultdict
from types import SimpleNamespace
from typing import Any

from .calibration import (
    CALIBRATION_STATUS_FAILED_NO_IMPROVEMENT,
    CALIBRATION_STATUS_FAILED_UNRELIABLE,
    CALIBRATION_STATUS_USABLE,
    DEFAULT_PARAMETER_BOUNDS,
    _build_parameter_uncertainty_payload,
    _build_residual_diagnostics,
    _get_parameter_value,
    _is_usable_calibration,
    _run_optimizer_sequence,
)
from .observation_model import compute_mae, compute_residuals, compute_rmse, observed_from_assessment, predict_biomarkers
from .simulation_bridge import run_patient_simulation
from .state_model import preview_state_from_assessment
from .therapy_schedule import build_therapy_schedule


COMPLETED_STATUS = "completed"
INSUFFICIENT_HISTORY_STATUS = "insufficient_history"
FAILED_STATUS = "failed"


def run_patient_backtest(patient, twin_state=None, *, minimum_history_points: int = 3) -> dict[str, Any]:
    assessments = list(patient.assessments.order_by("date"))
    if len(assessments) < int(minimum_history_points):
        return {
            "status": INSUFFICIENT_HISTORY_STATUS,
            "message": "At least three historical assessments are required for rolling backtesting.",
            "n_assessments": len(assessments),
            "n_folds": 0,
            "fold_rows": [],
            "overall_rmse": None,
            "overall_mae": None,
            "by_biomarker": {},
        }

    fold_rows = []
    residuals_by_biomarker: dict[str, list[float]] = defaultdict(list)
    all_residuals: list[float] = []
    all_normalized_residuals: list[float] = []
    train_sizes = []

    for fold_index in range(2, len(assessments)):
        training_assessments = assessments[:fold_index]
        held_out_assessment = assessments[fold_index]
        train_sizes.append(len(training_assessments))
        ephemeral_state = _ephemeral_state_from_assessment(training_assessments[0])
        fit_result = _fit_parameters_in_memory(
            patient=patient,
            start_state=ephemeral_state,
            assessments=training_assessments,
            therapies=patient.therapies.all(),
        )
        heldout_payload = _predict_heldout_assessment(
            patient=patient,
            start_state=ephemeral_state,
            fitted_parameters=fit_result["parameters"],
            held_out_assessment=held_out_assessment,
        )
        residual_payload = heldout_payload["residual_payload"]
        if residual_payload["residuals"]:
            for biomarker, value in residual_payload["residuals"].items():
                residuals_by_biomarker[biomarker].append(float(value))
                all_residuals.append(float(value))
            all_normalized_residuals.extend(float(value) for value in residual_payload["normalized_residuals"].values())

        fold_rows.append(
            {
                "fold_index": int(fold_index - 1),
                "status": COMPLETED_STATUS if residual_payload["residuals"] else FAILED_STATUS,
                "train_start_date": training_assessments[0].date.isoformat(),
                "train_end_date": training_assessments[-1].date.isoformat(),
                "test_date": held_out_assessment.date.isoformat(),
                "n_training_assessments": len(training_assessments),
                "fit_status": fit_result["optimizer"]["status"],
                "calibration_rmse": fit_result["diagnostics"].get("rmse"),
                "predicted": residual_payload["predicted"],
                "observed": residual_payload["observed"],
                "residuals": residual_payload["residuals"],
                "normalized_residuals": residual_payload["normalized_residuals"],
                "rmse": residual_payload["rmse"],
                "mae": residual_payload["mae"],
                "calibration_fit_window": {
                    "start_date": training_assessments[0].date.isoformat(),
                    "end_date": training_assessments[-1].date.isoformat(),
                    "n_assessments": len(training_assessments),
                },
            }
        )

    by_biomarker = {}
    for biomarker, values in sorted(residuals_by_biomarker.items()):
        by_biomarker[biomarker] = {
            "n_points": len(values),
            "rmse": compute_rmse({"residuals": {str(index): value for index, value in enumerate(values)}}),
            "mae": compute_mae({str(index): value for index, value in enumerate(values)}),
        }

    overall_rmse = None
    overall_mae = None
    if all_residuals:
        indexed = {str(index): value for index, value in enumerate(all_residuals)}
        overall_rmse = compute_rmse({"residuals": indexed})
        overall_mae = compute_mae(indexed)

    return {
        "status": COMPLETED_STATUS,
        "message": "Rolling-origin backtesting compares one-step-ahead held-out biomarker predictions against observed assessment points.",
        "n_assessments": len(assessments),
        "n_folds": len(fold_rows),
        "training_assessment_range": {
            "min": min(train_sizes) if train_sizes else 0,
            "max": max(train_sizes) if train_sizes else 0,
        },
        "overall_rmse": overall_rmse,
        "overall_mae": overall_mae,
        "normalized_residual_mae": sum(abs(value) for value in all_normalized_residuals) / len(all_normalized_residuals) if all_normalized_residuals else None,
        "by_biomarker": by_biomarker,
        "fold_rows": fold_rows,
        "limitations": [
            "Rolling backtesting measures held-out agreement at observed biomarker points only.",
            "Backtest diagnostics remain research validation signals and do not establish clinical validity or causal effect estimation.",
        ],
    }


def _fit_parameters_in_memory(*, patient, start_state, assessments, therapies, bounds: dict[str, tuple[float, float]] | None = None) -> dict[str, Any]:
    assessments = sorted(list(assessments), key=lambda item: item.date)
    parameter_bounds = dict(DEFAULT_PARAMETER_BOUNDS)
    if bounds:
        parameter_bounds.update(bounds)
    parameter_names = list(parameter_bounds)
    baseline_parameters = copy.deepcopy(dict(start_state.parameters or {}))
    initial_theta = [_get_parameter_value(start_state, name) for name in parameter_names]
    baseline_diagnostics = _build_residual_diagnostics(
        patient=patient,
        start_state=start_state,
        assessments=assessments,
        therapies=therapies,
        overrides=baseline_parameters,
    )

    best_attempt = None
    try:
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
    except Exception:
        best_attempt = None

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

    return {
        "parameters": chosen_parameters,
        "baseline_diagnostics": baseline_diagnostics,
        "diagnostics": chosen_diagnostics,
        "optimizer": {
            "optimizer": optimizer_name,
            "success": optimizer_success,
            "message": optimizer_message,
            "status": optimizer_status,
        },
        "parameter_uncertainty": _build_parameter_uncertainty_payload(
            optimizer=optimizer_name,
            optimizer_success=optimizer_success,
            optimizer_message=optimizer_message,
            baseline_diagnostics=baseline_diagnostics,
            calibrated_diagnostics=chosen_diagnostics,
            n_parameters=len(parameter_names),
            calibration_status=optimizer_status,
        ),
    }


def _predict_heldout_assessment(*, patient, start_state, fitted_parameters: dict[str, Any], held_out_assessment) -> dict[str, Any]:
    start_date = start_state.state_date
    horizon_days = max((held_out_assessment.date - start_date).days, 1)
    therapy_schedule = build_therapy_schedule(patient, start_date, held_out_assessment.date)
    simulation_result = run_patient_simulation(
        start_state,
        therapy_schedule=therapy_schedule,
        horizon_days=horizon_days,
        overrides=fitted_parameters,
    )
    predicted = predict_biomarkers(
        simulation_result["trajectory_frame"],
        start_state.assessment or held_out_assessment,
        max((held_out_assessment.date - start_date).days, 0),
        calibration_parameters=fitted_parameters,
    )
    observed = observed_from_assessment(held_out_assessment)
    residual_payload = compute_residuals(predicted, observed)
    return {
        "simulation_result": simulation_result,
        "predicted": predicted,
        "observed": observed,
        "residual_payload": residual_payload,
    }


def _ephemeral_state_from_assessment(assessment):
    preview = preview_state_from_assessment(assessment)
    return SimpleNamespace(
        id=None,
        patient=assessment.patient,
        patient_id=assessment.patient_id,
        assessment=assessment,
        assessment_id=assessment.id,
        state_date=assessment.date,
        state_vector=preview["state_vector"],
        parameters=preview["parameters"],
        parameter_uncertainty={},
        risk_score=preview["twin_payload"].get("risk_score"),
        created_by=None,
    )