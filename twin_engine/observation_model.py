from __future__ import annotations

from typing import Any


DEFAULT_OBSERVATION_PARAMETERS = {
    "alpha_M": 1.0e-9,
    "beta_M": 0.0,
    "alpha_F": 1.0,
    "gamma_F": 1.0,
    "T_ref": 1.0e9,
    "Hb_baseline": 12.0,
    "H_ref": 5.0e11,
    "eta_H": 0.35,
    "eps": 1.0e-9,
}


def default_observation_parameters(assessment=None, *, tumor_cells: float | None = None, healthy_cells: float | None = None) -> dict[str, float]:
    payload = dict(DEFAULT_OBSERVATION_PARAMETERS)

    if assessment is not None:
        hemoglobin = _maybe_float(getattr(assessment, "hemoglobin_g_dl", None))
        m_protein = _maybe_float(getattr(assessment, "m_protein_g_dl", None))
        flc_ratio = _maybe_float(getattr(assessment, "flc_ratio", None))
        if hemoglobin is not None:
            payload["Hb_baseline"] = hemoglobin
        if tumor_cells is not None and tumor_cells > 0 and m_protein is not None:
            payload["alpha_M"] = max(m_protein / tumor_cells, payload["alpha_M"])
        if flc_ratio is not None and flc_ratio > 0:
            payload["alpha_F"] = flc_ratio

    if tumor_cells is not None and tumor_cells > 0:
        payload["T_ref"] = float(tumor_cells)
    if healthy_cells is not None and healthy_cells > 0:
        payload["H_ref"] = float(healthy_cells)

    return payload


def predict_biomarkers(
    trajectory,
    baseline_assessment,
    target_assessment_date_time_offset,
    calibration_parameters: dict[str, Any] | None = None,
) -> dict[str, float]:
    state = _trajectory_state_at(trajectory, float(target_assessment_date_time_offset))
    tumor_cells = float(state["tumor_cells"])
    healthy_cells = float(state["healthy_cells"])

    parameters = default_observation_parameters(
        baseline_assessment,
        tumor_cells=tumor_cells,
        healthy_cells=healthy_cells,
    )
    if calibration_parameters:
        observation_params = calibration_parameters.get("observation", calibration_parameters)
        parameters.update({
            key: float(value)
            for key, value in observation_params.items()
            if key in parameters and value is not None
        })

    eps = float(parameters.get("eps", 1.0e-9))
    m_hat = float(parameters["alpha_M"]) * tumor_cells + float(parameters["beta_M"])
    flc_hat = float(parameters["alpha_F"]) * (
        tumor_cells / max(float(parameters["T_ref"]), eps)
    ) ** float(parameters["gamma_F"])
    hb_hat = float(parameters["Hb_baseline"]) * (
        healthy_cells / max(float(parameters["H_ref"]), eps)
    ) ** float(parameters["eta_H"])

    return {
        "m_protein_g_dl": m_hat,
        "flc_ratio": flc_hat,
        "hemoglobin_g_dl": hb_hat,
    }


def build_predicted_biomarker_projection(
    trajectory,
    baseline_assessment,
    *,
    prediction_day: float,
    calibration_parameters: dict[str, Any] | None = None,
    model_version: str | None = None,
    config_hash: str | None = None,
    milestone_days: dict[str, float] | None = None,
) -> dict[str, Any]:
    if baseline_assessment is None:
        empty = _decorate_prediction_payload({}, prediction_day, model_version=model_version, config_hash=config_hash)
        return {"end": empty, "milestones": {}}

    milestone_days = dict(milestone_days or {})
    end_prediction = predict_biomarkers(
        trajectory,
        baseline_assessment,
        prediction_day,
        calibration_parameters,
    )
    milestone_predictions = {
        label: _decorate_prediction_payload(
            predict_biomarkers(
                trajectory,
                baseline_assessment,
                day,
                calibration_parameters,
            ),
            day,
            model_version=model_version,
            config_hash=config_hash,
        )
        for label, day in milestone_days.items()
    }
    return {
        "end": _decorate_prediction_payload(
            end_prediction,
            prediction_day,
            model_version=model_version,
            config_hash=config_hash,
        ),
        "milestones": milestone_predictions,
    }


def observed_from_assessment(assessment) -> dict[str, float | None]:
    return {
        "m_protein_g_dl": _maybe_float(getattr(assessment, "m_protein_g_dl", None)),
        "flc_ratio": _maybe_float(getattr(assessment, "flc_ratio", None)),
        "hemoglobin_g_dl": _maybe_float(getattr(assessment, "hemoglobin_g_dl", None)),
    }


def compute_residuals(
    predicted: dict[str, float | None],
    observed: dict[str, float | None],
    weights: dict[str, float] | None = None,
) -> dict[str, Any]:
    weights = dict(weights or {})
    residuals: dict[str, float] = {}
    normalized_residuals: dict[str, float] = {}

    for name, predicted_value in predicted.items():
        observed_value = observed.get(name)
        if predicted_value is None or observed_value is None:
            continue
        residual = float(observed_value) - float(predicted_value)
        residuals[name] = residual
        normalized_residuals[name] = residual / max(abs(float(observed_value)), 1.0e-6)
        weights.setdefault(name, 1.0)

    return {
        "predicted": predicted,
        "observed": observed,
        "residuals": residuals,
        "normalized_residuals": normalized_residuals,
        "weights": weights,
        "rmse": compute_rmse(residuals, weights=weights),
        "mae": compute_mae(residuals, weights=weights),
    }


def compute_rmse(residuals, weights: dict[str, float] | None = None) -> float | None:
    residual_payload = residuals.get("residuals") if isinstance(residuals, dict) and "residuals" in residuals else residuals
    if not residual_payload:
        return None
    weights = dict(weights or {})
    numerator = 0.0
    denominator = 0.0
    for name, value in residual_payload.items():
        weight = float(weights.get(name, 1.0))
        numerator += weight * (float(value) ** 2)
        denominator += weight
    if denominator <= 0:
        return None
    return (numerator / denominator) ** 0.5


def compute_mae(residuals, weights: dict[str, float] | None = None) -> float | None:
    if not residuals:
        return None
    weights = dict(weights or {})
    numerator = 0.0
    denominator = 0.0
    for name, value in residuals.items():
        weight = float(weights.get(name, 1.0))
        numerator += weight * abs(float(value))
        denominator += weight
    if denominator <= 0:
        return None
    return numerator / denominator


def _trajectory_state_at(trajectory, target_offset_days: float) -> dict[str, float]:
    if hasattr(trajectory, "iloc") and hasattr(trajectory, "columns"):
        if "time" not in trajectory.columns:
            return {
                "tumor_cells": float(trajectory.iloc[-1]["tumor_cells"]),
                "healthy_cells": float(trajectory.iloc[-1]["healthy_cells"]),
            }
        deltas = (trajectory["time"] - float(target_offset_days)).abs()
        index = int(deltas.idxmin())
        row = trajectory.loc[index]
        return {
            "tumor_cells": float(row["tumor_cells"]),
            "healthy_cells": float(row["healthy_cells"]),
        }

    if isinstance(trajectory, dict):
        if isinstance(trajectory.get("time"), list):
            time_points = trajectory.get("time") or []
            if not time_points:
                raise ValueError("Trajectory time series is empty.")
            nearest_index = min(
                range(len(time_points)),
                key=lambda idx: abs(float(time_points[idx]) - float(target_offset_days)),
            )
            return {
                "tumor_cells": float((trajectory.get("tumor_cells") or [])[nearest_index]),
                "healthy_cells": float((trajectory.get("healthy_cells") or [])[nearest_index]),
            }
        if "tumor_cells" in trajectory and "healthy_cells" in trajectory:
            return {
                "tumor_cells": float(trajectory["tumor_cells"]),
                "healthy_cells": float(trajectory["healthy_cells"]),
            }

    raise ValueError("Unsupported trajectory payload for biomarker prediction.")


def _maybe_float(value: Any) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _decorate_prediction_payload(
    predicted: dict[str, float | None],
    prediction_day: float,
    *,
    model_version: str | None,
    config_hash: str | None,
) -> dict[str, Any]:
    return {
        "m_protein_g_dl": predicted.get("m_protein_g_dl"),
        "flc_ratio": predicted.get("flc_ratio"),
        "hemoglobin_g_dl": predicted.get("hemoglobin_g_dl"),
        "prediction_day": float(prediction_day),
        "source": "observation_model",
        "model_version": model_version,
        "config_hash": config_hash,
    }
