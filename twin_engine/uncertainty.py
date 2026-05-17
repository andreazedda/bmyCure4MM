from __future__ import annotations

import copy
import math
import time
from dataclasses import dataclass
from datetime import timedelta
from typing import Any

import numpy as np

from .calibration import DEFAULT_PARAMETER_BOUNDS
from .counterfactual import _determine_toxicity_constraint_penalty, _normalize_intervention_definition
from .simulation_bridge import build_solver_inputs_from_twin_state, run_patient_simulation
from .therapy_schedule import build_therapy_schedule
from .toxicity_dynamics import DEFAULT_TOXICITY_COEFFICIENTS, compute_toxicity_dynamics
from .toxicity_model import compute_toxicity_constraints


HEURISTIC_PARAMETER_UNCERTAINTY_SOURCE = "heuristic_perturbation_not_calibrated_distribution"
CALIBRATED_PARAMETER_UNCERTAINTY_SOURCE = "calibrated_parameter_standard_errors"

COMPLETED_STATUS = "completed"
INSUFFICIENT_PARAMETER_UNCERTAINTY_STATUS = "insufficient_parameter_uncertainty"
SIMULATION_FAILED_STATUS = "simulation_failed"


@dataclass(frozen=True)
class UncertaintyConfig:
    n_samples: int = 100
    random_seed: int = 17
    parameter_noise_scale: float = 0.10
    observation_noise_scale: float = 0.05
    toxicity_coefficient_noise_scale: float = 0.10
    max_runtime_seconds: float | None = None


def build_parameter_uncertainty_space(twin_state) -> dict[str, Any]:
    diagnostics = dict(getattr(twin_state, "parameter_uncertainty", {}) or {})
    flat_parameters = _flatten_numeric_parameters(dict(getattr(twin_state, "parameters", {}) or {}))
    calibrated_scales = dict(diagnostics.get("parameter_standard_errors") or {})
    source = CALIBRATED_PARAMETER_UNCERTAINTY_SOURCE if calibrated_scales else HEURISTIC_PARAMETER_UNCERTAINTY_SOURCE

    parameters = {}
    for name, baseline in sorted(flat_parameters.items()):
        scale = _parameter_standard_error(name, baseline, calibrated_scales)
        lower_bound, upper_bound = _parameter_bounds(name, baseline)
        parameters[name] = {
            "baseline": float(baseline),
            "sampling_std": float(scale),
            "lower_bound": lower_bound,
            "upper_bound": upper_bound,
            "sampling_distribution": "truncated_normal",
        }

    toxicity_coefficients = {
        name: {
            "baseline": float(value),
            "sampling_std": max(abs(float(value)) * 0.10, 0.01),
            "lower_bound": 0.0,
            "upper_bound": 1.5,
        }
        for name, value in sorted(DEFAULT_TOXICITY_COEFFICIENTS.items())
    }
    status = COMPLETED_STATUS if parameters else INSUFFICIENT_PARAMETER_UNCERTAINTY_STATUS
    return {
        "status": status,
        "parameter_uncertainty_source": source,
        "n_parameters": len(parameters),
        "parameters": parameters,
        "toxicity_coefficients": toxicity_coefficients,
    }


def sample_parameter_sets(twin_state, config: UncertaintyConfig) -> dict[str, Any]:
    space = build_parameter_uncertainty_space(twin_state)
    if space["status"] != COMPLETED_STATUS:
        return {
            "status": INSUFFICIENT_PARAMETER_UNCERTAINTY_STATUS,
            "parameter_uncertainty_source": space["parameter_uncertainty_source"],
            "n_samples": 0,
            "seed": int(config.random_seed),
            "samples": [],
        }

    rng = np.random.default_rng(int(config.random_seed))
    samples = []
    for sample_index in range(max(int(config.n_samples), 0)):
        flat_sample = {}
        for name, item in space["parameters"].items():
            std = _noise_scale_for_parameter(name, item["sampling_std"], config)
            value = rng.normal(loc=float(item["baseline"]), scale=std)
            flat_sample[name] = _clip_sample(value, item["lower_bound"], item["upper_bound"], item["baseline"])
        toxicity_coefficients = {}
        for name, item in space["toxicity_coefficients"].items():
            std = max(abs(float(item["baseline"])) * float(config.toxicity_coefficient_noise_scale), 0.01)
            toxicity_coefficients[name] = _clip_sample(
                rng.normal(loc=float(item["baseline"]), scale=std),
                item["lower_bound"],
                item["upper_bound"],
                item["baseline"],
            )
        samples.append(
            {
                "sample_index": int(sample_index),
                "parameters": _restore_flat_parameters(flat_sample),
                "toxicity_coefficients": toxicity_coefficients,
                "parameter_uncertainty_source": space["parameter_uncertainty_source"],
            }
        )

    return {
        "status": COMPLETED_STATUS,
        "parameter_uncertainty_source": space["parameter_uncertainty_source"],
        "n_samples": len(samples),
        "seed": int(config.random_seed),
        "samples": samples,
    }


def summarize_distribution(values) -> dict[str, Any]:
    cleaned = _clean_numeric_values(values)
    if not cleaned:
        return {
            "count": 0,
            "mean": None,
            "median": None,
            "std": None,
            "min": None,
            "max": None,
        }
    array = np.asarray(cleaned, dtype=float)
    return {
        "count": int(array.size),
        "mean": float(np.mean(array)),
        "median": float(np.median(array)),
        "std": float(np.std(array, ddof=0)),
        "min": float(np.min(array)),
        "max": float(np.max(array)),
    }


def compute_interval(values, lower=0.05, upper=0.95) -> dict[str, Any]:
    cleaned = _clean_numeric_values(values)
    if not cleaned:
        return {
            "lower": None,
            "upper": None,
            "width": None,
            "lower_quantile": float(lower),
            "upper_quantile": float(upper),
        }
    array = np.asarray(cleaned, dtype=float)
    lower_value, upper_value = np.quantile(array, [float(lower), float(upper)])
    return {
        "lower": float(lower_value),
        "upper": float(upper_value),
        "width": float(upper_value - lower_value),
        "lower_quantile": float(lower),
        "upper_quantile": float(upper),
    }


def classify_uncertainty_width(metric_name, interval_width, reference_scale) -> str:
    if interval_width is None or reference_scale is None:
        return "unavailable"
    scale = max(abs(float(reference_scale)), _reference_floor(metric_name))
    if scale <= 0.0:
        return "unavailable"
    ratio = abs(float(interval_width)) / scale
    if ratio <= 0.15:
        return "narrow"
    if ratio <= 0.40:
        return "moderate"
    return "wide"


def run_counterfactual_uncertainty(patient, twin_state, intervention_definition, horizon_days, config: UncertaintyConfig) -> dict[str, Any]:
    horizon = int(horizon_days)
    sampled_sets = sample_parameter_sets(twin_state, config)
    if sampled_sets["status"] != COMPLETED_STATUS:
        return {
            "status": INSUFFICIENT_PARAMETER_UNCERTAINTY_STATUS,
            "parameter_uncertainty_source": sampled_sets["parameter_uncertainty_source"],
            "requested_samples": int(config.n_samples),
            "n_completed_samples": 0,
            "seed": int(config.random_seed),
            "metric_summaries": {},
            "samples": [],
        }

    execution_definition = _normalize_intervention_definition(
        intervention_definition,
        twin_state.state_date,
        horizon,
    )
    start_date = twin_state.state_date
    end_date = start_date + timedelta(days=max(horizon - 1, 0))
    baseline_schedule = build_therapy_schedule(patient, start_date, end_date)
    baseline_solver_inputs = build_solver_inputs_from_twin_state(twin_state, baseline_schedule, horizon)
    toxicity_constraints = compute_toxicity_constraints(patient)

    point_metrics = _run_uncertainty_sample(
        patient=patient,
        twin_state=twin_state,
        execution_definition=execution_definition,
        horizon_days=horizon,
        baseline_solver_inputs=baseline_solver_inputs,
        toxicity_constraints=toxicity_constraints,
        sampled_parameters=dict(getattr(twin_state, "parameters", {}) or {}),
        sampled_toxicity_coefficients=None,
    )

    started_at = time.monotonic()
    completed_samples = []
    failed_samples = []
    for item in sampled_sets["samples"]:
        if config.max_runtime_seconds is not None and (time.monotonic() - started_at) >= float(config.max_runtime_seconds):
            break
        try:
            metrics = _run_uncertainty_sample(
                patient=patient,
                twin_state=twin_state,
                execution_definition=execution_definition,
                horizon_days=horizon,
                baseline_solver_inputs=baseline_solver_inputs,
                toxicity_constraints=toxicity_constraints,
                sampled_parameters=item["parameters"],
                sampled_toxicity_coefficients=item["toxicity_coefficients"],
            )
        except Exception as exc:
            failed_samples.append({"sample_index": item["sample_index"], "error": str(exc)})
            continue
        completed_samples.append(
            {
                "sample_index": item["sample_index"],
                "metrics": metrics,
            }
        )

    if not completed_samples:
        return {
            "status": SIMULATION_FAILED_STATUS,
            "parameter_uncertainty_source": sampled_sets["parameter_uncertainty_source"],
            "requested_samples": int(config.n_samples),
            "n_completed_samples": 0,
            "seed": int(config.random_seed),
            "metric_summaries": {},
            "samples": [],
            "failed_samples": failed_samples,
        }

    metric_summaries = {}
    metric_names = sorted(
        set(point_metrics.keys())
        | {
            metric_name
            for sample in completed_samples
            for metric_name in sample["metrics"].keys()
        }
    )

    for metric_name in metric_names:
        values = [sample["metrics"].get(metric_name) for sample in completed_samples]
        distribution = summarize_distribution(values)
        p05_p95 = compute_interval(values, lower=0.05, upper=0.95)
        p25_p75 = compute_interval(values, lower=0.25, upper=0.75)
        mean_value = distribution["mean"]
        std_value = distribution["std"]
        coefficient_of_variation = None
        if mean_value is not None and std_value is not None and abs(float(mean_value)) > 1.0e-9:
            coefficient_of_variation = abs(float(std_value) / float(mean_value))
        point_estimate = point_metrics.get(metric_name)
        reference_scale = _reference_scale(metric_name, point_estimate, distribution)
        metric_summaries[metric_name] = {
            "metric": metric_name,
            "point_estimate": point_estimate,
            "median": distribution["median"],
            "mean": mean_value,
            "p05": p05_p95["lower"],
            "p25": p25_p75["lower"],
            "p75": p25_p75["upper"],
            "p95": p05_p95["upper"],
            "interval_width": p05_p95["width"],
            "coefficient_of_variation": coefficient_of_variation,
            "uncertainty_classification": classify_uncertainty_width(metric_name, p05_p95["width"], reference_scale),
            "n_samples": len(completed_samples),
            "seed": int(config.random_seed),
            "status": COMPLETED_STATUS,
        }

    return {
        "status": COMPLETED_STATUS,
        "parameter_uncertainty_source": sampled_sets["parameter_uncertainty_source"],
        "requested_samples": int(config.n_samples),
        "n_completed_samples": len(completed_samples),
        "seed": int(config.random_seed),
        "metric_summaries": metric_summaries,
        "samples": completed_samples,
        "failed_samples": failed_samples,
        "runtime_seconds": float(time.monotonic() - started_at),
    }


def _run_uncertainty_sample(
    *,
    patient,
    twin_state,
    execution_definition: dict[str, Any],
    horizon_days: int,
    baseline_solver_inputs: dict[str, Any],
    toxicity_constraints: dict[str, Any],
    sampled_parameters: dict[str, Any],
    sampled_toxicity_coefficients: dict[str, float] | None,
) -> dict[str, float | None]:
    alternative_overrides = _merge_parameters(sampled_parameters, dict(execution_definition.get("parameter_overrides", {}) or {}))
    if execution_definition.get("drug_doses"):
        alternative_overrides["drug_doses"] = execution_definition["drug_doses"]
    if execution_definition.get("random_seed") is not None:
        alternative_overrides["seed"] = int(execution_definition["random_seed"])

    alternative_result = run_patient_simulation(
        twin_state,
        therapy_schedule=execution_definition.get("schedule"),
        horizon_days=int(horizon_days),
        overrides=alternative_overrides,
    )
    toxicity_dynamics = compute_toxicity_dynamics(
        patient,
        alternative_result["summary"].get("exposure_profiles") or {},
        coefficient_overrides=sampled_toxicity_coefficients,
    )
    penalty = _determine_toxicity_constraint_penalty(
        baseline_solver_inputs.get("raw_parameters") or {},
        alternative_result["solver_inputs"].get("raw_parameters") or {},
        toxicity_constraints or {},
    )
    utility_v2 = _research_utility_v2(
        summary=alternative_result["summary"],
        toxicity_penalty=penalty,
        toxicity_dynamics=toxicity_dynamics,
    )
    predicted = alternative_result["summary"].get("predicted_biomarkers") or {}
    return {
        "tumor_reduction": _safe_float(alternative_result["summary"].get("tumor_reduction")),
        "healthy_loss": _safe_float(alternative_result["summary"].get("healthy_loss")),
        "durability_index": _safe_float(alternative_result["summary"].get("durability_index")),
        "liver_toxicity_signal_0_1": _safe_float(toxicity_dynamics.get("liver_toxicity_signal_0_1")),
        "neutropenia_signal_0_1": _safe_float(toxicity_dynamics.get("neutropenia_signal_0_1")),
        "research_utility": _research_utility(
            summary=alternative_result["summary"],
            toxicity_penalty=penalty,
        ),
        "research_utility_v2": utility_v2,
        "predicted_m_protein_g_dl": _safe_float(predicted.get("m_protein_g_dl")),
        "predicted_flc_ratio": _safe_float(predicted.get("flc_ratio")),
        "predicted_hemoglobin_g_dl": _safe_float(predicted.get("hemoglobin_g_dl")),
    }


def _research_utility(*, summary: dict[str, Any], toxicity_penalty: float) -> float:
    disease_control_score = float(summary.get("tumor_reduction") or 0.0)
    healthy_preservation_score = 1.0 - float(summary.get("healthy_loss") or 0.0)
    durability_score = float(summary.get("durability_index") or 0.0)
    return disease_control_score + healthy_preservation_score + durability_score - float(toxicity_penalty)


def _research_utility_v2(*, summary: dict[str, Any], toxicity_penalty: float, toxicity_dynamics: dict[str, Any]) -> float:
    base = _research_utility(summary=summary, toxicity_penalty=toxicity_penalty)
    liver_signal = float(toxicity_dynamics.get("liver_toxicity_signal_0_1") or 0.0)
    neutropenia_signal = float(toxicity_dynamics.get("neutropenia_signal_0_1") or 0.0)
    return base - 0.5 * liver_signal - 0.5 * neutropenia_signal


def _flatten_numeric_parameters(payload: dict[str, Any], *, prefix: str = "") -> dict[str, float]:
    flattened = {}
    for key, value in sorted((payload or {}).items()):
        name = f"{prefix}.{key}" if prefix else str(key)
        if isinstance(value, dict):
            flattened.update(_flatten_numeric_parameters(value, prefix=name))
            continue
        if isinstance(value, bool):
            continue
        if isinstance(value, (int, float)):
            flattened[name] = float(value)
    return flattened


def _restore_flat_parameters(flattened: dict[str, float]) -> dict[str, Any]:
    restored: dict[str, Any] = {}
    for name, value in flattened.items():
        cursor = restored
        parts = str(name).split(".")
        for part in parts[:-1]:
            cursor = cursor.setdefault(part, {})
        cursor[parts[-1]] = float(value)
    return restored


def _parameter_standard_error(name: str, baseline: float, calibrated_scales: dict[str, Any]) -> float:
    lookup_names = [name, name.replace(".", "__")]
    for lookup_name in lookup_names:
        if calibrated_scales.get(lookup_name) is not None:
            try:
                return max(abs(float(calibrated_scales[lookup_name])), 1.0e-6)
            except (TypeError, ValueError):
                break
    return max(abs(float(baseline)) * 0.10, _reference_floor(name) * 0.05)


def _parameter_bounds(name: str, baseline: float) -> tuple[float | None, float | None]:
    lookup_name = name.replace(".", "__")
    if lookup_name in DEFAULT_PARAMETER_BOUNDS:
        lower, upper = DEFAULT_PARAMETER_BOUNDS[lookup_name]
        return float(lower), float(upper)
    if baseline >= 0.0:
        return 0.0, None
    return None, None


def _noise_scale_for_parameter(name: str, baseline_scale: float, config: UncertaintyConfig) -> float:
    if name.startswith("observation."):
        multiplier = float(config.observation_noise_scale) / 0.05 if float(config.observation_noise_scale) > 0 else 0.0
    else:
        multiplier = float(config.parameter_noise_scale) / 0.10 if float(config.parameter_noise_scale) > 0 else 0.0
    if multiplier <= 0.0:
        return 0.0
    return max(float(baseline_scale) * multiplier, 1.0e-6)


def _clip_sample(value: float, lower: float | None, upper: float | None, baseline: float) -> float:
    sampled = float(value)
    if lower is not None:
        sampled = max(sampled, float(lower))
    if upper is not None:
        sampled = min(sampled, float(upper))
    if not math.isfinite(sampled):
        return float(baseline)
    return sampled


def _merge_parameters(base: dict[str, Any], updates: dict[str, Any]) -> dict[str, Any]:
    merged = copy.deepcopy(dict(base or {}))
    for key, value in (updates or {}).items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _merge_parameters(merged.get(key) or {}, value)
        else:
            merged[key] = copy.deepcopy(value)
    return merged


def _clean_numeric_values(values) -> list[float]:
    cleaned = []
    for value in values or []:
        try:
            numeric = float(value)
        except (TypeError, ValueError):
            continue
        if math.isfinite(numeric):
            cleaned.append(numeric)
    return cleaned


def _reference_scale(metric_name: str, point_estimate: float | None, distribution: dict[str, Any]) -> float:
    if point_estimate is not None:
        return max(abs(float(point_estimate)), _reference_floor(metric_name))
    median = distribution.get("median")
    if median is not None:
        return max(abs(float(median)), _reference_floor(metric_name))
    return _reference_floor(metric_name)


def _reference_floor(metric_name: str) -> float:
    name = str(metric_name).lower()
    if name in {"healthy_loss", "durability_index", "liver_toxicity_signal_0_1", "neutropenia_signal_0_1"}:
        return 1.0
    if "utility" in name:
        return 1.0
    if "predicted_" in name:
        return 0.5
    return 0.1


def _safe_float(value) -> float | None:
    if value is None:
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None