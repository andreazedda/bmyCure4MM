from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path
from typing import Any, Dict

import yaml

from clinic.models import Assessment


def build_patient_twin(assessment: Assessment) -> Dict[str, float]:
    """
    Build a lightweight patient twin inferred from laboratory values.

    Risk score drives downstream PK/PD modifiers used by the simulator.
    """
    config = _load_twin_config()
    weights = config.get("weights", {})
    total_weight = float(sum(weights.values()) or 1.0)

    r_iss_component = _score_riss(assessment, config)
    ldh_component = _score_ldh(assessment, config)
    beta_component = _score_beta2m(assessment, config)
    flc_component = _score_flc_ratio(assessment, config)

    weighted_sum = (
        weights.get("riss", 0.0) * r_iss_component
        + weights.get("ldh", 0.0) * ldh_component
        + weights.get("beta2m", 0.0) * beta_component
        + weights.get("flc_ratio", 0.0) * flc_component
    )
    risk_score = max(0.0, min(1.0, weighted_sum / total_weight))

    growth_cfg = config.get("growth_mapping", {})
    carrying_cfg = config.get("carrying_capacity", {})
    immune_cfg = config.get("immune_compromise_index", {})

    twin_params = {
        "risk_score": risk_score,
        "tumor_growth_rate": _lerp(growth_cfg.get("tumor_growth_rate", {}), risk_score),
        "healthy_growth_rate": _lerp(growth_cfg.get("healthy_growth_rate", {}), 1.0 - risk_score),
        "carrying_capacity_tumor": _lerp(carrying_cfg.get("tumor", {}), risk_score),
        "carrying_capacity_healthy": _lerp(carrying_cfg.get("healthy", {}), 1.0 - risk_score),
        "immune_compromise_index": _lerp(immune_cfg, risk_score, default=1.0),
    }
    return twin_params


@lru_cache
def _load_twin_config() -> Dict[str, Any]:
    presets_dir = Path(__file__).resolve().parent / "presets"
    config_path = presets_dir / "twin_risk.yaml"
    if not config_path.exists():
        raise FileNotFoundError(f"Twin risk configuration missing: {config_path}")
    with config_path.open("r", encoding="utf-8") as handle:
        data = yaml.safe_load(handle) or {}
    return json.loads(json.dumps(data))  # ensure plain Python types


def _score_riss(assessment: Assessment, config: Dict[str, Any]) -> float:
    riss = (assessment.r_iss or "").upper()
    mapping = config.get("r_iss_map", {"I": 0.2, "II": 0.6, "III": 1.0})
    return float(mapping.get(riss, mapping.get("II", 0.6)))


def _score_ldh(assessment: Assessment, config: Dict[str, Any]) -> float:
    value = _maybe_float(assessment.ldH_u_l)
    if value is None:
        return 0.5
    cfg = config.get("ldh_breakpoints", {"normal_upper": 250, "high": 500})
    normal = float(cfg.get("normal_upper", 250))
    high = float(cfg.get("high", 500))
    if value <= normal:
        return 0.1
    if value >= high:
        return 1.0
    span = max(high - normal, 1e-3)
    return 0.1 + 0.9 * ((value - normal) / span)


def _score_beta2m(assessment: Assessment, config: Dict[str, Any]) -> float:
    value = _maybe_float(assessment.beta2m_mg_l)
    if value is None:
        return 0.5
    cfg = config.get("beta2m_range", {"min": 2.0, "max": 12.0})
    lo = float(cfg.get("min", 2.0))
    hi = float(cfg.get("max", 12.0))
    if value <= lo:
        return 0.1
    if value >= hi:
        return 1.0
    return (value - lo) / max(hi - lo, 1e-6)


def _score_flc_ratio(assessment: Assessment, config: Dict[str, Any]) -> float:
    value = _maybe_float(assessment.flc_ratio)
    if value is None:
        return 0.4
    cfg = config.get(
        "flc_ratio_range",
        {"normal_low": 0.26, "normal_high": 1.65, "max": 20.0, "min": 0.05},
    )
    low = float(cfg.get("normal_low", 0.26))
    high = float(cfg.get("normal_high", 1.65))
    extreme_high = float(cfg.get("max", 20.0))
    extreme_low = float(cfg.get("min", 0.05))
    if low <= value <= high:
        return 0.1
    if value > high:
        span = max(extreme_high - high, 1e-6)
        return min(1.0, 0.1 + 0.9 * ((value - high) / span))
    # value below low
    span = max(low - extreme_low, 1e-6)
    return min(1.0, 0.1 + 0.9 * ((low - value) / span))


def _lerp(config: Dict[str, Any], weight: float, default: float | None = None) -> float:
    if not config:
        if default is None:
            return float(weight)
        return default
    low = float(config.get("min", default if default is not None else 0.0))
    high = float(config.get("max", default if default is not None else 1.0))
    weight = max(0.0, min(1.0, weight))
    return low + (high - low) * weight


def _maybe_float(value: Any) -> float | None:
    if value is None or value == "":
        return None
    try:
        return float(value)
    except (TypeError, ValueError):
        return None
