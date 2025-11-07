from __future__ import annotations

import logging
import math
from functools import lru_cache
from pathlib import Path
from typing import Callable, Dict, Tuple

import yaml

logger = logging.getLogger(__name__)

ScheduleFn = Callable[[float], float]


def resolve(
    drug_doses: Dict[str, float],
    time_horizon: float,
    fallback_pk: Dict[str, Dict[str, float]] | None = None,
    fallback_pd: Dict[str, Dict[str, float]] | None = None,
) -> Tuple[Dict[str, Dict[str, float]], Dict[str, Dict[str, float]], Dict[str, ScheduleFn]]:
    """
    Build PK/PD dictionaries and dosing schedules for the requested drugs.

    Falls back to provided defaults when YAML presets are missing.
    """
    registry = _load_registry()
    pk_params: Dict[str, Dict[str, float]] = {}
    pd_params: Dict[str, Dict[str, float]] = {}
    dose_functions: Dict[str, ScheduleFn] = {}
    for drug, dose in drug_doses.items():
        profile = registry.get(drug.lower())
        if profile:
            pk_params[drug] = dict(profile["pk"])
            pd_params[drug] = dict(profile["pd"])
            schedule_cfg = profile.get("schedule", {"type": "continuous"})
            fn = _build_schedule(schedule_cfg, float(dose), time_horizon)
            if fn is not None:
                dose_functions[drug] = fn
        else:
            if fallback_pk and drug in fallback_pk:
                pk_params[drug] = dict(fallback_pk[drug])
            if fallback_pd and drug in fallback_pd:
                pd_params[drug] = dict(fallback_pd[drug])
    return pk_params, pd_params, dose_functions


def get_drug_profile(drug: str) -> Dict[str, object] | None:
    """Return the normalized profile for a drug, if available."""
    registry = _load_registry()
    return registry.get(drug.lower())


def list_profiles() -> Dict[str, Dict[str, object]]:
    """Expose a copy of the loaded registry for UIs/tests."""
    return {name: data.copy() for name, data in _load_registry().items()}


@lru_cache
def _load_registry() -> Dict[str, Dict[str, object]]:
    base_dir = Path(__file__).resolve().parent.parent / "presets" / "drugs"
    registry: Dict[str, Dict[str, object]] = {}
    if not base_dir.exists():
        return registry
    for config_path in sorted(base_dir.glob("*.yaml")):
        try:
            with config_path.open("r", encoding="utf-8") as handle:
                raw = yaml.safe_load(handle) or {}
        except yaml.YAMLError as exc:  # pragma: no cover - defensive for malformed YAML
            logger.warning("Drug preset %s failed to load: %s", config_path.name, exc)
            continue
        try:
            validated = _validate_profile(raw, source=config_path.name)
        except ValueError as exc:
            logger.warning("Drug preset %s skipped: %s", config_path.name, exc)
            continue
        registry[validated["name"]] = validated
    return registry


def _validate_profile(data: Dict[str, object], *, source: str) -> Dict[str, object]:
    if not isinstance(data, dict):
        raise ValueError("configuration must be a mapping")
    name = str(data.get("name", "")).strip().lower()
    if not name:
        raise ValueError("missing drug name")
    pk = data.get("pk") or {}
    pd = data.get("pd") or {}
    if not {"half_life", "Vd"} <= pk.keys():
        raise ValueError("pk requires half_life and Vd")
    if not {"Emax", "EC50"} <= pd.keys():
        raise ValueError("pd requires Emax and EC50")
    dose_range = data.get("dose_range") or {}
    range_min = float(dose_range.get("min", 0.0))
    range_max = float(dose_range.get("max", 0.0))
    if range_min <= 0 or range_max <= 0 or range_min > range_max:
        raise ValueError("dose_range must define positive min/max (min<=max)")
    schedule = _normalize_schedule(data.get("schedule") or {})
    profile = {
        "name": name,
        "display_name": data.get("display_name", name.title()),
        "pk": {"half_life": float(pk["half_life"]), "Vd": float(pk["Vd"])},
        "pd": {"Emax": float(pd["Emax"]), "EC50": float(pd["EC50"])},
        "dose_range": {"min": range_min, "max": range_max},
        "unit": data.get("unit", "mg"),
        "schedule": schedule,
        "schema": data.get("schema", schedule.get("type", "continuous")),
    }
    return profile


def _normalize_schedule(schedule: Dict[str, object]) -> Dict[str, object]:
    schedule_type = str(schedule.get("type", "continuous")).lower()
    normalized: Dict[str, object] = {"type": schedule_type}
    if schedule_type == "cycle":
        normalized["cycle_length"] = float(schedule.get("cycle_length", 28.0))
        normalized["days_on"] = float(schedule.get("days_on", normalized["cycle_length"]))
        normalized["administration_window_days"] = float(schedule.get("administration_window_days", 1.0))
    elif schedule_type == "weekly":
        normalized["days"] = [int(day) % 7 for day in schedule.get("days", [0, 3])]
        normalized["administration_window_days"] = float(schedule.get("administration_window_days", 0.5))
    elif schedule_type == "interval":
        normalized["interval_days"] = float(schedule.get("interval_days", 28.0))
        normalized["administration_window_days"] = float(schedule.get("administration_window_days", 1.0))
    elif schedule_type == "pulsed":
        normalized["days"] = [int(day) for day in schedule.get("days", [])]
        normalized["administration_window_days"] = float(schedule.get("administration_window_days", 1.0))
    else:
        normalized["type"] = "continuous"
        normalized["administration_window_days"] = float(schedule.get("administration_window_days", 1.0))
    return normalized


def _build_schedule(config: Dict[str, object], dose: float, time_horizon: float) -> ScheduleFn | None:
    schedule_type = config.get("type", "continuous")
    window = max(float(config.get("administration_window_days", 1.0)), 1e-3)

    if schedule_type == "cycle":
        cycle = max(float(config.get("cycle_length", time_horizon or 28.0)), 1e-3)
        days_on = min(max(float(config.get("days_on", cycle)), 0.0), cycle)

        def _cycle_fn(t: float) -> float:
            position = t % cycle
            return dose / window if position < days_on else 0.0

        return _cycle_fn

    if schedule_type == "weekly":
        days = [int(day) % 7 for day in config.get("days", [1, 4])]

        def _weekly_fn(t: float) -> float:
            day_index = int(math.floor(t)) % 7
            fractional = t - math.floor(t)
            if day_index in days and fractional < window:
                return dose / window
            return 0.0

        return _weekly_fn

    if schedule_type == "interval":
        interval = max(float(config.get("interval_days", 28.0)), 1e-3)

        def _interval_fn(t: float) -> float:
            if (t % interval) < window:
                return dose / window
            return 0.0

        return _interval_fn

    if schedule_type == "pulsed":
        days = [int(day) for day in config.get("days", [])]

        def _pulsed_fn(t: float) -> float:
            current_day = int(math.floor(t))
            if current_day in days and (t - current_day) < window:
                return dose / window
            return 0.0

        return _pulsed_fn

    if schedule_type == "continuous":
        def _continuous_fn(_t: float) -> float:
            return dose / max(time_horizon, 1e-6)

        return _continuous_fn

    return None
