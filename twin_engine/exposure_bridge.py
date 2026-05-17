from __future__ import annotations

import hashlib
import json
import math
from dataclasses import dataclass
from datetime import date, timedelta
from typing import Any, Callable


class ExposureScheduleError(ValueError):
    pass


@dataclass(frozen=True)
class ExposureProfile:
    horizon_days: int
    time_grid_days: list[int]
    drug: str
    dose_events: list[dict[str, Any]]
    daily_administered_dose_mg: list[float]
    cumulative_dose_mg: list[float]
    average_daily_dose_mg: float
    peak_administered_dose_mg: float
    dosing_frequency_days: float | None
    interruption_days: list[int]
    schedule_type: str
    schedule_label: str
    exposure_resolution_status: str
    total_cumulative_dose_mg: float
    interruption_fraction: float
    exposure_profile_hash: str

    def to_dict(self) -> dict[str, Any]:
        return {
            "horizon_days": self.horizon_days,
            "time_grid_days": list(self.time_grid_days),
            "drug": self.drug,
            "dose_events": list(self.dose_events),
            "daily_administered_dose_mg": list(self.daily_administered_dose_mg),
            "cumulative_dose_mg": list(self.cumulative_dose_mg),
            "total_cumulative_dose_mg": self.total_cumulative_dose_mg,
            "average_daily_dose_mg": self.average_daily_dose_mg,
            "peak_administered_dose_mg": self.peak_administered_dose_mg,
            "dosing_frequency_days": self.dosing_frequency_days,
            "interruption_days": list(self.interruption_days),
            "interruption_fraction": self.interruption_fraction,
            "schedule_type": self.schedule_type,
            "schedule_label": self.schedule_label,
            "exposure_resolution_status": self.exposure_resolution_status,
            "exposure_profile_hash": self.exposure_profile_hash,
        }


def build_exposure_profile(
    *,
    drug: str,
    dose_mg: float | int | None,
    horizon_days: int,
    schedule: dict[str, Any] | None = None,
    schedule_label: str = "research_intervention",
    phases: list[dict[str, Any]] | None = None,
) -> ExposureProfile:
    if int(horizon_days) <= 0:
        raise ExposureScheduleError("horizon_days must be positive")
    canonical_drug = str(drug or "").strip().lower()
    if not canonical_drug:
        raise ExposureScheduleError("drug is required")
    base_date = date(2000, 1, 1)
    if phases:
        items = []
        for phase in phases:
            start_day = int(phase.get("start_day", 0) or 0)
            duration = int(phase.get("duration_days", horizon_days - start_day) or 0)
            phase_schedule = phase.get("schedule") or schedule or {"type": "daily"}
            phase_dose = phase.get("dose_mg", dose_mg)
            item_start = base_date + timedelta(days=start_day)
            item_end = item_start + timedelta(days=max(duration - 1, 0))
            items.append(_schedule_item(canonical_drug, phase_dose, phase_schedule, item_start, item_end, schedule_label))
        schedule_payload = {
            "horizon_start_date": base_date.isoformat(),
            "start_date": base_date.isoformat(),
            "end_date": (base_date + timedelta(days=int(horizon_days) - 1)).isoformat(),
            "items": items,
            "interruptions": [],
        }
        return build_exposure_profiles_from_schedule(schedule_payload, int(horizon_days))[canonical_drug]

    item = _schedule_item(
        canonical_drug,
        dose_mg,
        schedule or {"type": "daily"},
        base_date,
        base_date + timedelta(days=int(horizon_days) - 1),
        schedule_label,
    )
    schedule_payload = {
        "horizon_start_date": base_date.isoformat(),
        "start_date": base_date.isoformat(),
        "end_date": (base_date + timedelta(days=int(horizon_days) - 1)).isoformat(),
        "items": [item],
        "interruptions": [],
    }
    return build_exposure_profiles_from_schedule(schedule_payload, int(horizon_days))[canonical_drug]


def build_legacy_scalar_exposure_profile(
    *,
    drug: str,
    scalar_dose_mg: float,
    horizon_days: int,
    schedule_label: str = "legacy_scalar_override",
) -> ExposureProfile:
    horizon = int(horizon_days)
    daily_value = float(scalar_dose_mg) / max(horizon, 1)
    return _finalize_profile(
        horizon_days=horizon,
        drug=str(drug).lower(),
        daily=[daily_value for _ in range(horizon)],
        interruption_days=[],
        schedule_type="legacy_scalar_override",
        schedule_label=schedule_label,
        exposure_resolution_status="legacy_scalar_average_only; temporal profile unavailable",
    )


def build_exposure_profiles_from_schedule(schedule_payload: dict[str, Any] | None, horizon_days: int) -> dict[str, ExposureProfile]:
    if not schedule_payload:
        return {}
    horizon = int(horizon_days)
    if horizon <= 0:
        raise ExposureScheduleError("horizon_days must be positive")

    horizon_start = _parse_date(
        schedule_payload.get("horizon_start_date")
        or schedule_payload.get("simulation_start_date")
        or schedule_payload.get("start_date")
    )
    daily_by_drug: dict[str, list[float]] = {}
    interruptions_by_drug: dict[str, set[int]] = {}
    schedule_types_by_drug: dict[str, set[str]] = {}
    labels_by_drug: dict[str, str] = {}

    for item in schedule_payload.get("items", []) or []:
        start_offset = max((_parse_date(item.get("start_date")) - horizon_start).days, 0)
        end_offset = min((_parse_date(item.get("end_date")) - horizon_start).days, horizon - 1)
        if end_offset < 0 or start_offset >= horizon or end_offset < start_offset:
            continue
        for drug, dose_entry in _dose_entries_for_item(item).items():
            canonical_drug = str(drug or "").strip().lower()
            if not canonical_drug:
                continue
            dose = _extract_numeric_dose(dose_entry)
            if dose is None:
                continue
            daily = daily_by_drug.setdefault(canonical_drug, [0.0 for _ in range(horizon)])
            interruptions = interruptions_by_drug.setdefault(canonical_drug, set())
            schedule_spec = _schedule_spec_for_item(item, dose_entry)
            schedule_type = str(schedule_spec.get("type") or "daily").lower()
            schedule_types_by_drug.setdefault(canonical_drug, set()).add(schedule_type)
            labels_by_drug.setdefault(canonical_drug, str(item.get("regimen_name") or item.get("label") or "research_schedule"))
            if schedule_type == "hold" or float(dose) <= 0.0:
                for day in range(start_offset, end_offset + 1):
                    daily[day] = 0.0
                    interruptions.add(day)
                continue
            admin_days = _administration_days(schedule_spec, start_offset, end_offset, item)
            for day in admin_days:
                daily[day] += float(dose)

    _apply_interruption_days(schedule_payload, horizon_start, horizon, daily_by_drug, interruptions_by_drug)

    profiles = {}
    for drug, daily in daily_by_drug.items():
        schedule_types = schedule_types_by_drug.get(drug, {"unavailable"})
        schedule_type = "phased" if len(schedule_types) > 1 else next(iter(schedule_types))
        profiles[drug] = _finalize_profile(
            horizon_days=horizon,
            drug=drug,
            daily=daily,
            interruption_days=sorted(interruptions_by_drug.get(drug, set())),
            schedule_type=schedule_type,
            schedule_label=labels_by_drug.get(drug, "research_schedule"),
            exposure_resolution_status="time_resolved_profile_available",
        )
    return profiles


def compare_exposure_profiles(
    profile_a: ExposureProfile | dict[str, Any] | None,
    profile_b: ExposureProfile | dict[str, Any] | None,
    *,
    tolerance: float = 1.0e-8,
) -> dict[str, Any]:
    if profile_a is None or profile_b is None:
        return {
            "same_average_exposure": False,
            "different_temporal_profile": False,
            "exposure_profile_distance": None,
            "max_daily_dose_difference": None,
            "interruption_fraction": None,
            "comparison_status": "profile_unavailable",
        }
    payload_a = _profile_payload(profile_a)
    payload_b = _profile_payload(profile_b)
    daily_a = [float(value) for value in payload_a.get("daily_administered_dose_mg", [])]
    daily_b = [float(value) for value in payload_b.get("daily_administered_dose_mg", [])]
    horizon = max(len(daily_a), len(daily_b), 1)
    daily_a = daily_a + [0.0] * (horizon - len(daily_a))
    daily_b = daily_b + [0.0] * (horizon - len(daily_b))
    diffs = [abs(a - b) for a, b in zip(daily_a, daily_b)]
    distance = sum(diffs) / horizon
    same_average = abs(float(payload_a.get("average_daily_dose_mg") or 0.0) - float(payload_b.get("average_daily_dose_mg") or 0.0)) <= tolerance
    return {
        "same_average_exposure": same_average,
        "different_temporal_profile": distance > tolerance,
        "exposure_profile_distance": float(distance),
        "max_daily_dose_difference": float(max(diffs) if diffs else 0.0),
        "interruption_fraction": float(payload_b.get("interruption_fraction") or 0.0),
        "interruption_fraction_difference": abs(float(payload_a.get("interruption_fraction") or 0.0) - float(payload_b.get("interruption_fraction") or 0.0)),
        "comparison_status": "compared",
        "profile_a_hash": payload_a.get("exposure_profile_hash"),
        "profile_b_hash": payload_b.get("exposure_profile_hash"),
    }


def build_exposure_dose_function(profile: ExposureProfile | dict[str, Any]) -> Callable[[float], float]:
    payload = _profile_payload(profile)
    daily = [float(value) for value in payload.get("daily_administered_dose_mg", [])]

    def _dose_at_time(t: float) -> float:
        day_index = int(math.floor(max(float(t), 0.0)))
        if day_index < 0 or day_index >= len(daily):
            return 0.0
        return max(float(daily[day_index]), 0.0)

    return _dose_at_time


def exposure_profiles_to_payload(profiles: dict[str, ExposureProfile]) -> dict[str, dict[str, Any]]:
    return {drug: profile.to_dict() for drug, profile in sorted(profiles.items())}


def exposure_profile_from_payload(payload: dict[str, Any] | None) -> ExposureProfile | None:
    if not payload:
        return None
    return ExposureProfile(
        horizon_days=int(payload.get("horizon_days") or 0),
        time_grid_days=[int(value) for value in payload.get("time_grid_days", [])],
        drug=str(payload.get("drug") or ""),
        dose_events=list(payload.get("dose_events") or []),
        daily_administered_dose_mg=[float(value) for value in payload.get("daily_administered_dose_mg", [])],
        cumulative_dose_mg=[float(value) for value in payload.get("cumulative_dose_mg", [])],
        average_daily_dose_mg=float(payload.get("average_daily_dose_mg") or 0.0),
        peak_administered_dose_mg=float(payload.get("peak_administered_dose_mg") or 0.0),
        dosing_frequency_days=payload.get("dosing_frequency_days"),
        interruption_days=[int(value) for value in payload.get("interruption_days", [])],
        schedule_type=str(payload.get("schedule_type") or "unavailable"),
        schedule_label=str(payload.get("schedule_label") or ""),
        exposure_resolution_status=str(payload.get("exposure_resolution_status") or "unavailable"),
        total_cumulative_dose_mg=float(payload.get("total_cumulative_dose_mg") or 0.0),
        interruption_fraction=float(payload.get("interruption_fraction") or 0.0),
        exposure_profile_hash=str(payload.get("exposure_profile_hash") or ""),
    )


def _schedule_item(drug: str, dose_mg: float | int | None, schedule: dict[str, Any], start: date, end: date, label: str) -> dict[str, Any]:
    return {
        "regimen_name": label,
        "start_date": start.isoformat(),
        "end_date": end.isoformat(),
        "cycle_length_days": schedule.get("cycle_length_days") or schedule.get("cycle_length") or schedule.get("interval_days") or schedule.get("every_days"),
        "days_on": schedule.get("days_on"),
        "components": [{"drug": drug, "supported_by_solver": True}],
        "doses": {drug: {"dose": float(dose_mg or 0.0), "unit": "mg", "schedule": dict(schedule)}},
    }


def _dose_entries_for_item(item: dict[str, Any]) -> dict[str, Any]:
    doses = item.get("doses") or {}
    if doses:
        return dict(doses)
    entries = {}
    for component in item.get("components", []) or []:
        drug = component.get("drug")
        if drug:
            entries[str(drug)] = {"dose": 0.0, "schedule": {"type": "hold"}}
    return entries


def _schedule_spec_for_item(item: dict[str, Any], dose_entry: Any) -> dict[str, Any]:
    if isinstance(dose_entry, dict) and isinstance(dose_entry.get("schedule"), dict):
        raw = dict(dose_entry.get("schedule") or {})
    else:
        raw = {}
    if not raw:
        raw = dict(item.get("schedule") or {}) if isinstance(item.get("schedule"), dict) else {}
    if not raw:
        raw = _derive_schedule_from_cycle_fields(item)
    schedule_type = str(raw.get("type") or "daily").lower()
    if schedule_type == "alternate_day":
        raw = {"type": "interval", "every_days": 2}
    elif schedule_type in {"every_n_days", "every"}:
        raw["type"] = "interval"
    elif schedule_type not in {"daily", "interval", "cycle", "hold"}:
        raise ExposureScheduleError(f"Unsupported exposure schedule type: {schedule_type}")
    return raw


def _derive_schedule_from_cycle_fields(item: dict[str, Any]) -> dict[str, Any]:
    cycle_length = item.get("cycle_length_days")
    days_on = item.get("days_on") or []
    if cycle_length and days_on:
        cycle_length_int = int(cycle_length)
        days = [int(day) for day in days_on]
        if cycle_length_int == 1 and days == [1]:
            return {"type": "daily"}
        if len(days) == 1 and days[0] == 1:
            return {"type": "interval", "every_days": cycle_length_int}
        return {"type": "cycle", "cycle_length_days": cycle_length_int, "days_on": days}
    return {"type": "daily"}


def _administration_days(schedule_spec: dict[str, Any], start_offset: int, end_offset: int, item: dict[str, Any]) -> list[int]:
    schedule_type = str(schedule_spec.get("type") or "daily").lower()
    if schedule_type == "daily":
        return list(range(start_offset, end_offset + 1))
    if schedule_type == "interval":
        every_days = int(schedule_spec.get("every_days") or schedule_spec.get("interval_days") or item.get("cycle_length_days") or 1)
        every_days = max(every_days, 1)
        return [day for day in range(start_offset, end_offset + 1) if (day - start_offset) % every_days == 0]
    if schedule_type == "cycle":
        cycle_length = int(schedule_spec.get("cycle_length_days") or schedule_spec.get("cycle_length") or item.get("cycle_length_days") or 28)
        days_on = _normalize_days_on(schedule_spec.get("days_on", item.get("days_on", [])), cycle_length)
        return [day for day in range(start_offset, end_offset + 1) if (((day - start_offset) % cycle_length) + 1) in days_on]
    if schedule_type == "hold":
        return []
    raise ExposureScheduleError(f"Unsupported exposure schedule type: {schedule_type}")


def _normalize_days_on(raw_days_on: Any, cycle_length: int) -> set[int]:
    if isinstance(raw_days_on, int):
        return set(range(1, min(int(raw_days_on), cycle_length) + 1))
    if isinstance(raw_days_on, list):
        return {int(day) for day in raw_days_on}
    if raw_days_on:
        return {int(raw_days_on)}
    return set(range(1, cycle_length + 1))


def _apply_interruption_days(
    schedule_payload: dict[str, Any],
    horizon_start: date,
    horizon_days: int,
    daily_by_drug: dict[str, list[float]],
    interruptions_by_drug: dict[str, set[int]],
) -> None:
    for interruption in schedule_payload.get("interruptions", []) or []:
        start = _parse_date(interruption.get("start_date"))
        end_value = interruption.get("end_date") or interruption.get("start_date")
        end = _parse_date(end_value)
        start_offset = max((start - horizon_start).days, 0)
        end_offset = min((end - horizon_start).days, horizon_days - 1)
        if end_offset < 0 or start_offset >= horizon_days or end_offset < start_offset:
            continue
        drug = str(interruption.get("drug") or "").strip().lower()
        target_drugs = [drug] if drug and drug in daily_by_drug else list(daily_by_drug)
        for target_drug in target_drugs:
            interruptions = interruptions_by_drug.setdefault(target_drug, set())
            daily = daily_by_drug[target_drug]
            for day in range(start_offset, end_offset + 1):
                daily[day] = 0.0
                interruptions.add(day)


def _finalize_profile(
    *,
    horizon_days: int,
    drug: str,
    daily: list[float],
    interruption_days: list[int],
    schedule_type: str,
    schedule_label: str,
    exposure_resolution_status: str,
) -> ExposureProfile:
    horizon = int(horizon_days)
    daily_values = [float(value) for value in daily[:horizon]] + [0.0] * max(horizon - len(daily), 0)
    cumulative = []
    total = 0.0
    dose_events = []
    admin_days = []
    for day, dose in enumerate(daily_values):
        total += float(dose)
        cumulative.append(float(total))
        if float(dose) > 0.0:
            admin_days.append(day)
            dose_events.append({"day": day, "dose_mg": float(dose), "drug": drug, "schedule_label": schedule_label})
    frequency = None
    if len(admin_days) >= 2:
        intervals = [later - earlier for earlier, later in zip(admin_days, admin_days[1:])]
        frequency = float(sum(intervals) / len(intervals))
    payload = {
        "horizon_days": horizon,
        "time_grid_days": list(range(horizon)),
        "drug": drug,
        "dose_events": dose_events,
        "daily_administered_dose_mg": daily_values,
        "cumulative_dose_mg": cumulative,
        "total_cumulative_dose_mg": float(total),
        "average_daily_dose_mg": float(total / max(horizon, 1)),
        "peak_administered_dose_mg": float(max(daily_values) if daily_values else 0.0),
        "dosing_frequency_days": frequency,
        "interruption_days": sorted(set(int(day) for day in interruption_days if 0 <= int(day) < horizon)),
        "interruption_fraction": float(len(set(interruption_days)) / max(horizon, 1)),
        "schedule_type": schedule_type,
        "schedule_label": schedule_label,
        "exposure_resolution_status": exposure_resolution_status,
    }
    digest = hashlib.sha256(json.dumps(payload, sort_keys=True, separators=(",", ":"), default=str).encode("utf-8")).hexdigest()
    return ExposureProfile(exposure_profile_hash=digest, **payload)


def _extract_numeric_dose(dose_entry: Any) -> float | None:
    if dose_entry is None or dose_entry == "":
        return None
    if isinstance(dose_entry, dict):
        dose_entry = dose_entry.get("dose")
    try:
        return float(dose_entry)
    except (TypeError, ValueError):
        return None


def _parse_date(value: Any) -> date:
    if isinstance(value, date):
        return value
    if value is None or value == "":
        raise ExposureScheduleError("schedule date is required")
    return date.fromisoformat(str(value))


def _profile_payload(profile: ExposureProfile | dict[str, Any]) -> dict[str, Any]:
    return profile.to_dict() if isinstance(profile, ExposureProfile) else dict(profile)