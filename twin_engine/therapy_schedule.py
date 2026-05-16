from __future__ import annotations

from datetime import timedelta
from typing import Any

from django.db.models import Q


SUPPORTED_DRUG_ALIASES = {
    "len": "lenalidomide",
    "lenalidomide": "lenalidomide",
    "revlimid": "lenalidomide",
    "bor": "bortezomib",
    "bortezomib": "bortezomib",
    "velcade": "bortezomib",
    "dara": "daratumumab",
    "daratumumab": "daratumumab",
    "darzalex": "daratumumab",
    "carfilzomib": "carfilzomib",
    "carfil": "carfilzomib",
    "kyprolis": "carfilzomib",
    "dex": "dexamethasone",
    "dexamethasone": "dexamethasone",
}

SIMULATOR_DRUGS = {"lenalidomide", "bortezomib", "daratumumab", "carfilzomib"}


def build_therapy_schedule(patient, start_date, end_date) -> dict[str, Any]:
    therapies = list(
        patient.therapies.select_related("regimen").filter(
            start_date__lte=end_date,
        ).filter(
            Q(end_date__isnull=True) | Q(end_date__gte=start_date)
        ).order_by("start_date")
    )

    items: list[dict[str, Any]] = []
    missing_doses: list[dict[str, Any]] = []

    for therapy in therapies:
        item_end = therapy.end_date or end_date
        item = {
            "therapy_id": therapy.id,
            "regimen_id": therapy.regimen_id,
            "regimen_name": therapy.regimen.name,
            "start_date": max(therapy.start_date, start_date).isoformat(),
            "end_date": min(item_end, end_date).isoformat(),
            "cycle_length_days": therapy.cycle_length_days,
            "days_on": list(therapy.days_on or []),
            "components": parse_regimen_components(therapy.regimen),
            "doses": therapy.doses or {},
            "source_quality": therapy.source_quality,
            "schedule_notes": therapy.schedule_notes,
            "provenance": therapy.provenance or {},
        }
        items.append(item)

        for component in item["components"]:
            if not component["supported_by_solver"]:
                continue
            if component["drug"] not in (therapy.doses or {}):
                missing_doses.append(
                    {
                        "therapy_id": therapy.id,
                        "regimen_name": therapy.regimen.name,
                        "drug": component["drug"],
                        "reason": "missing structured dose",
                    }
                )

    schedule = {
        "patient_id": patient.id,
        "start_date": start_date.isoformat(),
        "end_date": end_date.isoformat(),
        "items": items,
        "missing_doses": missing_doses,
    }
    schedule["validation"] = validate_schedule(schedule)
    return schedule


def convert_patient_therapies_to_drug_doses(schedule, strict: bool = True) -> dict[str, Any]:
    if isinstance(schedule, dict):
        items = schedule.get("items", [])
    else:
        items = schedule

    aggregated = {drug: 0.0 for drug in SIMULATOR_DRUGS}
    missing_doses: list[dict[str, Any]] = []

    for item in items:
        cycle_length = item.get("cycle_length_days")
        days_on = [int(day) for day in (item.get("days_on") or [])]
        start_date = _parse_date(item.get("start_date"))
        end_date = _parse_date(item.get("end_date"))
        doses = item.get("doses") or {}

        for component in item.get("components", []):
            drug = component.get("drug")
            if drug not in SIMULATOR_DRUGS:
                continue

            dose_entry = doses.get(drug)
            dose_value = _extract_numeric_dose(dose_entry)
            if dose_value is None:
                missing_doses.append(
                    {
                        "therapy_id": item.get("therapy_id"),
                        "regimen_name": item.get("regimen_name"),
                        "drug": drug,
                        "reason": "dose missing or non-numeric",
                    }
                )
                continue

            if not cycle_length or not days_on:
                missing_doses.append(
                    {
                        "therapy_id": item.get("therapy_id"),
                        "regimen_name": item.get("regimen_name"),
                        "drug": drug,
                        "reason": "cycle_length_days or days_on missing",
                    }
                )
                continue

            administrations = _count_administrations(start_date, end_date, cycle_length, days_on)
            aggregated[drug] += float(dose_value) * float(administrations)

    if strict and missing_doses:
        details = ", ".join(
            f"therapy={item['therapy_id']} drug={item['drug']} ({item['reason']})"
            for item in missing_doses
        )
        raise ValueError("Structured therapy schedule incomplete: " + details)

    return {
        "drug_doses": aggregated,
        "missing_doses": missing_doses,
    }


def parse_regimen_components(regimen) -> list[dict[str, Any]]:
    raw_components = str(getattr(regimen, "components", "") or "")
    parsed: list[dict[str, Any]] = []
    for raw_component in raw_components.split(","):
        component = raw_component.strip()
        if not component:
            continue
        canonical = SUPPORTED_DRUG_ALIASES.get(component.lower(), component.lower())
        parsed.append(
            {
                "label": component,
                "drug": canonical,
                "supported_by_solver": canonical in SIMULATOR_DRUGS,
            }
        )
    return parsed


def validate_schedule(schedule) -> dict[str, Any]:
    missing_doses = list(schedule.get("missing_doses", []))
    for item in schedule.get("items", []):
        if item.get("cycle_length_days") and item.get("cycle_length_days") <= 0:
            missing_doses.append(
                {
                    "therapy_id": item.get("therapy_id"),
                    "reason": "cycle_length_days must be positive",
                }
            )
    return {
        "is_valid": len(missing_doses) == 0,
        "missing_doses": missing_doses,
    }


def _extract_numeric_dose(dose_entry) -> float | None:
    if dose_entry is None or dose_entry == "":
        return None
    if isinstance(dose_entry, dict):
        dose_entry = dose_entry.get("dose")
    try:
        return float(dose_entry)
    except (TypeError, ValueError):
        return None


def _parse_date(value):
    if hasattr(value, "year"):
        return value
    from datetime import date

    return date.fromisoformat(str(value))


def _count_administrations(start_date, end_date, cycle_length_days: int, days_on: list[int]) -> int:
    administrations = 0
    current_day = start_date
    normalized_days = {int(day) for day in days_on}
    while current_day <= end_date:
        cycle_day = ((current_day - start_date).days % int(cycle_length_days)) + 1
        if cycle_day in normalized_days:
            administrations += 1
        current_day += timedelta(days=1)
    return administrations
