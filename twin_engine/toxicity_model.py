from __future__ import annotations

from typing import Any

from .models import AdverseEvent, LongitudinalLabResult, TherapyInterruption


def extract_toxicity_history(patient) -> dict[str, Any]:
    return {
        "labs": list(
            LongitudinalLabResult.objects.filter(
                patient=patient,
                analyte__in=[
                    LongitudinalLabResult.ANALYTE_AST,
                    LongitudinalLabResult.ANALYTE_ALT,
                    LongitudinalLabResult.ANALYTE_NEU,
                    LongitudinalLabResult.ANALYTE_WBC,
                    LongitudinalLabResult.ANALYTE_PLT,
                ],
            ).order_by("date", "analyte")
        ),
        "events": list(patient.adverse_events.order_by("date", "event_type")),
        "interruptions": list(patient.therapy_interruptions.order_by("start_date")),
        "therapies": list(patient.therapies.select_related("regimen").order_by("start_date")),
    }


def summarize_liver_toxicity(patient) -> dict[str, Any]:
    history = extract_toxicity_history(patient)
    ast_entries = [item for item in history["labs"] if item.analyte == LongitudinalLabResult.ANALYTE_AST and item.value is not None]
    alt_entries = [item for item in history["labs"] if item.analyte == LongitudinalLabResult.ANALYTE_ALT and item.value is not None]

    max_ast_entry = max(ast_entries, key=lambda item: float(item.value), default=None)
    max_alt_entry = max(alt_entries, key=lambda item: float(item.value), default=None)
    peak_candidates = [item for item in (max_ast_entry, max_alt_entry) if item is not None]
    peak_entry = max(peak_candidates, key=lambda item: float(item.value), default=None)

    hyper_events = [item for item in history["events"] if item.event_type == AdverseEvent.TYPE_HYPERTRANSAMINASEMIA]
    hepatic_steatosis = any(item.event_type == AdverseEvent.TYPE_HEPATIC_STEATOSIS for item in history["events"])
    recurrence_after_restart: bool | str
    if len(hyper_events) >= 2:
        recurrence_after_restart = True
    elif hyper_events:
        recurrence_after_restart = "unknown"
    else:
        recurrence_after_restart = "unknown"

    return {
        "max_ast": float(max_ast_entry.value) if max_ast_entry is not None else None,
        "max_alt": float(max_alt_entry.value) if max_alt_entry is not None else None,
        "peak_date": peak_entry.date.isoformat() if peak_entry is not None else None,
        "recurrence_after_restart": recurrence_after_restart,
        "hepatic_steatosis": hepatic_steatosis if hepatic_steatosis else "unknown",
        "toxicity_model_status": "descriptive_only",
    }


def summarize_neutropenia_history(patient) -> dict[str, Any]:
    history = extract_toxicity_history(patient)
    neu_entries = [item for item in history["labs"] if item.analyte == LongitudinalLabResult.ANALYTE_NEU and item.value is not None]
    min_neu_entry = min(neu_entries, key=lambda item: float(item.value), default=None)
    absolute_neutropenia_event = any(item.event_type == AdverseEvent.TYPE_NEUTROPENIA for item in history["events"])
    nivestim_support = any(
        "nivestim" in (therapy.regimen.name or "").lower()
        or "filgrastim" in (therapy.regimen.components or "").lower()
        for therapy in history["therapies"]
    )
    return {
        "min_neu": float(min_neu_entry.value) if min_neu_entry is not None else None,
        "absolute_neutropenia_event": absolute_neutropenia_event,
        "nivestim_support": nivestim_support,
        "toxicity_model_status": "descriptive_only",
    }


def summarize_infection_history(patient) -> dict[str, Any]:
    history = extract_toxicity_history(patient)
    pneumonia_event = any(item.event_type == AdverseEvent.TYPE_PNEUMONIA for item in history["events"])
    upper_respiratory_infection_event = any(
        item.event_type == AdverseEvent.TYPE_UPPER_RESPIRATORY_INFECTION for item in history["events"]
    )
    if upper_respiratory_infection_event:
        upper_respiratory = True
    elif pneumonia_event or any(item.event_type in {AdverseEvent.TYPE_PNEUMONIA, AdverseEvent.TYPE_OTHER} for item in history["events"]):
        upper_respiratory = False
    else:
        upper_respiratory = "unknown"
    return {
        "pneumonia_event": pneumonia_event,
        "upper_respiratory_infection_event": upper_respiratory,
        "toxicity_model_status": "descriptive_only",
    }


def compute_toxicity_constraints(patient) -> dict[str, Any]:
    liver = summarize_liver_toxicity(patient)
    neutropenia = summarize_neutropenia_history(patient)
    infection = summarize_infection_history(patient)
    interruptions = list(patient.therapy_interruptions.order_by("start_date"))

    lenalidomide_toxicity_limited = any(
        (interruption.drug or "").lower() == "lenalidomide"
        and interruption.reason in {
            TherapyInterruption.REASON_HYPERTRANSAMINASEMIA,
            TherapyInterruption.REASON_NEUTROPENIA,
            TherapyInterruption.REASON_INFECTION,
        }
        for interruption in interruptions
    )

    return {
        "liver": liver,
        "neutropenia": neutropenia,
        "infection": infection,
        "lenalidomide_toxicity_limited": lenalidomide_toxicity_limited,
        "toxicity_targets_modeled": False,
        "reason": "descriptive toxicity constraints only; no predictive AST/ALT/NEU model yet",
        "toxicity_model_status": "descriptive_only",
    }