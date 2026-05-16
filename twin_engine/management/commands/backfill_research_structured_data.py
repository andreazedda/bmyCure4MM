from __future__ import annotations

import json
from datetime import date

from django.core.management.base import BaseCommand, CommandError
from django.db import transaction

from clinic.models import Assessment, Patient, PatientTherapy
from twin_engine.models import AdverseEvent, LongitudinalLabResult, TherapyInterruption


CASE_LABEL = "MM_RESEARCH_CASE_001"
COMMON_PROVENANCE = {
    "source": "uploaded_clinical_documentation",
    "extraction_status": "manual_from_documents",
    "case_label": CASE_LABEL,
}

LAB_ROWS = [
    {"date": "2023-10-03", "items": [("HB", 10.5, "g/dL"), ("WBC", 3860, "/uL"), ("PLT", 212000, "/uL"), ("CREATININE", 0.69, "mg/dL")]},
    {"date": "2024-01-05", "items": [("HB", 11.2, "g/dL"), ("NEU", 3130, "/uL"), ("PLT", 180000, "/uL"), ("CREATININE", 0.68, "mg/dL"), ("FLC_RATIO", 1.78, "ratio")]},
    {"date": "2024-03-18", "items": [("HB", 10.9, "g/dL"), ("FLC_RATIO", 1.56, "ratio")]},
    {"date": "2024-05-20", "items": [("HB", 11.8, "g/dL"), ("WBC", 3520, "/uL"), ("NEU", 1280, "/uL"), ("CREATININE", 0.58, "mg/dL"), ("M_PROTEIN", 0.1, "g/dL"), ("KAPPA_FLC", 82, "mg/L"), ("LAMBDA_FLC", 33.93, "mg/L"), ("FLC_RATIO", 2.43, "ratio")]},
    {"date": "2024-07-08", "items": [("HB", 11.0, "g/dL"), ("WBC", 2720, "/uL"), ("NEU", 710, "/uL"), ("PLT", 159000, "/uL"), ("M_PROTEIN", 0.1, "g/dL"), ("FLC_RATIO", 1.94, "ratio")]},
    {"date": "2024-08-21", "items": [("HB", 11.1, "g/dL"), ("NEU", 1240, "/uL"), ("CREATININE", 0.60, "mg/dL"), ("M_PROTEIN", 0.1, "g/dL"), ("FLC_RATIO", 2.06, "ratio")]},
    {"date": "2024-10-02", "items": [("HB", 11.5, "g/dL"), ("NEU", 1260, "/uL"), ("CREATININE", 0.61, "mg/dL"), ("M_PROTEIN", 0.1, "g/dL"), ("KAPPA_FLC", 66, "mg/L"), ("LAMBDA_FLC", 28, "mg/L"), ("FLC_RATIO", 2.30, "ratio")]},
    {"date": "2024-11-13", "items": [("M_PROTEIN", 0.1, "g/dL"), ("FLC_RATIO", 2.46, "ratio")]},
    {"date": "2025-02-03", "items": [("HB", 11.4, "g/dL"), ("M_PROTEIN", 0.1, "g/dL"), ("FLC_RATIO", 2.19, "ratio")]},
    {"date": "2025-03-04", "items": [("HB", 11.3, "g/dL"), ("CREATININE", 0.61, "mg/dL"), ("M_PROTEIN", 0.0, "g/dL"), ("FLC_RATIO", 2.73, "ratio"), ("AST", 25, "U/L"), ("ALT", 43, "U/L")]},
    {"date": "2025-04-23", "items": [("HB", 11.3, "g/dL"), ("CREATININE", 0.57, "mg/dL"), ("FLC_RATIO", 3.08, "ratio")]},
    {"date": "2025-05-07", "items": [("AST", 406, "U/L"), ("ALT", 360, "U/L")]},
    {"date": "2025-05-08", "items": [("AST", 76, "U/L"), ("ALT", 209, "U/L")]},
    {"date": "2025-05-13", "items": [("AST", 19, "U/L"), ("ALT", 27, "U/L")]},
]

EVENT_ROWS = [
    {"date": "2025-01-10", "event_type": AdverseEvent.TYPE_PNEUMONIA, "suspected_drug": "", "observed_values": {"left_basal_pneumonia": True}, "action_taken": "Supportive anti-infective management documented.", "outcome": "Pneumonia treated.", "source_quality": "clinical_record"},
    {"date": "2025-01-10", "event_type": AdverseEvent.TYPE_NEUTROPENIA, "suspected_drug": "", "observed_values": {"absolute_neutropenia": True}, "action_taken": "Growth-factor support documented.", "outcome": "Temporary interruption and support.", "source_quality": "clinical_record"},
    {"date": "2025-01-10", "event_type": AdverseEvent.TYPE_DIARRHEA, "suspected_drug": "", "observed_values": {}, "action_taken": "Supportive care documented.", "outcome": "Intercurrent event noted.", "source_quality": "clinical_record"},
    {"date": "2025-05-07", "event_type": AdverseEvent.TYPE_HYPERTRANSAMINASEMIA, "suspected_drug": "lenalidomide", "observed_values": {"AST": 406, "ALT": 360}, "action_taken": "Lenalidomide interruption; supportive medication review.", "outcome": "Improved on repeat evaluation.", "source_quality": "clinical_record"},
    {"date": "2025-05-26", "event_type": AdverseEvent.TYPE_HYPERTRANSAMINASEMIA, "suspected_drug": "lenalidomide", "observed_values": {"alt_altered": True}, "action_taken": "Reduced to alternate-day lenalidomide.", "outcome": "Ongoing monitoring planned.", "source_quality": "clinical_record"},
    {"date": "2025-06-18", "event_type": AdverseEvent.TYPE_HEPATIC_STEATOSIS, "suspected_drug": "", "observed_values": {"high_grade_steatosis": True, "common_bile_duct_mm": 7, "spleen_mm": 65, "right_renal_cyst_mm": 18}, "action_taken": "Ultrasound documented; ongoing review.", "outcome": "Steatosis documented.", "source_quality": "clinical_record"},
    {"date": "2025-06-18", "event_type": AdverseEvent.TYPE_PAIN, "suspected_drug": "", "observed_values": {}, "action_taken": "Fentanyl patch 25 to 37 mcg/h; possible 50 mcg/h if persistent and tolerated.", "outcome": "Pain-management escalation documented.", "source_quality": "clinical_record"},
]

INTERRUPTION_ROWS = [
    {"start_date": "2024-10-14", "end_date": "2024-10-21", "drug": "lenalidomide", "reason": TherapyInterruption.REASON_HYPERTRANSAMINASEMIA, "evidence": {"documented_issue": "hypertransaminasemia"}, "action_taken": "Interrupted then resumed after normalization."},
    {"start_date": "2025-01-10", "end_date": "2025-01-28", "drug": "lenalidomide", "reason": TherapyInterruption.REASON_INFECTION, "evidence": {"absolute_neutropenia": True, "infection_event": "pneumonia"}, "action_taken": "Held during pneumonia and neutropenia support period."},
    {"start_date": "2025-05-07", "end_date": "2025-05-19", "drug": "lenalidomide", "reason": TherapyInterruption.REASON_HYPERTRANSAMINASEMIA, "evidence": {"AST": 406, "ALT": 360}, "action_taken": "Interrupted until liver values normalized."},
]


class Command(BaseCommand):
    help = "Backfill structured longitudinal labs, adverse events, and therapy interruptions for a pseudonymized research case."

    def add_arguments(self, parser):
        parser.add_argument("--patient-id", type=int, required=True)
        parser.add_argument("--dry-run", action="store_true")

    def handle(self, *args, **options):
        patient = Patient.objects.filter(pk=options["patient_id"]).first()
        if patient is None:
            raise CommandError("Patient not found")

        payload = _execute_backfill(patient, dry_run=options["dry_run"])
        self.stdout.write(json.dumps(payload, indent=2, default=str))


def _execute_backfill(patient, *, dry_run: bool) -> dict[str, object]:
    counters = {
        "patient_id": patient.id,
        "dry_run": dry_run,
        "labs_created": 0,
        "labs_updated": 0,
        "events_created": 0,
        "events_updated": 0,
        "interruptions_created": 0,
        "interruptions_updated": 0,
    }

    def _record_lab(result_date: str, analyte: str, value: float | None, unit: str):
        assessment = Assessment.objects.filter(patient=patient, date=date.fromisoformat(result_date)).first()
        existing = LongitudinalLabResult.objects.filter(patient=patient, date=result_date, analyte=analyte).first()
        defaults = {
            "assessment": assessment,
            "value": value,
            "unit": unit,
            "source_quality": LongitudinalLabResult.SOURCE_QUALITY_EXTRACTED_DOCUMENT,
            "provenance": COMMON_PROVENANCE,
            "notes": "Backfilled structured longitudinal lab result for pseudonymized research case.",
        }
        if dry_run:
            counters["labs_created" if existing is None else "labs_updated"] += 1
            return
        _, created = LongitudinalLabResult.objects.update_or_create(
            patient=patient,
            date=result_date,
            analyte=analyte,
            defaults=defaults,
        )
        counters["labs_created" if created else "labs_updated"] += 1

    def _record_event(item: dict[str, object]):
        existing = AdverseEvent.objects.filter(patient=patient, date=item["date"], event_type=item["event_type"]).first()
        defaults = {
            "grade": item.get("grade", ""),
            "suspected_drug": item.get("suspected_drug", ""),
            "observed_values": item.get("observed_values", {}),
            "action_taken": item.get("action_taken", ""),
            "outcome": item.get("outcome", ""),
            "provenance": COMMON_PROVENANCE,
        }
        if dry_run:
            counters["events_created" if existing is None else "events_updated"] += 1
            return
        _, created = AdverseEvent.objects.update_or_create(
            patient=patient,
            date=item["date"],
            event_type=item["event_type"],
            defaults=defaults,
        )
        counters["events_created" if created else "events_updated"] += 1

    def _record_interruption(item: dict[str, object]):
        therapy = _resolve_patient_therapy(patient, item["drug"], item["start_date"], item.get("end_date"))
        existing = TherapyInterruption.objects.filter(
            patient=patient,
            drug=item["drug"],
            start_date=item["start_date"],
            end_date=item.get("end_date"),
        ).first()
        defaults = {
            "patient_therapy": therapy,
            "reason": item["reason"],
            "evidence": {**COMMON_PROVENANCE, **(item.get("evidence", {}) or {})},
            "action_taken": item.get("action_taken", ""),
            "source_quality": "clinical_record",
        }
        if dry_run:
            counters["interruptions_created" if existing is None else "interruptions_updated"] += 1
            return
        _, created = TherapyInterruption.objects.update_or_create(
            patient=patient,
            drug=item["drug"],
            start_date=item["start_date"],
            end_date=item.get("end_date"),
            defaults=defaults,
        )
        counters["interruptions_created" if created else "interruptions_updated"] += 1

    def _run():
        for row in LAB_ROWS:
            for analyte, value, unit in row["items"]:
                _record_lab(row["date"], analyte, value, unit)
        for item in EVENT_ROWS:
            _record_event(item)
        for item in INTERRUPTION_ROWS:
            _record_interruption(item)

    if dry_run:
        _run()
        return counters

    with transaction.atomic():
        _run()
    return counters


def _resolve_patient_therapy(patient, drug: str, start_date: str, end_date: str | None):
    start = date.fromisoformat(start_date)
    end = date.fromisoformat(end_date) if end_date else None
    therapies = patient.therapies.select_related("regimen").order_by("start_date")
    for therapy in therapies:
        therapy_end = therapy.end_date or end or start
        if therapy.start_date > (end or start):
            continue
        if therapy_end < start:
            continue
        components = (therapy.regimen.components or "").lower()
        dose_keys = {str(key).lower() for key in (therapy.doses or {}).keys()}
        if drug.lower() in components or drug.lower() in dose_keys:
            return therapy
    return None