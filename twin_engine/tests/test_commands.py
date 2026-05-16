from __future__ import annotations

from datetime import date
from io import StringIO

from django.core.management import call_command
from django.test import TestCase

from clinic.models import Patient
from twin_engine.management.commands.backfill_research_structured_data import (
    EVENT_ROWS,
    INTERRUPTION_ROWS,
    LAB_ROWS,
)
from twin_engine.models import AdverseEvent, LongitudinalLabResult, TherapyInterruption
from twin_engine.toxicity_model import compute_toxicity_constraints


class BackfillResearchStructuredDataCommandTests(TestCase):
    def setUp(self) -> None:
        self.patient = Patient.objects.create(
            mrn="MM-BACKFILL-001",
            first_name="Research",
            last_name="PatientBF",
            birth_date=date(1962, 8, 3),
            sex="F",
            diagnosis_date=date(2022, 5, 20),
        )

    def test_dry_run_creates_nothing(self) -> None:
        output = StringIO()
        call_command(
            "backfill_research_structured_data",
            patient_id=self.patient.id,
            dry_run=True,
            stdout=output,
        )
        self.assertEqual(LongitudinalLabResult.objects.filter(patient=self.patient).count(), 0)
        self.assertEqual(AdverseEvent.objects.filter(patient=self.patient).count(), 0)
        self.assertEqual(TherapyInterruption.objects.filter(patient=self.patient).count(), 0)

    def test_real_run_is_idempotent_and_matches_expected_counts(self) -> None:
        expected_labs = sum(len(row["items"]) for row in LAB_ROWS)
        expected_events = len(EVENT_ROWS)
        expected_interruptions = len(INTERRUPTION_ROWS)

        call_command("backfill_research_structured_data", patient_id=self.patient.id)
        call_command("backfill_research_structured_data", patient_id=self.patient.id)

        self.assertEqual(LongitudinalLabResult.objects.filter(patient=self.patient).count(), expected_labs)
        self.assertEqual(AdverseEvent.objects.filter(patient=self.patient).count(), expected_events)
        self.assertEqual(TherapyInterruption.objects.filter(patient=self.patient).count(), expected_interruptions)

    def test_toxicity_summary_uses_backfilled_values(self) -> None:
        call_command("backfill_research_structured_data", patient_id=self.patient.id)

        summary = compute_toxicity_constraints(self.patient)
        self.assertEqual(summary["liver"]["max_ast"], 406.0)
        self.assertEqual(summary["liver"]["max_alt"], 360.0)
        self.assertEqual(summary["liver"]["peak_date"], "2025-05-07")
        self.assertTrue(summary["neutropenia"]["absolute_neutropenia_event"])
        self.assertTrue(summary["infection"]["pneumonia_event"])
        self.assertTrue(summary["lenalidomide_toxicity_limited"])
        self.assertEqual(summary["toxicity_model_status"], "descriptive_only")