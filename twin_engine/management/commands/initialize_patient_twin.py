from __future__ import annotations

import json

from django.core.management.base import BaseCommand, CommandError

from clinic.models import Assessment, Patient
from twin_engine.state_model import initialize_from_assessment, preview_state_from_assessment


class Command(BaseCommand):
    help = "Initialize a persisted research twin state from a patient assessment."

    def add_arguments(self, parser):
        parser.add_argument("--patient-id", type=int, required=True)
        parser.add_argument("--assessment-id", type=int, required=True)
        parser.add_argument("--dry-run", action="store_true")

    def handle(self, *args, **options):
        patient = Patient.objects.filter(pk=options["patient_id"]).first()
        if patient is None:
            raise CommandError("Patient not found")
        assessment = Assessment.objects.filter(pk=options["assessment_id"], patient=patient).first()
        if assessment is None:
            raise CommandError("Assessment not found for the requested patient")

        if options["dry_run"]:
            payload = preview_state_from_assessment(assessment)
            self.stdout.write(json.dumps({
                "dry_run": True,
                "patient_id": patient.id,
                "assessment_id": assessment.id,
                **payload,
            }, indent=2, default=str))
            return

        state = initialize_from_assessment(assessment)
        self.stdout.write(json.dumps({
            "dry_run": False,
            "patient_id": patient.id,
            "assessment_id": assessment.id,
            "state_id": state.id,
            "state_date": state.state_date.isoformat(),
            "method": state.method,
            "model_version": state.model_version,
        }, indent=2))
