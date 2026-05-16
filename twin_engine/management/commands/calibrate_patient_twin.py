from __future__ import annotations

import json
from datetime import date

from django.core.management.base import BaseCommand, CommandError

from clinic.models import Patient
from twin_engine.calibration import calibrate_patient_parameters
from twin_engine.state_model import get_current_twin_state, initialize_from_assessment


class Command(BaseCommand):
    help = "Calibrate a research twin state from longitudinal patient history."

    def add_arguments(self, parser):
        parser.add_argument("--patient-id", type=int, required=True)
        parser.add_argument("--start-date")
        parser.add_argument("--end-date")
        parser.add_argument("--dry-run", action="store_true")

    def handle(self, *args, **options):
        patient = Patient.objects.prefetch_related("assessments", "therapies").filter(pk=options["patient_id"]).first()
        if patient is None:
            raise CommandError("Patient not found")

        assessments = list(patient.assessments.order_by("date"))
        if options.get("start_date"):
            start_date = date.fromisoformat(options["start_date"])
            assessments = [item for item in assessments if item.date >= start_date]
        if options.get("end_date"):
            end_date = date.fromisoformat(options["end_date"])
            assessments = [item for item in assessments if item.date <= end_date]

        if not assessments:
            raise CommandError("No assessments available in the selected interval")

        current_state = get_current_twin_state(patient) or initialize_from_assessment(assessments[0])

        if options["dry_run"]:
            self.stdout.write(json.dumps({
                "dry_run": True,
                "patient_id": patient.id,
                "current_state_id": current_state.id,
                "assessment_count": len(assessments),
                "therapy_count": patient.therapies.count(),
            }, indent=2))
            return

        result = calibrate_patient_parameters(patient, current_state, assessments, patient.therapies.all())
        self.stdout.write(json.dumps({
            "dry_run": False,
            "patient_id": patient.id,
            "state_id": result["state"].id,
            "optimizer": result["optimizer"],
            "residual_count": len(result["residuals"]),
            "calibration_status": result["state"].parameter_uncertainty.get("calibration_status"),
            "objective_before": result["state"].parameter_uncertainty.get("objective_before"),
            "objective_after": result["state"].parameter_uncertainty.get("objective_after"),
            "rmse_before": result["state"].parameter_uncertainty.get("rmse_before"),
            "rmse_after": result["state"].parameter_uncertainty.get("rmse_after"),
            "mae_before": result["state"].parameter_uncertainty.get("mae_before"),
            "mae_after": result["state"].parameter_uncertainty.get("mae_after"),
        }, indent=2, default=str))
