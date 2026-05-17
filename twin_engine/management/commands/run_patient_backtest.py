from __future__ import annotations

import json

from django.core.management.base import BaseCommand, CommandError

from clinic.models import Patient
from twin_engine.backtesting import run_patient_backtest
from twin_engine.provenance import CURRENT_MODEL_VERSION, record_simulation_metadata
from twin_engine.state_model import get_current_twin_state


class Command(BaseCommand):
    help = "Run rolling-origin backtesting for a patient without modifying source clinical records."

    def add_arguments(self, parser):
        parser.add_argument("--patient-id", type=int, required=True)
        parser.add_argument("--minimum-history-points", type=int, default=3)

    def handle(self, *args, **options):
        patient = Patient.objects.filter(pk=options["patient_id"]).first()
        if patient is None:
            raise CommandError("Patient not found")
        current_state = get_current_twin_state(patient)
        result = run_patient_backtest(patient, current_state, minimum_history_points=options["minimum_history_points"])
        record = record_simulation_metadata(
            twin_state=current_state,
            model_version=CURRENT_MODEL_VERSION,
            solver_name="rolling_origin_backtest",
            input_payload={"patient_id": patient.id, "minimum_history_points": options["minimum_history_points"]},
            solver_parameters={"diagnostic_summary": result},
            output_payload=result,
        )
        self.stdout.write(json.dumps({"patient_id": patient.id, "metadata_id": record.id, "backtest": result}, indent=2))
