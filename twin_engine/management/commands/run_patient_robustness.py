from __future__ import annotations

import json

from django.core.management.base import BaseCommand, CommandError

from clinic.models import Patient
from twin_engine.management.commands.run_patient_uncertainty import _latest_runs_by_label
from twin_engine.provenance import record_simulation_metadata
from twin_engine.robustness import compute_robust_scenario_ranking
from twin_engine.state_model import get_current_twin_state


class Command(BaseCommand):
    help = "Compute robust scenario ranking from stored uncertainty summaries."

    def add_arguments(self, parser):
        parser.add_argument("--patient-id", type=int, required=True)

    def handle(self, *args, **options):
        patient = Patient.objects.filter(pk=options["patient_id"]).first()
        if patient is None:
            raise CommandError("Patient not found")
        runs = _latest_runs_by_label(patient)
        result = compute_robust_scenario_ranking(runs)
        current_state = get_current_twin_state(patient)
        metadata = record_simulation_metadata(
            twin_state=current_state,
            model_id="counterfactual_model",
            solver_name="robust_scenario_ranking",
            input_payload={"patient_id": patient.id, "counterfactual_run_ids": [run.id for run in runs]},
            solver_parameters={"diagnostic_summary": result},
            output_payload=result,
        )
        self.stdout.write(json.dumps({"patient_id": patient.id, "metadata_id": metadata.id, "robustness": result}, indent=2))
