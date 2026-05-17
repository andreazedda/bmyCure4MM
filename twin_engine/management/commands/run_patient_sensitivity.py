from __future__ import annotations

import json

from django.core.management.base import BaseCommand, CommandError

from clinic.models import Patient
from twin_engine.management.commands.run_patient_uncertainty import _latest_runs_by_label
from twin_engine.provenance import CURRENT_MODEL_VERSION, record_simulation_metadata
from twin_engine.sensitivity import run_counterfactual_sensitivity


class Command(BaseCommand):
    help = "Run one-at-a-time sensitivity diagnostics for the latest completed scenarios of a patient."

    def add_arguments(self, parser):
        parser.add_argument("--patient-id", type=int, required=True)
        parser.add_argument("--horizon-days", type=int, required=True)

    def handle(self, *args, **options):
        patient = Patient.objects.filter(pk=options["patient_id"]).first()
        if patient is None:
            raise CommandError("Patient not found")
        runs = _latest_runs_by_label(patient)
        if not runs:
            raise CommandError("No completed counterfactual runs found")

        results = []
        for run in runs:
            sensitivity = run_counterfactual_sensitivity(patient, run.base_twin_state, run.intervention_definition or {}, options["horizon_days"])
            metrics = dict(run.comparison_metrics or {})
            metrics["sensitivity"] = sensitivity
            run.comparison_metrics = metrics
            run.save(update_fields=["comparison_metrics"])
            metadata = record_simulation_metadata(
                counterfactual_run=run,
                twin_state=run.base_twin_state,
                model_version=CURRENT_MODEL_VERSION,
                solver_name="counterfactual_sensitivity",
                input_payload={"counterfactual_run_id": run.id, "horizon_days": options["horizon_days"]},
                solver_parameters={"diagnostic_summary": {"status": sensitivity.get("status"), "top_drivers": sensitivity.get("top_drivers") or []}},
                output_payload=sensitivity,
            )
            results.append({"scenario_label": (run.intervention_definition or {}).get("label") or f"Run {run.id}", "counterfactual_run_id": run.id, "metadata_id": metadata.id, "sensitivity": sensitivity})
        self.stdout.write(json.dumps({"patient_id": patient.id, "results": results}, indent=2))
