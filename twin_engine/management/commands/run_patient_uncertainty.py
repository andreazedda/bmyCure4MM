from __future__ import annotations

import json

from django.core.management.base import BaseCommand, CommandError

from clinic.models import Patient
from twin_engine.models import CounterfactualRun
from twin_engine.provenance import CURRENT_MODEL_VERSION, record_simulation_metadata
from twin_engine.uncertainty import UncertaintyConfig, run_counterfactual_uncertainty


def _latest_runs_by_label(patient) -> list[CounterfactualRun]:
    latest = {}
    queryset = patient.counterfactual_runs.filter(status=CounterfactualRun.STATUS_COMPLETED).order_by("-id")
    for run in queryset:
        label = str((run.intervention_definition or {}).get("label") or f"Run {run.id}")
        latest.setdefault(label, run)
    return list(latest.values())


class Command(BaseCommand):
    help = "Run scenario uncertainty diagnostics for the latest completed scenarios of a patient."

    def add_arguments(self, parser):
        parser.add_argument("--patient-id", type=int, required=True)
        parser.add_argument("--horizon-days", type=int, required=True)
        parser.add_argument("--samples", type=int, default=100)
        parser.add_argument("--seed", type=int, default=17)

    def handle(self, *args, **options):
        patient = Patient.objects.filter(pk=options["patient_id"]).first()
        if patient is None:
            raise CommandError("Patient not found")
        runs = _latest_runs_by_label(patient)
        if not runs:
            raise CommandError("No completed counterfactual runs found")

        config = UncertaintyConfig(n_samples=options["samples"], random_seed=options["seed"])
        results = []
        for run in runs:
            uncertainty = run_counterfactual_uncertainty(patient, run.base_twin_state, run.intervention_definition or {}, options["horizon_days"], config)
            metrics = dict(run.comparison_metrics or {})
            metrics["uncertainty"] = uncertainty
            run.comparison_metrics = metrics
            run.save(update_fields=["comparison_metrics"])
            metadata = record_simulation_metadata(
                counterfactual_run=run,
                twin_state=run.base_twin_state,
                model_version=CURRENT_MODEL_VERSION,
                solver_name="counterfactual_uncertainty",
                input_payload={"counterfactual_run_id": run.id, "horizon_days": options["horizon_days"], "samples": options["samples"], "seed": options["seed"]},
                solver_parameters={"diagnostic_summary": {"status": uncertainty.get("status"), "metric_summaries": uncertainty.get("metric_summaries") or {}}},
                output_payload=uncertainty,
                random_seed=options["seed"],
            )
            results.append({"scenario_label": (run.intervention_definition or {}).get("label") or f"Run {run.id}", "counterfactual_run_id": run.id, "metadata_id": metadata.id, "uncertainty": uncertainty})
        self.stdout.write(json.dumps({"patient_id": patient.id, "results": results}, indent=2))
