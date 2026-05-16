from __future__ import annotations

import json
from pathlib import Path

from django.core.management.base import BaseCommand, CommandError

from clinic.models import Patient
from twin_engine.counterfactual import run_counterfactual
from twin_engine.state_model import get_current_twin_state


class Command(BaseCommand):
    help = "Run a patient-specific what-if counterfactual research simulation."

    def add_arguments(self, parser):
        parser.add_argument("--patient-id", type=int, required=True)
        parser.add_argument("--base-state-id", type=int)
        parser.add_argument("--intervention-json", required=True)
        parser.add_argument("--horizon-days", type=int, required=True)
        parser.add_argument("--dry-run", action="store_true")

    def handle(self, *args, **options):
        patient = Patient.objects.filter(pk=options["patient_id"]).first()
        if patient is None:
            raise CommandError("Patient not found")

        intervention_path = Path(options["intervention_json"])
        if not intervention_path.exists():
            raise CommandError("Intervention JSON file not found")
        intervention_definition = json.loads(intervention_path.read_text(encoding="utf-8"))

        base_state = patient.twin_states.filter(pk=options.get("base_state_id")).first() if options.get("base_state_id") else get_current_twin_state(patient)
        if base_state is None:
            raise CommandError("Base twin state not found")

        if options["dry_run"]:
            self.stdout.write(json.dumps({
                "dry_run": True,
                "patient_id": patient.id,
                "base_state_id": base_state.id,
                "horizon_days": options["horizon_days"],
                "intervention_definition": intervention_definition,
            }, indent=2))
            return

        run = run_counterfactual(
            patient,
            base_state,
            intervention_definition,
            options["horizon_days"],
        )
        self.stdout.write(json.dumps({
            "dry_run": False,
            "patient_id": patient.id,
            "counterfactual_run_id": run.id,
            "status": run.status,
            "predicted_biomarkers": (run.simulation_summary or {}).get("predicted_biomarkers"),
            "simulation_summary": run.simulation_summary,
            "comparison_metrics": run.comparison_metrics,
            "provenance_id": run.metadata_records.order_by("-id").values_list("id", flat=True).first(),
            "report_artifact": run.report_artifact,
            "trajectory_artifact": run.trajectory_artifact,
        }, indent=2))
