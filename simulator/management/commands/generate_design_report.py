from __future__ import annotations

from pathlib import Path

from django.core.management.base import BaseCommand

from simulator.design.reporting import run_design_report, write_report_json


class Command(BaseCommand):
    help = "Generate a demo MM design/simulation report JSON (tradeoffs + dynamics + toxicity rationale)."

    def add_arguments(self, parser):
        parser.add_argument("--seed", type=int, default=7)
        parser.add_argument("--steps", type=int, default=18)
        parser.add_argument("--dt-days", type=float, default=7.0)
        parser.add_argument(
            "--out",
            type=str,
            default="artifacts/design_report.json",
            help="Output path (relative to project root by default)",
        )

    def handle(self, *args, **options):
        seed = int(options["seed"])
        steps = int(options["steps"])
        dt_days = float(options["dt_days"])
        out = Path(options["out"]).as_posix()

        report = run_design_report(seed=seed, steps=steps, dt_days=dt_days)
        path = write_report_json(report, out)

        self.stdout.write(self.style.SUCCESS(f"Wrote design report: {path}"))
        self.stdout.write(f"Chosen therapy: {report.rationale.get('selected_therapy')}")
