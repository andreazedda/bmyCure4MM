from __future__ import annotations

import json

from django.core.management.base import (
    BaseCommand,
    CommandError,
)

from clinic.models import Patient
from twin_engine.research_dataset_import import (
    DatasetImportError,
    import_dataset,
    load_dataset_bundle,
)


class Command(BaseCommand):
    help = (
        "Validate and import a versioned structured research "
        "dataset for a pseudonymized patient."
    )

    def add_arguments(self, parser):
        parser.add_argument(
            "--patient-id",
            type=int,
            required=True,
        )
        parser.add_argument(
            "--dataset",
            required=True,
            help="Path to private dataset.json",
        )
        parser.add_argument(
            "--manifest",
            help=(
                "Path to manifest.json. Defaults to the "
                "dataset file's sibling manifest.json."
            ),
        )
        parser.add_argument(
            "--dry-run",
            action="store_true",
        )

    def handle(self, *args, **options):
        patient = Patient.objects.filter(
            pk=options["patient_id"]
        ).first()

        if patient is None:
            raise CommandError("Patient not found")

        try:
            bundle = load_dataset_bundle(
                options["dataset"],
                options.get("manifest"),
            )

            result = import_dataset(
                patient,
                bundle,
                dry_run=options["dry_run"],
            )

        except DatasetImportError as exc:
            raise CommandError(str(exc)) from exc

        self.stdout.write(
            json.dumps(
                result,
                indent=2,
                sort_keys=True,
                default=str,
            )
        )

        if result["conflicts"] and not options["dry_run"]:
            raise CommandError(
                "Dataset import blocked by semantic conflicts; "
                "no database mutations were applied."
            )
