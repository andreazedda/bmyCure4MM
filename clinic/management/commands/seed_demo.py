from __future__ import annotations

import datetime

from django.core.management.base import BaseCommand

from clinic import models


class Command(BaseCommand):
    help = "Seed demo patients and regimens for quick testing."

    def handle(self, *args, **options):
        patient, _ = models.Patient.objects.get_or_create(
            mrn="MM001",
            defaults={
                "first_name": "John",
                "last_name": "Doe",
                "birth_date": datetime.date(1970, 1, 1),
                "sex": "M",
                "diagnosis_date": datetime.date(2022, 1, 1),
            },
        )
        regimen, _ = models.Regimen.objects.get_or_create(
            name="VRd",
            defaults={
                "line": "First-line",
                "components": "Bortezomib, Lenalidomide, Dexamethasone",
                "intent": "Induction",
            },
        )
        dara_regimen, _ = models.Regimen.objects.get_or_create(
            name="Dara-VRd",
            defaults={
                "line": "First-line",
                "components": "Daratumumab, Bortezomib, Lenalidomide, Dexamethasone",
                "intent": "Induction",
            },
        )
        self.stdout.write(
            self.style.SUCCESS(
                f"Seeded patient {patient.mrn} and regimens {regimen.name}, {dara_regimen.name}"
            )
        )
