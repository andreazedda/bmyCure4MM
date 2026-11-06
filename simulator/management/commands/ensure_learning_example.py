"""Fill an empty learning catalogue using the existing synthetic Scenario model."""

from django.core.management.base import BaseCommand
from django.db import transaction

from accounts.services import locked_policy
from simulator.models import Scenario


class Command(BaseCommand):
    help = "Create one independent synthetic learning example only when the scenario catalogue is empty."

    def handle(self, *args, **options):
        # Reuse the existing cross-process DB write lock; preserve intentional inactive catalogues.
        with transaction.atomic():
            locked_policy()
            if Scenario.objects.exists():
                self.stdout.write("LEARNING_CATALOGUE=PRESERVED")
                return
            Scenario.objects.create(
                title="Synthetic model exploration / Esplorazione sintetica",
                summary=(
                    "Independent synthetic example using the existing model defaults. "
                    "Vary one numeric input at a time and compare simulated outputs. "
                    "Esempio sintetico indipendente: modifica un input alla volta e confronta gli output simulati. "
                    "No real patient data; research only / Nessun dato reale di pazienti; solo ricerca."
                ),
                clinical_stage="newly_diagnosed",
                risk_stratification="Synthetic / Sintetico",
                active=True,
            )
        self.stdout.write("LEARNING_CATALOGUE=INITIALIZED_SYNTHETIC")
