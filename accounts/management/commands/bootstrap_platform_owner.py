from django.core.exceptions import ValidationError
from django.core.management.base import BaseCommand, CommandError

from accounts.services import bootstrap_owner


class Command(BaseCommand):
    help = "Provision the configured owner once after migrations. Never prints credentials."

    def handle(self, **options):
        try:
            _, created = bootstrap_owner()
        except (ValidationError, ValueError) as exc:
            raise CommandError("Owner bootstrap rejected: identity/configuration conflict.") from exc
        self.stdout.write("OWNER_BOOTSTRAP=CREATED" if created else "OWNER_BOOTSTRAP=PRESERVED")
