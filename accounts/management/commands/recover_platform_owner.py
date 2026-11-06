from django.core.exceptions import ValidationError
from django.core.management.base import BaseCommand, CommandError
from django.db import transaction
from django.utils import timezone

from accounts.models import AccessProfile, Capability
from accounts.services import audit, issue_capability, locked_policy


class Command(BaseCommand):
    help = "Explicit break-glass recovery. Revokes sessions/tokens before issuing a new owner link."

    def add_arguments(self, parser):
        parser.add_argument("--revoke-and-reissue", action="store_true", required=True)
        parser.add_argument(
            "--reset-authenticator",
            action="store_true",
            help="Explicit identity-verified recovery after authenticator loss.",
        )

    @transaction.atomic
    def handle(self, **options):
        locked_policy()
        try:
            identity = AccessProfile.objects.select_related("user").get(is_owner=True)
            identity.session_version += 1
            identity.access_state = AccessProfile.AccessState.PENDING
            if options["reset_authenticator"]:
                identity.totp_secret = ""
                identity.totp_confirmed = False
                identity.totp_last_step = -1
                audit("owner_authenticator_recovery", identity)
            identity.save()
            identity.user.set_unusable_password()
            identity.user.is_active = False
            identity.user.save(update_fields=["password", "is_active"])
            Capability.objects.filter(identity=identity, consumed_at=None).update(revoked_at=timezone.now())
            issue_capability(identity, "owner")
            audit("owner_break_glass_recovery", identity)
        except (AccessProfile.DoesNotExist, ValidationError) as exc:
            raise CommandError("Owner recovery unavailable; check governed identity configuration.") from exc
        self.stdout.write("OWNER_RECOVERY=QUEUED; OWNER_FIRST_LOGIN=PENDING_USER_ACTION")
