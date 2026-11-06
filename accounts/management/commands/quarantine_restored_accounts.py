from django.core.management.base import BaseCommand
from django.db import transaction
from django.utils import timezone

from accounts.models import AccessPolicy, AccessProfile, Capability, Notification
from accounts.services import audit, locked_policy


class Command(BaseCommand):
    help = "Post-restore quarantine: close pilot gates and invalidate restored sessions and capabilities."

    def add_arguments(self, parser):
        parser.add_argument("--post-restore-quarantine", action="store_true", required=True)

    @transaction.atomic
    def handle(self, **options):
        locked_policy()
        AccessPolicy.objects.filter(pk=1).update(requests_open=False, activations_open=False)
        for identity in AccessProfile.objects.select_related("user"):
            identity.session_version += 1
            identity.access_state = "PENDING_ACTIVATION" if identity.is_owner else "SUSPENDED"
            identity.save(update_fields=["session_version", "access_state"])
            if identity.user_id:
                identity.user.is_active = False
                identity.user.set_unusable_password()
                identity.user.save(update_fields=["is_active", "password"])
            audit("post_restore_quarantine", identity)
        Capability.objects.filter(revoked_at=None).update(revoked_at=timezone.now())
        Notification.objects.exclude(encrypted_payload="").update(
            encrypted_payload="", state="failed", error_code="RESTORE_QUARANTINE", attempts=5
        )
        self.stdout.write("RESTORE_QUARANTINE=PASS; OWNER_FIRST_LOGIN=PENDING_USER_ACTION")
