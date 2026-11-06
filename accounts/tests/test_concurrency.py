from concurrent.futures import ThreadPoolExecutor

from django.contrib.auth import get_user_model
from django.core.exceptions import ValidationError
from django.db import connection
from django.test import TransactionTestCase, override_settings

from accounts import services
from accounts.models import AccessProfile, Notification
from accounts.notifications import dispatch_notification

from .helpers import TEST_SETTINGS, LifecycleMixin, token_for


@override_settings(**TEST_SETTINGS)
class RegistrationConcurrencyTests(LifecycleMixin, TransactionTestCase):
    def setUp(self):
        self.initialize()

    def concurrent(self, operation):
        def run(_):
            try:
                return operation()
            except ValidationError:
                return "already consumed"
            finally:
                connection.close()

        with ThreadPoolExecutor(max_workers=2) as executor:
            return list(executor.map(run, range(2)))

    def test_duplicate_request_approval_activation_and_delivery(self):
        self.concurrent(
            lambda: services.request_access(
                {"email": "race@example.test", "nickname": "Synthetic race", "reason": "Concurrency test"}
            )
        )
        identity = AccessProfile.objects.get(email="race@example.test")
        self.assertEqual(Notification.objects.filter(identity=identity, kind="verify").count(), 1)
        services.confirm_capability(token_for(identity, "verify"), "verify")
        self.concurrent(lambda: services.decide_request(self.owner.user, identity.pk, approve=True))
        identity.refresh_from_db()
        self.assertEqual(get_user_model().objects.filter(email=identity.email).count(), 1)
        self.assertEqual(Notification.objects.filter(identity=identity, kind="activate").count(), 1)
        raw = token_for(identity, "activate")
        self.concurrent(lambda: services.confirm_capability(raw, "activate", self.password))
        identity.refresh_from_db()
        self.assertEqual(identity.access_state, "ACTIVE")
        self.assertEqual(identity.session_version, 2)
        Notification.objects.exclude(identity=identity, kind="activate").update(
            state="sent", encrypted_payload=""
        )
        # The consumed activation is discarded, never sent by a racing dispatcher.
        self.concurrent(dispatch_notification)
        item = Notification.objects.get(identity=identity, kind="activate")
        self.assertEqual(item.state, "failed")
        self.assertEqual(item.error_code, "CAPABILITY_NO_LONGER_VALID")
