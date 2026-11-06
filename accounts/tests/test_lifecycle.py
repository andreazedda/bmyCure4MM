import io
from datetime import timedelta
from unittest.mock import patch

from django.contrib.auth import get_user_model
from django.core import mail
from django.core.exceptions import ValidationError
from django.core.management import call_command
from django.db import IntegrityError, transaction
from django.test import Client, TestCase, override_settings
from django.utils import timezone

from accounts import services
from accounts.forms import RequestForm
from accounts.models import AccessProfile as Identity
from accounts.models import Capability, Notification
from accounts.notifications import dispatch_notification

from .helpers import TEST_SETTINGS, LifecycleMixin, token_for


@override_settings(**TEST_SETTINGS)
class AccessLifecycleTests(LifecycleMixin, TestCase):
    def setUp(self):
        self.initialize()

    def test_owner_bootstrap_preserves_credentials_and_revocation(self):
        before = self.owner.user.password
        self.owner.access_state = Identity.AccessState.REVOKED
        self.owner.save()
        owner, created = services.bootstrap_owner()
        self.assertFalse(created)
        self.assertEqual(owner.user.password, before)
        self.assertEqual(owner.access_state, Identity.AccessState.REVOKED)
        self.assertEqual(Identity.objects.filter(is_owner=True).count(), 1)

    def test_owner_bootstrap_identity_conflicts_fail_closed(self):
        with override_settings(PLATFORM_OWNER_EMAIL="different@example.test"):
            with self.assertRaises(ValidationError):
                services.bootstrap_owner()
        self.assertEqual(get_user_model().objects.count(), 1)

    def test_bootstrap_does_not_promote_a_matching_legacy_user(self):
        with override_settings(
            PLATFORM_OWNER_USERNAME="unbound", PLATFORM_OWNER_EMAIL="unbound@example.test"
        ):
            self.owner.is_owner = False
            self.owner.save()
            user = get_user_model().objects.create_user("unbound", email="unbound@example.test")
            with self.assertRaises(ValidationError):
                services.bootstrap_owner()
            user.refresh_from_db()
            self.assertFalse(user.is_staff)

    def test_command_output_contains_no_owner_identity_or_token(self):
        output = io.StringIO()
        call_command("bootstrap_platform_owner", stdout=output)
        self.assertEqual(output.getvalue().strip(), "OWNER_BOOTSTRAP=PRESERVED")
        self.assertNotIn(self.owner.email, output.getvalue())

    def test_explicit_recovery_rotates_without_bootstrap_reactivating(self):
        self.sign_in(self.owner)
        call_command("recover_platform_owner", "--revoke-and-reissue", stdout=io.StringIO())
        self.assertEqual(self.client.get("/manage/access-requests/").status_code, 302)
        self.owner.refresh_from_db()
        self.assertEqual(self.owner.access_state, Identity.AccessState.PENDING)
        services.bootstrap_owner()
        self.assertEqual(
            Capability.objects.filter(
                identity=self.owner, purpose="owner", consumed_at=None, revoked_at=None
            ).count(),
            1,
        )

    def test_public_request_email_verify_review_activation_login(self):
        response = self.client.post(
            "/access/request/",
            {
                "email": " Alice+pilot@Example.Test ",
                "nickname": "Alice",
                "reason": "Explore synthetic models",
                "acknowledgement": "on",
            },
        )
        self.assertEqual(response.status_code, 302)
        identity = Identity.objects.get(email="alice+pilot@example.test")
        self.assertEqual(identity.request_state, Identity.RequestState.UNVERIFIED)
        self.assertIsNone(identity.user_id)
        self.assertFalse(Notification.objects.filter(kind="verified_request").exists())
        raw = token_for(identity, "verify")
        self.assertEqual(self.client.get("/access/verify/").status_code, 200)
        self.assertIsNone(Capability.objects.get(identity=identity, purpose="verify").consumed_at)
        self.assertEqual(self.client.post("/access/verify/", {"token": raw}).status_code, 302)
        identity.refresh_from_db()
        self.assertEqual(identity.request_state, Identity.RequestState.REVIEW)
        self.assertTrue(Notification.objects.filter(kind="verified_request").exists())
        services.decide_request(self.owner.user, identity.pk, approve=True)
        identity.refresh_from_db()
        self.assertFalse(identity.user.is_staff)
        self.assertFalse(identity.user.is_superuser)
        self.assertFalse(identity.user.is_active)
        self.assertFalse(identity.user.has_usable_password())
        response = self.client.post(
            "/access/activate/",
            {
                "token": token_for(identity, "activate"),
                "password": self.password,
                "password_confirm": self.password,
            },
        )
        self.assertEqual(response.status_code, 302)
        self.sign_in(identity)
        self.assertContains(self.client.get("/sim/getting-started/"), "getting started")
        self.assertEqual(self.client.post("/accounts/logout/").status_code, 302)
        self.assertEqual(self.client.get("/sim/getting-started/").status_code, 302)

    def test_public_form_cannot_claim_owner_or_privilege(self):
        for key in ("role", "is_staff", "is_superuser", "owner", "quota", "workspace_id", "user_id"):
            form = RequestForm(
                {
                    "email": "injected@example.test",
                    "nickname": "Test",
                    "reason": "Research",
                    "acknowledgement": True,
                    key: "1",
                }
            )
            self.assertFalse(form.is_valid(), key)
        self.assertIsNone(
            services.request_access(
                {"email": self.owner.email.upper(), "nickname": "Attacker", "reason": "Request"}
            )
        )
        self.assertEqual(Identity.objects.count(), 1)

    def test_email_normalization_preserves_dots_and_plus(self):
        self.assertEqual(services.normalize_email(" A.B+Pilot@Example.Test "), "a.b+pilot@example.test")
        self.assertNotEqual(
            services.normalize_email("a.b@example.test"), services.normalize_email("ab@example.test")
        )
        with transaction.atomic():
            get_user_model().objects.create_user("other", email="another@example.test")
        with self.assertRaises(IntegrityError), transaction.atomic():
            get_user_model().objects.create_user("duplicate", email=" Another@Example.Test ")

    def test_pending_unverified_rejected_cannot_login_or_approve(self):
        identity = services.request_access(
            {"email": "pending@example.test", "nickname": "Pending", "reason": "Test"}
        )
        with self.assertRaises(ValidationError):
            services.decide_request(self.owner.user, identity.pk, approve=True)
        services.confirm_capability(token_for(identity, "verify"), "verify")
        services.decide_request(self.owner.user, identity.pk, approve=False)
        identity.refresh_from_db()
        self.assertEqual(identity.request_state, Identity.RequestState.REJECTED)
        response = self.client.post("/accounts/login/", {"email": identity.email, "password": self.password})
        self.assertContains(response, "Unable to sign in")
        self.assertEqual(self.client.get("/sim/getting-started/").status_code, 302)

    def test_capability_wrong_purpose_expiry_and_replay(self):
        identity = services.request_access(
            {"email": "token@example.test", "nickname": "Token", "reason": "Test"}
        )
        raw = token_for(identity, "verify")
        with self.assertRaises(ValidationError):
            services.confirm_capability(raw, "activate", self.password)
        Capability.objects.filter(identity=identity).update(expires_at=timezone.now() - timedelta(seconds=1))
        with self.assertRaises(ValidationError):
            services.confirm_capability(raw, "verify")
        services.issue_capability(identity, "verify")
        raw = token_for(identity, "verify")
        services.confirm_capability(raw, "verify")
        with self.assertRaises(ValidationError):
            services.confirm_capability(raw, "verify")

    def test_activation_resend_revokes_old_link_and_deduplicates_approval(self):
        identity = self.make_member(activate=False)
        old = token_for(identity, "activate")
        count = get_user_model().objects.count()
        services.decide_request(self.owner.user, identity.pk, approve=True)
        self.assertEqual(get_user_model().objects.count(), count)
        services.resend_activation(self.owner.user, identity.pk)
        with self.assertRaises(ValidationError):
            services.confirm_capability(old, "activate", self.password)
        services.confirm_capability(token_for(identity, "activate"), "activate", self.password)

    def test_activation_rechecks_expiry_and_policy(self):
        identity = self.make_member(activate=False)
        raw = token_for(identity, "activate")
        self.policy.activations_open = False
        self.policy.save()
        with self.assertRaises(ValidationError):
            services.confirm_capability(raw, "activate", self.password)
        self.policy.activations_open = True
        self.policy.save()
        identity.expires_at = timezone.now() - timedelta(seconds=1)
        identity.save()
        with self.assertRaises(ValidationError):
            services.confirm_capability(raw, "activate", self.password)

    def test_revoked_suspended_expired_reset_cannot_reactivate(self):
        for index, state in enumerate(("REVOKED", "SUSPENDED", "EXPIRED")):
            identity = self.make_member(email=f"status{index}@example.test")
            services.issue_capability(identity, "reset")
            raw = token_for(identity, "reset")
            identity.access_state = state
            identity.save()
            with self.assertRaises(ValidationError):
                services.confirm_capability(raw, "reset", self.password)
            response = self.client.post("/password/reset/", {"email": identity.email})
            self.assertEqual(response.status_code, 200)
            self.assertEqual(Capability.objects.filter(identity=identity, purpose="reset").count(), 1)

    def test_status_requires_requester_session(self):
        identity = self.make_member()
        response = Client().get(f"/access/request/status/?email={identity.email}&id={identity.pk}")
        self.assertNotContains(response, identity.nickname)
        self.assertContains(response, "only in the browser")

    def test_owner_can_renew_suspended_access_but_cannot_restore_revocation(self):
        identity = self.make_member()
        self.sign_in(identity)
        services.change_access(self.owner.user, identity.pk, "SUSPENDED")
        self.assertEqual(self.client.get("/feedback/").status_code, 302)
        renewed = services.renew_access(self.owner.user, identity.pk, days=7)
        self.assertEqual(renewed.access_state, "ACTIVE")
        self.sign_in(renewed)
        self.assertEqual(self.client.get("/feedback/").status_code, 200)
        services.change_access(self.owner.user, identity.pk, "REVOKED")
        with self.assertRaises(ValidationError):
            services.renew_access(self.owner.user, identity.pk)

    def test_email_failure_and_console_retry_do_not_activate(self):
        identity = self.make_member(activate=False)
        Notification.objects.exclude(identity=identity, kind="activate").update(attempts=5)
        with patch("accounts.notifications.EmailMessage.send", side_effect=OSError("synthetic failure")):
            self.assertTrue(dispatch_notification())
        item = Notification.objects.get(identity=identity, kind="activate")
        self.assertEqual(item.state, "failed")
        self.assertEqual(item.error_code, "MAIL_BACKEND_FAILED")
        identity.refresh_from_db()
        self.assertEqual(identity.access_state, Identity.AccessState.PENDING)
        self.sign_in(self.owner)
        self.client.post("/manage/notifications/retry/", {"notification_id": item.pk})
        self.assertTrue(dispatch_notification())
        item.refresh_from_db()
        self.assertEqual(item.state, "sent")
        self.assertEqual(item.encrypted_payload, "")
        self.assertEqual(len(mail.outbox), 1)
        self.assertFalse(identity.user.is_active)

    def test_no_raw_token_persisted_or_echoed(self):
        identity = self.make_member(activate=False)
        raw = token_for(identity, "activate")
        self.assertFalse(Capability.objects.filter(digest=raw).exists())
        item = Notification.objects.get(identity=identity, kind="activate")
        self.assertNotIn(raw, item.encrypted_payload)
        response = self.client.post(
            "/access/activate/", {"token": raw, "password": "short", "password_confirm": "different"}
        )
        self.assertNotContains(response, raw)
        self.assertEqual(response["Referrer-Policy"], "same-origin")
        self.assertIn("no-store", response["Cache-Control"])

    def test_rate_budget_and_database_failure_fail_closed(self):
        for _ in range(5):
            self.assertTrue(
                services.rate_limit("synthetic", self.client.request().wsgi_request, "same", limit=5)
            )
        self.assertFalse(
            services.rate_limit("synthetic", self.client.request().wsgi_request, "same", limit=5)
        )
        from django.db import DatabaseError

        with patch("accounts.services.rate_limit", side_effect=DatabaseError):
            self.assertEqual(
                self.client.post(
                    "/accounts/login/", {"email": self.owner.email, "password": self.password}
                ).status_code,
                503,
            )
