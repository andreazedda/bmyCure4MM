"""Account UX regressions using the real lifecycle and ordinary password login."""

import uuid

from django.core.exceptions import ValidationError
from django.test import TestCase, override_settings
from django.utils import timezone

from accounts import services
from accounts.models import AccessProfile, SupportTicket

from .helpers import TEST_SETTINGS, LifecycleMixin


@override_settings(**TEST_SETTINGS)
class AccountRecoveryUXTests(LifecycleMixin, TestCase):
    def setUp(self):
        self.initialize()
        self.alice = self.make_member()

    def test_invalid_reply_preserves_text_and_does_not_publish_it(self):
        ticket = SupportTicket.objects.create(
            identity=self.alice, category="bug", message="Synthetic report", submission_key=uuid.uuid4()
        )
        self.sign_in(self.owner)
        draft = "<script>draft only</script>" + "x" * 2001
        response = self.client.post(
            "/manage/feedback/",
            {
                "ticket_id": ticket.pk,
                "category": "bug",
                "status": "answered",
                "reply": draft,
                "internal_note": "Keep this draft note",
            },
        )
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "&lt;script&gt;draft only&lt;/script&gt;")
        self.assertContains(response, "Keep this draft note")
        self.assertContains(response, "data-error-summary")
        ticket.refresh_from_db()
        self.assertEqual(ticket.reply, "")
        self.assertEqual(ticket.status, "open")
        self.client.post("/accounts/logout/")
        self.sign_in(self.alice)
        self.assertNotContains(self.client.get("/feedback/"), "Keep this draft note")

    def test_older_tickets_and_requests_remain_reachable(self):
        SupportTicket.objects.bulk_create(
            [
                SupportTicket(
                    identity=self.alice,
                    category="bug",
                    message=f"Synthetic ticket {i}",
                    submission_key=uuid.uuid4(),
                )
                for i in range(23)
            ]
        )
        self.sign_in(self.alice)
        first = self.client.get("/feedback/")
        second = self.client.get("/feedback/?page=2")
        self.assertContains(first, "?page=2")
        self.assertEqual(len(first.context["tickets"]), 20)
        self.assertEqual(len(second.context["tickets"]), 3)
        self.assertEqual(
            set(t.pk for t in first.context["tickets"]) | set(t.pk for t in second.context["tickets"]),
            set(SupportTicket.objects.values_list("pk", flat=True)),
        )
        self.client.post("/accounts/logout/")
        self.sign_in(self.owner)
        AccessProfile.objects.bulk_create(
            [
                AccessProfile(
                    email=f"pending{i}@example.test",
                    nickname=f"Synthetic pending {i}",
                    verified_at=timezone.now(),
                    request_state="PENDING_REVIEW",
                )
                for i in range(23)
            ]
        )
        page = self.client.get("/manage/access-requests/?page=2")
        self.assertGreater(len(page.context["rows"]), 0)
        self.assertEqual(page.context["rows"].paginator.count, 24)
        self.assertEqual(len(self.client.get("/manage/feedback/?page=2").context["tickets"]), 3)

    def test_revoke_requires_a_second_explicit_post(self):
        self.sign_in(self.owner)
        response = self.client.post(
            "/manage/access-requests/", {"identity_id": self.alice.pk, "action": "revoke"}
        )
        self.assertContains(response, "Confirm revocation")
        self.alice.refresh_from_db()
        self.assertEqual(self.alice.access_state, "ACTIVE")
        self.client.get("/manage/access-requests/")
        self.alice.refresh_from_db()
        self.assertEqual(self.alice.access_state, "ACTIVE")
        self.assertEqual(
            self.client.post(
                "/manage/access-requests/",
                {"identity_id": self.alice.pk, "action": "revoke", "confirm_revoke": "true"},
            ).status_code,
            302,
        )
        self.alice.refresh_from_db()
        self.assertEqual(self.alice.access_state, "REVOKED")
        with self.assertRaises(ValidationError):
            services.change_access(self.owner.user, self.alice.pk, "SUSPENDED")
        with self.assertRaises(ValidationError):
            services.renew_access(self.owner.user, self.alice.pk, 14)
        self.alice.refresh_from_db()
        self.assertEqual(self.alice.access_state, "REVOKED")

    def test_rate_limit_has_a_recoverable_localized_page(self):
        for _ in range(5):
            self.client.post("/access/request/", {"email": "invalid"})
        response = self.client.post("/access/request/", {"email": "invalid"})
        self.assertEqual(response.status_code, 429)
        self.assertEqual(response["Retry-After"], "3600")
        self.assertContains(response, "Troppi tentativi", status_code=429)
        self.assertContains(response, "/accounts/login/", status_code=429)

    def test_learning_form_description_references_resolve(self):
        from html.parser import HTMLParser

        from django.urls import reverse

        from simulator.models import Scenario

        scenario = Scenario.objects.create(
            title="Synthetic accessibility check",
            clinical_stage="newly_diagnosed",
            summary="Independent synthetic form test",
            risk_stratification="Standard",
        )
        self.sign_in(self.alice)
        response = self.client.get(reverse("simulator:scenario_detail", args=[scenario.pk]))
        self.assertEqual(response.status_code, 200)

        class References(HTMLParser):
            def __init__(self):
                super().__init__()
                self.ids = set()
                self.references = []

            def handle_starttag(self, tag, attrs):
                values = dict(attrs)
                if values.get("id"):
                    self.ids.add(values["id"])
                for reference in values.get("aria-describedby", "").split():
                    self.references.append((values.get("id"), reference))

        document = References()
        document.feed(response.content.decode())
        self.assertEqual([item for item in document.references if item[1] not in document.ids], [])
