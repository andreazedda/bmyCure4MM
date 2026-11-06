import base64
import uuid

from django.test import Client, TestCase, override_settings

from accounts import services
from accounts.models import SupportTicket

from .helpers import TEST_SETTINGS, LifecycleMixin


@override_settings(**TEST_SETTINGS)
class ExistingSiteTests(LifecycleMixin, TestCase):
    def setUp(self):
        self.initialize()
        self.alice = self.make_member()
        self.bob = self.make_member("bob@example.test")

    def test_login_opens_the_existing_guide_and_navigation(self):
        self.sign_in(self.alice)
        self.assertContains(self.client.get("/sim/getting-started/"), "Explore the learning scenarios")
        self.assertContains(self.client.get("/sim/"), "Suggestions / tickets")
        self.assertFalse(self.alice.user.is_staff or self.alice.user.is_superuser)

    def test_registration_does_not_grant_clinical_api_or_basic_auth_access(self):
        self.sign_in(self.alice)
        for path in ("/patients/", "/patients/1/", "/research/", "/api/patients/", "/simulator/manage/"):
            self.assertIn(self.client.get(path).status_code, (403, 404), path)
        basic = base64.b64encode(f"{self.alice.user.username}:{self.password}".encode()).decode()
        self.assertIn(
            Client().get("/api/patients/", HTTP_AUTHORIZATION=f"Basic {basic}").status_code, (401, 403)
        )

    def test_ticket_reply_is_private_and_survives_relogin(self):
        self.sign_in(self.alice)
        key = uuid.uuid4()
        data = {
            "category": "feature",
            "message": "Please clarify how to start using the platform.",
            "submission_key": key,
        }
        self.assertEqual(self.client.post("/feedback/", data).status_code, 302)
        self.assertEqual(self.client.post("/feedback/", data).status_code, 302)
        self.assertEqual(SupportTicket.objects.count(), 1)
        ticket = SupportTicket.objects.get()
        self.client.post("/accounts/logout/")
        self.sign_in(self.owner)
        response = self.client.post(
            "/manage/feedback/",
            {
                "ticket_id": ticket.pk,
                "category": "feature",
                "status": "answered",
                "reply": "Suggestion received and reviewed.",
                "internal_note": "Owner internal note",
            },
        )
        self.assertEqual(response.status_code, 302)
        self.client.post("/accounts/logout/")
        self.sign_in(self.alice)
        response = self.client.get("/feedback/")
        self.assertContains(response, "Suggestion received and reviewed.")
        self.assertNotContains(response, "Owner internal note")
        self.client.post("/accounts/logout/")
        self.sign_in(self.bob)
        self.assertNotContains(self.client.get("/feedback/"), ticket.message)

    def test_owner_review_requires_authenticator_confirmation(self):
        self.client.post("/accounts/login/", {"email": self.owner.email, "password": self.password})
        self.assertRedirects(self.client.get("/manage/access-requests/"), "/accounts/security/")
        self.assertEqual(self.client.post("/accounts/security/", {"code": "000000"}).status_code, 200)

    def test_owner_governance_does_not_grant_implicit_clinical_access(self):
        self.sign_in(self.owner)
        self.assertEqual(self.client.get("/manage/access-requests/").status_code, 200)
        for path in ("/patients/", "/research/", "/api/patients/"):
            self.assertIn(self.client.get(path).status_code, (403, 404), path)

    def test_documentation_navigation_is_not_shared_between_accounts(self):
        self.alice.nickname = "PrivateAliceNavigation"
        self.alice.save()
        self.sign_in(self.alice)
        first = self.client.get("/docs/")
        self.assertContains(first, self.alice.nickname)
        self.assertIn("no-store", first["Cache-Control"])
        self.client.post("/accounts/logout/")
        self.sign_in(self.bob)
        self.assertNotContains(self.client.get("/docs/"), self.alice.nickname)

    def test_registered_user_runs_existing_engine_and_cannot_export_another_users_run(self):
        from django.urls import reverse

        from simulator.models import Scenario, SimulationAttempt

        scenario = Scenario.objects.create(
            title="Synthetic registration integration example",
            clinical_stage="newly_diagnosed",
            summary="Independent synthetic software test",
            risk_stratification="Standard",
        )
        self.sign_in(self.alice)
        response = self.client.post(
            reverse("simulator:scenario_simulate", args=[scenario.pk]),
            {
                "baseline_tumor_cells": 1e9,
                "baseline_healthy_cells": 5e11,
                "lenalidomide_dose": 25,
                "bortezomib_dose": 1.3,
                "daratumumab_dose": 16,
                "time_horizon": 60,
                "tumor_growth_rate": 0.023,
                "healthy_growth_rate": 0.015,
                "interaction_strength": 0.05,
            },
        )
        self.assertContains(response, "Simulated outputs")
        attempt = SimulationAttempt.objects.get(user=self.alice.user)
        self.assertIn("tumor_reduction", attempt.results_summary)
        export = reverse("simulator:attempt_export", args=[attempt.pk])
        self.assertEqual(self.client.get(export).status_code, 200)
        csv_response = self.client.get(export + "?format=csv")
        self.assertEqual(csv_response.status_code, 200)
        self.assertEqual(csv_response["Content-Type"], "text/csv")
        self.assertIn("no-store", csv_response["Cache-Control"])
        self.assertIn(b"time", b"".join(csv_response.streaming_content))
        self.assertEqual(
            self.client.get(reverse("simulator:attempt_plot_embed", args=[attempt.pk])).status_code, 200
        )
        self.assertContains(response, reverse("simulator:attempt_plot_embed", args=[attempt.pk]))
        self.client.post("/accounts/logout/")
        self.sign_in(self.bob)
        self.assertEqual(self.client.get(export).status_code, 403)
        self.assertEqual(self.client.get(export + "?format=csv").status_code, 403)
        self.assertEqual(
            self.client.get(reverse("simulator:attempt_plot_embed", args=[attempt.pk])).status_code, 403
        )
        self.assertNotContains(
            self.client.get(reverse("simulator:scenario_detail", args=[scenario.pk])), "1 simulation run(s)"
        )

        self.client.post("/accounts/logout/")
        self.sign_in(self.alice)
        services.change_access(self.owner.user, self.alice.pk, "SUSPENDED")
        self.assertEqual(self.client.get(export + "?format=csv").status_code, 302)
        self.assertEqual(
            self.client.get(reverse("simulator:attempt_plot_embed", args=[attempt.pk])).status_code, 302
        )
