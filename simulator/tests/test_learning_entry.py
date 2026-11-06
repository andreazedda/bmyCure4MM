from io import StringIO

from cryptography.fernet import Fernet
from django.core.management import call_command
from django.test import TestCase, override_settings

from accounts.tests.helpers import TEST_SETTINGS, LifecycleMixin
from simulator.models import Scenario


@override_settings(
    **TEST_SETTINGS,
    ACCOUNT_TOKEN_ENCRYPTION_KEY=Fernet.generate_key().decode(),
    EMAIL_BACKEND="django.core.mail.backends.locmem.EmailBackend",
)
class LearningEntryTests(LifecycleMixin, TestCase):
    def setUp(self):
        self.initialize()

    def test_empty_catalogue_initialization_preserves_existing_and_inactive_content(self):
        call_command("ensure_learning_example", stdout=StringIO())
        example = Scenario.objects.get()
        self.assertTrue(example.active)
        self.assertIsNone(example.patient_age)
        self.assertEqual(example.lab_snapshot, {})
        self.assertFalse(example.recommended_regimens.exists())
        example.active = False
        example.summary = "Owner-managed description"
        example.save()
        call_command("ensure_learning_example", stdout=StringIO())
        self.assertEqual(Scenario.objects.count(), 1)
        example.refresh_from_db()
        self.assertFalse(example.active)
        self.assertEqual(example.summary, "Owner-managed description")

    def test_visitor_and_registered_navigation_do_not_advertise_private_patient_routes(self):
        for path in ("/learn/", "/sim/", "/simulator/"):
            response = self.client.get(path)
            self.assertContains(response, "Request access")
            self.assertNotContains(response, 'href="/patients/"')
            self.assertNotContains(response, "Mario Rossi")
        self.sign_in(self.make_member())
        response = self.client.get("/learn/")
        self.assertContains(response, "Send a ticket")
        self.assertNotContains(response, "Use the admin")
        self.assertNotContains(response, 'href="/patients/"')

    def test_registration_cannot_discover_or_submit_legacy_assessment_context(self):
        from datetime import date

        from clinic.models import Assessment, Patient
        from simulator.models import SimulationAttempt

        patient = Patient.objects.create(
            mrn="DEMO-AUDIT-SYNTHETIC",
            first_name="SyntheticPrivateMarker",
            last_name="AuditOnly",
            birth_date=date(1980, 1, 1),
            diagnosis_date=date(2025, 1, 1),
            sex="O",
        )
        assessment = Assessment.objects.create(patient=patient, date=date(2025, 1, 1))
        call_command("ensure_learning_example", stdout=StringIO())
        scenario = Scenario.objects.get()
        for identity in (self.make_member(), self.owner):
            self.sign_in(identity)
            for prefix in ("/sim/", "/simulator/", "/learn/"):
                response = self.client.get(f"{prefix}{scenario.pk}/")
                self.assertEqual(response.status_code, 200)
                self.assertNotContains(response, "SyntheticPrivateMarker")
                self.assertNotContains(response, "Load demo data (DEMO)")
                self.assertNotContains(response, 'name="twin_assessment_id"')
                denied = self.client.post(
                    f"{prefix}scenarios/{scenario.pk}/simulate/",
                    {"use_twin": "on", "twin_assessment_id": assessment.pk},
                )
                self.assertEqual(denied.status_code, 400)
            self.client.post("/accounts/logout/")
        self.assertEqual(SimulationAttempt.objects.count(), 0)
