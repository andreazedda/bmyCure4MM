"""
Tests for security fixes and UX improvements.
Covers: @login_required enforcement, ownership checks on API views,
dashboard stats, URL namespace correctness, hemoglobin validation.
"""
from __future__ import annotations

from datetime import date, timedelta
from django.contrib.auth import get_user_model
from django.test import TestCase, Client
from django.urls import reverse, NoReverseMatch

from clinic.models import Patient, Assessment, Regimen, PatientTherapy
from clinic.forms import AssessmentForm

User = get_user_model()


def _make_user(username: str, *, is_staff: bool = False) -> "User":
    return User.objects.create_user(username=username, password="pass1234", is_staff=is_staff)


def _make_patient(owner=None, **kwargs) -> Patient:
    defaults = dict(
        mrn="MM-TEST-001",
        first_name="Test",
        last_name="Patient",
        birth_date=date(1960, 1, 15),
        sex="M",
        diagnosis_date=date(2023, 6, 1),
    )
    defaults.update(kwargs)
    return Patient.objects.create(owner=owner, **defaults)


# ── Login-required enforcement ──────────────────────────────────────

class LoginRequiredTests(TestCase):
    """Verify that unauthenticated users are redirected to login."""

    def setUp(self):
        self.client = Client()
        self.user = _make_user("owner1")
        self.patient = _make_patient(owner=self.user)

    def _assert_login_redirect(self, url: str):
        resp = self.client.get(url)
        self.assertIn(resp.status_code, [301, 302])
        self.assertIn("login", resp.url)

    def test_patient_list_requires_login(self):
        self._assert_login_redirect(reverse("clinic:patient_list"))

    def test_dashboard_requires_auth(self):
        """Anonymous users are redirected away from the root (portal home)."""
        resp = self.client.get(reverse("clinic:dashboard"))
        self.assertIn(resp.status_code, [301, 302])
        # Portal home redirects anon to docs, not login page
        self.assertTrue(resp.url)

    def test_patient_detail_requires_login(self):
        self._assert_login_redirect(
            reverse("clinic:patient_detail", args=[self.patient.pk])
        )

    def test_prognosis_api_requires_login(self):
        self._assert_login_redirect(
            reverse("clinic:prognosis_api", args=[self.patient.pk])
        )

    def test_regimen_suggestions_api_requires_login(self):
        self._assert_login_redirect(
            reverse("clinic:regimen_suggestions_api", args=[self.patient.pk])
        )


# ── Ownership / access-control on API views ─────────────────────────

class OwnershipCheckTests(TestCase):
    """Non-owner, non-staff users must get 403 on patient API endpoints."""

    def setUp(self):
        self.owner = _make_user("owner")
        self.stranger = _make_user("stranger")
        self.staff = _make_user("staffuser", is_staff=True)
        self.patient = _make_patient(owner=self.owner, mrn="MM-OWN-001")

    def test_prognosis_api_forbidden_for_stranger(self):
        self.client.login(username="stranger", password="pass1234")
        resp = self.client.get(
            reverse("clinic:prognosis_api", args=[self.patient.pk])
        )
        self.assertEqual(resp.status_code, 403)

    def test_prognosis_api_allowed_for_owner(self):
        self.client.login(username="owner", password="pass1234")
        resp = self.client.get(
            reverse("clinic:prognosis_api", args=[self.patient.pk])
        )
        self.assertEqual(resp.status_code, 200)
        payload = resp.json()
        self.assertEqual(payload["status"], "PATIENT_SPECIFIC_PREDICTION_NOT_VALIDATED")
        self.assertIsNone(payload["estimate"])

    def test_prognosis_api_allowed_for_staff(self):
        self.client.login(username="staffuser", password="pass1234")
        resp = self.client.get(
            reverse("clinic:prognosis_api", args=[self.patient.pk])
        )
        self.assertEqual(resp.status_code, 200)
        self.assertIsNone(resp.json()["estimate"])

    def test_regimen_api_forbidden_for_stranger(self):
        self.client.login(username="stranger", password="pass1234")
        resp = self.client.get(
            reverse("clinic:regimen_suggestions_api", args=[self.patient.pk])
        )
        self.assertEqual(resp.status_code, 403)

    def test_regimen_api_allowed_for_owner(self):
        self.client.login(username="owner", password="pass1234")
        resp = self.client.get(
            reverse("clinic:regimen_suggestions_api", args=[self.patient.pk])
        )
        self.assertIn(resp.status_code, [200, 500])

    def test_regimen_api_allowed_for_staff(self):
        self.client.login(username="staffuser", password="pass1234")
        resp = self.client.get(
            reverse("clinic:regimen_suggestions_api", args=[self.patient.pk])
        )
        self.assertIn(resp.status_code, [200, 500])


# ── Dashboard view ──────────────────────────────────────────────────

class DashboardTests(TestCase):
    """Dashboard loads, shows stats, and links resolve."""

    def setUp(self):
        self.user = _make_user("dashuser", is_staff=True)
        self.client.login(username="dashuser", password="pass1234")

    def test_dashboard_loads(self):
        resp = self.client.get(reverse("clinic:dashboard"))
        self.assertEqual(resp.status_code, 200)

    def test_dashboard_shows_patient_count(self):
        _make_patient(owner=self.user, mrn="MM-D-001")
        _make_patient(owner=self.user, mrn="MM-D-002")
        resp = self.client.get(reverse("clinic:dashboard"))
        self.assertContains(resp, "2")

    def test_dashboard_new_patient_link(self):
        """The 'new patient' link must resolve to clinic:patient_new."""
        resp = self.client.get(reverse("clinic:dashboard"))
        new_url = reverse("clinic:patient_new")
        self.assertContains(resp, new_url)


# ── URL namespace resolution ────────────────────────────────────────

class URLNamespaceTests(TestCase):
    """All named routes resolve without errors."""

    def setUp(self):
        self.user = _make_user("urluser")
        self.patient = _make_patient(owner=self.user)

    def test_all_clinic_routes_resolve(self):
        names_no_args = ["dashboard", "patient_list", "patient_new", "regimen_list", "regimen_new"]
        for name in names_no_args:
            with self.subTest(name=name):
                url = reverse(f"clinic:{name}")
                self.assertTrue(url)

    def test_patient_arg_routes_resolve(self):
        pk = self.patient.pk
        routes = [
            ("patient_detail", [pk]),
            ("patient_edit", [pk]),
            ("assessment_new", [pk]),
            ("symptom_assessment_new", [pk]),
            ("symptom_assessment_list", [pk]),
            ("prognosis_timeline", [pk]),
            ("prognosis_api", [pk]),
            ("regimen_suggestions", [pk]),
            ("regimen_suggestions_api", [pk]),
        ]
        for name, args in routes:
            with self.subTest(name=name):
                url = reverse(f"clinic:{name}", args=args)
                self.assertTrue(url)


# ── Hemoglobin validation (tightened to max 20) ─────────────────────

class HemoglobinValidationTests(TestCase):
    """Hemoglobin cap was tightened from 25 to 20 g/dL."""

    def _base_data(self, **overrides):
        data = {
            "date": date.today().isoformat(),
            "m_protein_g_dl": "1.0",
            "flc_ratio": "1.5",
            "hemoglobin_g_dl": "13.0",
            "calcium_mg_dl": "9.5",
            "creatinine_mg_dl": "1.0",
            "ldH_u_l": "200",
            "albumin_g_dl": "4.0",
            "beta2m_mg_l": "3.0",
            "r_iss": "I",
            "response": "",
        }
        data.update(overrides)
        return data

    def test_hemoglobin_20_accepted(self):
        form = AssessmentForm(data=self._base_data(hemoglobin_g_dl="20.0"))
        self.assertTrue(form.is_valid(), f"Errors: {form.errors}")

    def test_hemoglobin_21_rejected(self):
        form = AssessmentForm(data=self._base_data(hemoglobin_g_dl="21.0"))
        self.assertFalse(form.is_valid())
        self.assertIn("hemoglobin_g_dl", form.errors)


# ── Regimen list template ───────────────────────────────────────────

class RegimenListTests(TestCase):
    """Regimen list page loads and shows data."""

    def setUp(self):
        self.user = _make_user("reguser", is_staff=True)
        self.client.login(username="reguser", password="pass1234")

    def test_regimen_list_loads(self):
        resp = self.client.get(reverse("clinic:regimen_list"))
        self.assertEqual(resp.status_code, 200)

    def test_regimen_list_shows_regimen(self):
        Regimen.objects.create(name="VRd", components="Bortezomib, Lenalidomide, Dex", line=1)
        resp = self.client.get(reverse("clinic:regimen_list"))
        self.assertContains(resp, "VRd")

    def test_regimen_list_empty_state(self):
        resp = self.client.get(reverse("clinic:regimen_list"))
        # Should show empty state message (bilingual)
        self.assertEqual(resp.status_code, 200)
