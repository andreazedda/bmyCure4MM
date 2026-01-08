from __future__ import annotations

import tempfile
from datetime import date

from django.contrib.auth import get_user_model
from django.test import Client, TestCase, override_settings
from django.urls import reverse

from clinic.models import Assessment, Patient
from simulator.models import Scenario, SimulationAttempt


class TwinPR1IntegrationTests(TestCase):
    def setUp(self) -> None:
        self.client = Client()
        User = get_user_model()
        self.user = User.objects.create_user(username="tester", password="pass")
        self.other_user = User.objects.create_user(username="other", password="pass")

        self.patient = Patient.objects.create(
            mrn="MRN-001",
            owner=self.user,
            first_name="Ada",
            last_name="Lovelace",
            birth_date=date(1970, 1, 1),
            sex="F",
            diagnosis_date=date(2020, 1, 1),
        )
        self.assessment = Assessment.objects.create(
            patient=self.patient,
            date=date(2025, 1, 1),
            r_iss="II",
            ldH_u_l=250,
            beta2m_mg_l=3.2,
            flc_ratio=2.0,
        )
        self.scenario = Scenario.objects.create(
            title="Test scenario",
            summary="A minimal scenario for integration testing.",
            active=True,
        )

        self.other_patient = Patient.objects.create(
            mrn="MRN-002",
            owner=self.other_user,
            first_name="Grace",
            last_name="Hopper",
            birth_date=date(1975, 1, 1),
            sex="F",
            diagnosis_date=date(2021, 1, 1),
        )
        self.other_assessment = Assessment.objects.create(
            patient=self.other_patient,
            date=date(2025, 2, 1),
            r_iss="III",
            ldH_u_l=600,
            beta2m_mg_l=10.0,
            flc_ratio=15.0,
        )

        self.null_owner_patient = Patient.objects.create(
            mrn="MRN-003",
            owner=None,
            first_name="Null",
            last_name="Owner",
            birth_date=date(1980, 1, 1),
            sex="M",
            diagnosis_date=date(2022, 1, 1),
        )
        self.null_owner_assessment = Assessment.objects.create(
            patient=self.null_owner_patient,
            date=date(2025, 3, 1),
            r_iss="I",
            ldH_u_l=200,
            beta2m_mg_l=2.2,
            flc_ratio=1.0,
        )

        self.demo_patient = Patient.objects.create(
            mrn="DEMO-TEST-001",
            owner=None,
            first_name="Demo",
            last_name="Patient",
            birth_date=date(1985, 1, 1),
            sex="F",
            diagnosis_date=date(2023, 1, 1),
        )
        self.demo_assessment = Assessment.objects.create(
            patient=self.demo_patient,
            date=date(2025, 4, 1),
            r_iss="II",
            ldH_u_l=300,
            beta2m_mg_l=4.0,
            flc_ratio=3.0,
        )

    def test_twin_preview_requires_login(self) -> None:
        url = reverse("simulator:twin_preview")
        response = self.client.get(url, {"id": self.assessment.pk})
        self.assertEqual(response.status_code, 302)

    def test_twin_preview_returns_payload_for_owner(self) -> None:
        self.client.login(username="tester", password="pass")
        url = reverse("simulator:twin_preview")
        response = self.client.get(url, {"id": self.assessment.pk})
        self.assertEqual(response.status_code, 200)
        data = response.json()
        self.assertIn("assessment", data)
        self.assertIn("inputs", data["assessment"])
        self.assertIn("ldh_u_l", data["assessment"]["inputs"])
        self.assertIn("twin", data)
        self.assertIn("risk_score", data["twin"])

    def test_twin_preview_denies_non_owner(self) -> None:
        self.client.login(username="tester", password="pass")
        url = reverse("simulator:twin_preview")
        response = self.client.get(url, {"id": self.other_assessment.pk})
        self.assertEqual(response.status_code, 404)

    def test_twin_preview_hides_null_owner_from_non_staff(self) -> None:
        self.client.login(username="tester", password="pass")
        url = reverse("simulator:twin_preview")
        response = self.client.get(url, {"id": self.null_owner_assessment.pk})
        self.assertEqual(response.status_code, 404)

    def test_twin_preview_allows_demo_for_non_staff(self) -> None:
        self.client.login(username="tester", password="pass")
        url = reverse("simulator:twin_preview")
        response = self.client.get(url, {"id": self.demo_assessment.pk})
        self.assertEqual(response.status_code, 200)
        data = response.json()
        self.assertIn("twin", data)
        self.assertIn("risk_score", data["twin"])

    def test_twin_preview_allows_null_owner_for_staff(self) -> None:
        User = get_user_model()
        staff = User.objects.create_user(username="staff", password="pass", is_staff=True)
        self.client.login(username="staff", password="pass")
        url = reverse("simulator:twin_preview")
        response = self.client.get(url, {"id": self.null_owner_assessment.pk})
        self.assertEqual(response.status_code, 200)

    @override_settings(PREDLAB_V2=True)
    def test_simulate_scenario_writes_twin_artifact(self) -> None:
        self.client.login(username="tester", password="pass")

        with tempfile.TemporaryDirectory() as tmpdir:
            with override_settings(MEDIA_ROOT=tmpdir, MEDIA_URL="/media/"):
                url = reverse("simulator:scenario_simulate", args=[self.scenario.pk])
                payload = {
                    "baseline_tumor_cells": 1.0e9,
                    "baseline_healthy_cells": 5.0e11,
                    "lenalidomide_dose": 25.0,
                    "bortezomib_dose": 1.3,
                    "daratumumab_dose": 16.0,
                    "time_horizon": 30,
                    "tumor_growth_rate": 0.023,
                    "healthy_growth_rate": 0.015,
                    "interaction_strength": 0.05,
                    "preset": "VRd",
                    "creatinine_clearance": 90.0,
                    "neuropathy_grade": 0,
                    "anc": 2.0,
                    "platelets": 150,
                    "cohort_size": 1,
                    "use_twin": "on",
                    "twin_assessment_id": str(self.assessment.pk),
                    "twin_biology_mode": "auto",
                }
                response = self.client.post(url, data=payload)
                self.assertEqual(response.status_code, 200)

                attempt = SimulationAttempt.objects.filter(scenario=self.scenario, user=self.user).latest(
                    "submitted"
                )
                self.assertIn("twin_params.json", attempt.results)

    def test_resolve_solver_inputs_strips_auto(self) -> None:
        resolved = {
            "baseline_tumor_cells": "auto",
            "baseline_healthy_cells": "AUTO",
            "time_horizon": "auto",
            "tumor_growth_rate": "auto",
            "healthy_growth_rate": "auto",
            "interaction_strength": "auto",
            "immune_compromise_index": "auto",
            "carrying_capacity_tumor": "auto",
            "carrying_capacity_healthy": "auto",
            "lenalidomide_dose": "auto",
            "bortezomib_dose": "auto",
            "daratumumab_dose": "auto",
        }
        out = SimulationAttempt._resolve_solver_inputs(resolved)
        self.assertIsInstance(out["baseline_tumor_cells"], float)
        self.assertIsInstance(out["growth_rates"]["tumor"], float)
        self.assertTrue(all(isinstance(v, float) for v in out["drug_doses"].values()))

    def test_resolve_solver_inputs_raises_actionable_field_error(self) -> None:
        resolved = {
            "baseline_tumor_cells": 1.0,
            "baseline_healthy_cells": 1.0,
            "time_horizon": 1.0,
            "tumor_growth_rate": "not-a-number",
            "healthy_growth_rate": 0.1,
            "interaction_strength": 0.1,
            "immune_compromise_index": 1.0,
            "carrying_capacity_tumor": None,
            "carrying_capacity_healthy": None,
            "lenalidomide_dose": 1.0,
            "bortezomib_dose": 1.0,
            "daratumumab_dose": 1.0,
        }
        with self.assertRaises(ValueError) as ctx:
            SimulationAttempt._resolve_solver_inputs(resolved)
        self.assertIn("growth_rates.tumor", str(ctx.exception))
