from __future__ import annotations

import tempfile
from datetime import date

from django.contrib.auth import get_user_model
from django.test import TestCase, override_settings
from django.urls import reverse

from clinic.models import Assessment, Patient
from simulator.models import Scenario


class TwinEngineRegressionTests(TestCase):
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("regression", password="pass1234")
        self.client.force_login(self.user)
        self.patient = Patient.objects.create(
            mrn="MM-REG-001",
            owner=self.user,
            first_name="Regression",
            last_name="Patient",
            birth_date=date(1970, 1, 1),
            sex="M",
            diagnosis_date=date(2024, 1, 1),
        )
        self.assessment = Assessment.objects.create(
            patient=self.patient,
            date=date(2025, 1, 1),
            m_protein_g_dl=1.0,
            flc_ratio=2.0,
            hemoglobin_g_dl=12.0,
            r_iss="II",
            beta2m_mg_l=3.0,
            ldH_u_l=220,
        )
        self.scenario = Scenario.objects.create(
            title="Regression scenario",
            summary="Existing educational scenario path should keep working.",
            active=True,
        )

    def test_existing_scenario_simulation_still_works(self) -> None:
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
        with tempfile.TemporaryDirectory() as tmpdir:
            with override_settings(MEDIA_ROOT=tmpdir, MEDIA_URL="/media/"):
                response = self.client.post(reverse("simulator:scenario_simulate", args=[self.scenario.pk]), data=payload)
        self.assertEqual(response.status_code, 200)
        self.assertIn("Simulation Results", response.content.decode())

    def test_existing_twin_preview_still_works(self) -> None:
        response = self.client.get(reverse("simulator:twin_preview"), {"id": self.assessment.pk})
        self.assertEqual(response.status_code, 200)
        self.assertIn("risk_score", response.json()["twin"])

    def test_existing_clinic_patient_detail_still_works(self) -> None:
        response = self.client.get(reverse("clinic:patient_detail", args=[self.patient.pk]))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Latest simulation")
        self.assertContains(response, "Research Twin / What-if Simulation")
