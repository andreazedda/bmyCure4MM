from __future__ import annotations

from datetime import date

from django.contrib.auth import get_user_model
from django.test import TestCase

from clinic.models import Assessment, Patient
from twin_engine.sensitivity import run_counterfactual_sensitivity
from twin_engine.state_model import initialize_from_assessment


class SensitivityDiagnosticsTests(TestCase):
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("sensitivity-user", password="pass1234")
        self.patient = Patient.objects.create(
            mrn="MM-SENS-001",
            owner=self.user,
            first_name="Synthetic",
            last_name="Sensitivity",
            birth_date=date(1966, 4, 4),
            sex="F",
            diagnosis_date=date(2024, 1, 1),
        )
        self.assessment = Assessment.objects.create(
            patient=self.patient,
            date=date(2025, 1, 1),
            m_protein_g_dl=1.1,
            flc_ratio=2.1,
            hemoglobin_g_dl=11.7,
            r_iss="II",
            beta2m_mg_l=3.0,
            ldH_u_l=220,
        )
        self.state = initialize_from_assessment(self.assessment, user=self.user)

    def test_sensitivity_returns_ranked_driver_rows(self) -> None:
        intervention = {
            "label": "LEN_5MG_DAILY_14D",
            "classification": "mechanistic_simulation_only",
            "intervention": {
                "drug": "lenalidomide",
                "dose_mg": 5.0,
                "schedule": {"type": "daily"},
                "duration_days": 14,
                "start_day": 0,
            },
        }
        result = run_counterfactual_sensitivity(self.patient, self.state, intervention, 14)
        self.assertEqual(result["status"], "completed")
        self.assertGreater(len(result["rows"]), 0)
        self.assertIn("rank", result["rows"][0])
        self.assertIn("research_utility_v2", result["baseline_metrics"])
