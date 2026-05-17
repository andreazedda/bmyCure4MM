from __future__ import annotations

from datetime import date

from django.contrib.auth import get_user_model
from django.test import TestCase

from clinic.models import Assessment, Patient
from twin_engine.backtesting import run_patient_backtest
from twin_engine.state_model import initialize_from_assessment


class BacktestingDiagnosticsTests(TestCase):
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("backtest-user", password="pass1234")
        self.patient = Patient.objects.create(
            mrn="MM-BT-001",
            owner=self.user,
            first_name="Synthetic",
            last_name="Backtest",
            birth_date=date(1966, 2, 2),
            sex="F",
            diagnosis_date=date(2024, 1, 1),
        )
        self.assessments = [
            Assessment.objects.create(patient=self.patient, date=date(2025, 1, 1), m_protein_g_dl=1.2, flc_ratio=2.4, hemoglobin_g_dl=11.8),
            Assessment.objects.create(patient=self.patient, date=date(2025, 1, 15), m_protein_g_dl=1.1, flc_ratio=2.2, hemoglobin_g_dl=11.9),
            Assessment.objects.create(patient=self.patient, date=date(2025, 1, 29), m_protein_g_dl=0.95, flc_ratio=2.0, hemoglobin_g_dl=12.0),
        ]
        initialize_from_assessment(self.assessments[0], user=self.user)

    def test_backtest_requires_minimum_history(self) -> None:
        patient = Patient.objects.create(
            mrn="MM-BT-002",
            owner=self.user,
            first_name="Synthetic",
            last_name="Short",
            birth_date=date(1968, 3, 3),
            sex="M",
            diagnosis_date=date(2024, 1, 1),
        )
        Assessment.objects.create(patient=patient, date=date(2025, 1, 1), m_protein_g_dl=1.0, flc_ratio=2.0, hemoglobin_g_dl=12.0)
        Assessment.objects.create(patient=patient, date=date(2025, 1, 15), m_protein_g_dl=0.9, flc_ratio=1.9, hemoglobin_g_dl=12.1)
        result = run_patient_backtest(patient)
        self.assertEqual(result["status"], "insufficient_history")

    def test_backtest_returns_fold_rows_and_aggregate_metrics(self) -> None:
        result = run_patient_backtest(self.patient)
        self.assertEqual(result["status"], "completed")
        self.assertEqual(result["n_folds"], 1)
        self.assertEqual(result["fold_rows"][0]["train_end_date"], "2025-01-15")
        self.assertIn("m_protein_g_dl", result["by_biomarker"])
