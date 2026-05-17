from __future__ import annotations

from datetime import date
from types import SimpleNamespace

from django.contrib.auth import get_user_model
from django.test import TestCase

from clinic.models import Assessment, Patient
from twin_engine.state_model import initialize_from_assessment
from twin_engine.uncertainty import (
    HEURISTIC_PARAMETER_UNCERTAINTY_SOURCE,
    UncertaintyConfig,
    build_parameter_uncertainty_space,
    classify_uncertainty_width,
    compute_interval,
    run_counterfactual_uncertainty,
    sample_parameter_sets,
    summarize_distribution,
)


class UncertaintyDiagnosticsTests(TestCase):
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("uncertainty-user", password="pass1234")
        self.patient = Patient.objects.create(
            mrn="MM-UNC-001",
            owner=self.user,
            first_name="Synthetic",
            last_name="Uncertainty",
            birth_date=date(1968, 1, 1),
            sex="F",
            diagnosis_date=date(2024, 1, 1),
        )
        self.assessment = Assessment.objects.create(
            patient=self.patient,
            date=date(2025, 1, 1),
            m_protein_g_dl=1.1,
            flc_ratio=2.1,
            hemoglobin_g_dl=11.8,
            r_iss="II",
            beta2m_mg_l=3.1,
            ldH_u_l=215,
        )

    def test_sample_parameter_sets_are_reproducible_with_seed(self) -> None:
        twin_state = SimpleNamespace(
            parameters={
                "tumor_growth_rate": 0.023,
                "healthy_growth_rate": 0.015,
                "immune_compromise_index": 1.1,
                "observation": {"alpha_M": 1.0e-9},
            },
            parameter_uncertainty={},
        )
        config = UncertaintyConfig(n_samples=3, random_seed=123)
        sample_a = sample_parameter_sets(twin_state, config)
        sample_b = sample_parameter_sets(twin_state, config)
        self.assertEqual(sample_a["samples"], sample_b["samples"])

    def test_compute_interval_and_summary_distribution(self) -> None:
        distribution = summarize_distribution([1.0, 2.0, 3.0, 4.0])
        interval = compute_interval([1.0, 2.0, 3.0, 4.0], lower=0.25, upper=0.75)
        self.assertEqual(distribution["count"], 4)
        self.assertAlmostEqual(distribution["median"], 2.5)
        self.assertAlmostEqual(interval["lower"], 1.75)
        self.assertAlmostEqual(interval["upper"], 3.25)
        self.assertAlmostEqual(interval["width"], 1.5)

    def test_classify_uncertainty_width_works(self) -> None:
        self.assertEqual(classify_uncertainty_width("tumor_reduction", 0.05, 1.0), "narrow")
        self.assertEqual(classify_uncertainty_width("tumor_reduction", 0.25, 1.0), "moderate")
        self.assertEqual(classify_uncertainty_width("tumor_reduction", 0.75, 1.0), "wide")

    def test_heuristic_uncertainty_source_is_used_without_calibrated_distribution(self) -> None:
        twin_state = SimpleNamespace(
            parameters={"tumor_growth_rate": 0.023, "immune_compromise_index": 1.1},
            parameter_uncertainty={},
        )
        uncertainty_space = build_parameter_uncertainty_space(twin_state)
        self.assertEqual(uncertainty_space["parameter_uncertainty_source"], HEURISTIC_PARAMETER_UNCERTAINTY_SOURCE)

    def test_run_counterfactual_uncertainty_returns_metric_summaries(self) -> None:
        twin_state = initialize_from_assessment(self.assessment, user=self.user)
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
        result = run_counterfactual_uncertainty(
            self.patient,
            twin_state,
            intervention,
            14,
            UncertaintyConfig(n_samples=4, random_seed=42),
        )
        self.assertEqual(result["status"], "completed")
        self.assertEqual(result["parameter_uncertainty_source"], HEURISTIC_PARAMETER_UNCERTAINTY_SOURCE)
        self.assertIn("research_utility_v2", result["metric_summaries"])
        self.assertEqual(result["metric_summaries"]["research_utility_v2"]["status"], "completed")