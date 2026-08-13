from __future__ import annotations

import json
import tempfile
from datetime import date
from pathlib import Path
from unittest.mock import patch

from django.contrib.auth import get_user_model
from django.test import Client, TestCase, override_settings
from django.urls import reverse

from clinic.models import Assessment, Patient
from twin_engine.cockpit import (
    build_assessment_recommendations,
    build_lab_chart_data,
    build_scenario_rows,
    list_latest_completed_runs_by_label,
    summarize_checks,
    write_local_feedback,
)
from twin_engine.developer_checks import detect_schedule_collapse, run_developer_checks
from twin_engine.exposure_bridge import build_exposure_profile
from twin_engine.models import CounterfactualRun, LongitudinalLabResult, SimulationRunMetadata
from twin_engine.privacy import scan_text_for_sensitive_markers
from twin_engine.state_model import initialize_from_assessment


class ResearchCockpitTests(TestCase):
    def setUp(self) -> None:
        User = get_user_model()
        self.owner = User.objects.create_user("cockpit-owner", password="pass1234")
        self.stranger = User.objects.create_user("cockpit-stranger", password="pass1234")
        self.staff = User.objects.create_user("cockpit-staff", password="pass1234", is_staff=True)
        self.client = Client()
        self.patient = Patient.objects.create(
            mrn="MM-COCKPIT-001",
            owner=self.owner,
            first_name="Private",
            last_name="Person",
            birth_date=date(1964, 4, 4),
            sex="F",
            diagnosis_date=date(2024, 1, 1),
        )
        self.assessment_low = Assessment.objects.create(
            patient=self.patient,
            date=date(2024, 1, 5),
            flc_ratio=2.1,
            hemoglobin_g_dl=11.7,
        )
        self.assessment_recommended = Assessment.objects.create(
            patient=self.patient,
            date=date(2024, 5, 20),
            m_protein_g_dl=1.3,
            flc_ratio=2.8,
            hemoglobin_g_dl=11.5,
            creatinine_mg_dl=0.9,
        )
        self.assessment_later_tie = Assessment.objects.create(
            patient=self.patient,
            date=date(2024, 10, 2),
            m_protein_g_dl=1.1,
            flc_ratio=2.4,
            hemoglobin_g_dl=11.9,
            creatinine_mg_dl=1.0,
        )
        LongitudinalLabResult.objects.create(
            patient=self.patient,
            assessment=self.assessment_recommended,
            date=date(2024, 5, 20),
            analyte=LongitudinalLabResult.ANALYTE_M_PROTEIN,
            value=1.25,
            unit="g/dL",
        )
        LongitudinalLabResult.objects.create(
            patient=self.patient,
            date=date(2024, 6, 20),
            analyte=LongitudinalLabResult.ANALYTE_ALT,
            value=42,
            unit="U/L",
        )
        self.state = initialize_from_assessment(self.assessment_recommended, user=self.owner)

    def test_cockpit_route_renders_core_sections_without_direct_identifiers(self) -> None:
        self.client.login(username="cockpit-owner", password="pass1234")
        response = self.client.get(reverse("twin_engine:research_cockpit", args=[self.patient.id]))
        self.assertEqual(response.status_code, 200)
        for label in [
            "Research What-if Cockpit",
            "Data availability",
            "Twin Inputs over time",
            "Calibration quality",
            "What-if scenarios",
            "Toxicity constraints",
            "Causality status",
            "Developer Console",
            "Utility v2",
        ]:
            self.assertContains(response, label)
        self.assertContains(response, f"Research Patient {self.patient.id}")
        self.assertNotContains(response, self.patient.first_name)
        self.assertNotContains(response, self.patient.last_name)

    def test_non_owner_cannot_access_cockpit(self) -> None:
        self.client.login(username="cockpit-stranger", password="pass1234")
        response = self.client.get(reverse("twin_engine:research_cockpit", args=[self.patient.id]))
        self.assertEqual(response.status_code, 403)

    def test_initialize_get_recommends_earliest_highest_completeness_assessment(self) -> None:
        guidance = build_assessment_recommendations(self.patient)
        self.assertEqual(guidance["recommended"]["assessment"].date, date(2024, 5, 20))
        self.client.login(username="cockpit-owner", password="pass1234")
        response = self.client.get(reverse("twin_engine:initialize_twin_from_assessment", args=[self.patient.id]))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Initialize Twin")
        self.assertContains(response, "2024-05-20")

    def test_lab_chart_prefers_structured_lab_values(self) -> None:
        chart_data = build_lab_chart_data(self.patient)
        disease_series = chart_data["disease"]["series"]
        m_protein = next(item for item in disease_series if item["analyte"] == LongitudinalLabResult.ANALYTE_M_PROTEIN)
        self.assertEqual(m_protein["points"][0]["y"], 1.25)
        self.assertEqual(m_protein["points"][0]["source"], "Unknown")

    def test_latest_completed_runs_per_label_are_sorted_by_utility(self) -> None:
        old_run = self._create_run("SAME_LABEL", utility=0.5)
        latest_run = self._create_run("SAME_LABEL", utility=1.5)
        high_run = self._create_run("HIGH_LABEL", utility=2.0)
        latest = list_latest_completed_runs_by_label(self.patient)
        self.assertEqual(latest["SAME_LABEL"], latest_run)
        self.assertNotIn(old_run, latest.values())
        rows = build_scenario_rows(self.patient, latest, [])
        self.assertEqual([row["run"] for row in rows], [high_run, latest_run])

    def test_developer_console_staff_only(self) -> None:
        self.client.login(username="cockpit-owner", password="pass1234")
        response = self.client.get(reverse("twin_engine:developer_console"))
        self.assertEqual(response.status_code, 403)
        self.client.login(username="cockpit-staff", password="pass1234")
        response = self.client.get(reverse("twin_engine:developer_console"))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Developer Console")
        self.assertContains(response, "privacy checks")

    def test_developer_checks_have_actionable_shape(self) -> None:
        summary = summarize_checks(run_developer_checks(self.patient))
        self.assertGreater(summary["total"], 0)
        for checks in summary["groups"].values():
            for check in checks:
                self.assertIn("status", check)
                self.assertIn("detail", check)
                self.assertIn("next_action", check)

    def test_cockpit_renders_validation_uncertainty_and_robustness_section(self) -> None:
        run = self._create_run(
            "VALIDATED_SCENARIO",
            utility=1.2,
            simulation_summary={
                "label": "research simulation",
                "predicted_biomarkers": {"m_protein_g_dl": 1.0},
                "classification": {"counterfactual_class": "mechanistic"},
            },
        )
        uncertainty = {
            "status": "completed",
            "parameter_uncertainty_source": "heuristic_perturbation_not_calibrated_distribution",
            "metric_summaries": {
                "research_utility_v2": {"status": "completed", "p05": 0.9, "median": 1.1, "p95": 1.2, "uncertainty_classification": "moderate"},
                "tumor_reduction": {"p05": 0.3, "median": 0.4, "p95": 0.5},
                "healthy_loss": {"p05": 0.1, "median": 0.2, "p95": 0.3},
                "durability_index": {"p05": 0.2, "median": 0.3, "p95": 0.4},
                "liver_toxicity_signal_0_1": {"p05": 0.1, "median": 0.2, "p95": 0.3},
                "neutropenia_signal_0_1": {"p05": 0.1, "median": 0.2, "p95": 0.3},
            },
        }
        sensitivity = {
            "status": "completed",
            "top_drivers": [{"parameter": "tumor_growth_rate", "max_abs_utility_v2_delta": 0.12, "sensitivity_classification": "moderate"}],
            "unstable_parameters": ["tumor_growth_rate"],
        }
        SimulationRunMetadata.objects.create(
            counterfactual_run=run,
            twin_state=self.state,
            model_version="research-twin-v1",
            solver_name="counterfactual_uncertainty",
            solver_parameters={"diagnostic_summary": uncertainty},
            input_hash="u",
        )
        SimulationRunMetadata.objects.create(
            counterfactual_run=run,
            twin_state=self.state,
            model_version="research-twin-v1",
            solver_name="counterfactual_sensitivity",
            solver_parameters={"diagnostic_summary": sensitivity},
            input_hash="s",
        )
        SimulationRunMetadata.objects.create(
            twin_state=self.state,
            model_version="research-twin-v1",
            solver_name="rolling_origin_backtest",
            solver_parameters={"diagnostic_summary": {"status": "completed", "n_folds": 2, "overall_rmse": 0.2, "overall_mae": 0.15}},
            input_hash="x",
        )
        SimulationRunMetadata.objects.create(
            twin_state=self.state,
            model_version="research-twin-v1",
            solver_name="robust_scenario_ranking",
            solver_parameters={"diagnostic_summary": {"status": "completed", "n_scenarios": 1, "n_aligned_samples": 3, "rows": [{"scenario_label": "VALIDATED_SCENARIO", "point_rank": 1, "robust_rank": 1, "probability_best": 1.0, "p05": 0.9, "p95": 1.2, "robustness_classification": "robust"}], "unstable_rank_flag": False}},
            input_hash="y",
        )

        self.client.login(username="cockpit-owner", password="pass1234")
        response = self.client.get(reverse("twin_engine:research_cockpit", args=[self.patient.id]))
        self.assertEqual(response.status_code, 200)
        self.assertContains(response, "Validation, uncertainty, and robustness")
        self.assertContains(response, "Rolling backtest")
        self.assertContains(response, "VALIDATED_SCENARIO")

    def test_schedule_collapse_classifies_average_exposure_timing_difference(self) -> None:
        daily_profile = build_exposure_profile(
            drug="lenalidomide",
            dose_mg=5,
            horizon_days=10,
            schedule={"type": "daily"},
        ).to_dict()
        alternate_day_profile = build_exposure_profile(
            drug="lenalidomide",
            dose_mg=10,
            horizon_days=10,
            schedule={"type": "interval", "every_days": 2},
        ).to_dict()
        run_a = self._create_run(
            "LEN_5MG_DAILY",
            utility=1.0,
            simulation_summary={
                "label": "research simulation",
                "predicted_biomarkers": {"m_protein_g_dl": 1.0},
                "classification": {"counterfactual_class": "mechanistic"},
                "alternative": {"exposure_profiles": {"lenalidomide": daily_profile}},
            },
        )
        run_b = self._create_run(
            "LEN_10MG_ALT",
            utility=1.0,
            simulation_summary={
                "label": "research simulation",
                "predicted_biomarkers": {"m_protein_g_dl": 1.0},
                "classification": {"counterfactual_class": "mechanistic"},
                "alternative": {"exposure_profiles": {"lenalidomide": alternate_day_profile}},
            },
        )
        shared_trajectory = {
            "alternative_trajectory": {
                "tumor_cells": [1.0, 0.8, 0.6],
                "healthy_cells": [1.0, 1.0, 1.0],
            }
        }
        with patch(
            "twin_engine.developer_checks._load_trajectory",
            side_effect=lambda run: shared_trajectory if run.id in {run_a.id, run_b.id} else None,
        ):
            classifications = detect_schedule_collapse(self.patient)
        self.assertEqual(len(classifications), 1)
        self.assertEqual(classifications[0]["classification"], "AVERAGE_EXPOSURE_COLLAPSE")

    def test_privacy_scanner_flags_direct_identifier_markers(self) -> None:
        sample_name = "Ros" + "sana " + "Agu" + "eci"
        sample_tax_code = "ABC" + "DEF" + "12" + "G34" + "H567" + "I"
        marker = "Codice" + " Fiscale"
        findings = scan_text_for_sensitive_markers(f"{sample_name} {marker} {sample_tax_code}")
        self.assertGreaterEqual(len(findings), 2)
        self.assertTrue(all(item["status"] == "fail" for item in findings))

    def test_local_feedback_writes_under_local_private(self) -> None:
        with tempfile.TemporaryDirectory() as tmpdir:
            with override_settings(BASE_DIR=Path(tmpdir)):
                path = write_local_feedback(self.staff, {"patient_id": self.patient.id, "message": "Review artifact"})
                self.assertIn("local_private", path.parts)
                payload = json.loads(path.read_text(encoding="utf-8").strip())
                self.assertEqual(payload["payload"]["message"], "Review artifact")

    def test_pre_push_script_contains_required_commands(self) -> None:
        script = Path("scripts/pre_push_research_safety_check.sh").read_text(encoding="utf-8")
        self.assertIn("manage.py check", script)
        self.assertIn("manage.py test twin_engine.tests", script)
        self.assertIn("clinic.tests.test_patient_crud", script)

    def _create_run(self, label: str, *, utility: float, simulation_summary: dict | None = None) -> CounterfactualRun:
        return CounterfactualRun.objects.create(
            patient=self.patient,
            base_twin_state=self.state,
            intervention_definition={"label": label},
            simulation_summary=simulation_summary
            or {
                "label": "research simulation",
                "predicted_biomarkers": {"m_protein_g_dl": 1.0},
                "classification": {"counterfactual_class": "mechanistic"},
            },
            comparison_metrics={"research_utility": utility, "toxicity_constraint_penalty": 0.0},
            status=CounterfactualRun.STATUS_COMPLETED,
            created_by=self.owner,
        )
