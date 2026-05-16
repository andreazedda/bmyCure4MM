from __future__ import annotations

import copy
import tempfile
from datetime import date
from pathlib import Path
from unittest.mock import patch

from django.contrib.auth import get_user_model
from django.test import TestCase, override_settings

from clinic.models import Assessment, Patient, PatientTherapy, Regimen
from twin_engine.calibration import CALIBRATION_STATUS_FAILED_UNRELIABLE, calibrate_patient_parameters
from twin_engine.counterfactual import run_counterfactual
from twin_engine.models import ObservationResidual, TherapyInterruption
from twin_engine.observation_model import compute_residuals, predict_biomarkers
from twin_engine.provenance import hash_json
from twin_engine.simulation_bridge import build_solver_inputs_from_twin_state, run_patient_simulation
from twin_engine.state_model import initialize_from_assessment
from twin_engine.therapy_schedule import build_therapy_schedule, convert_patient_therapies_to_drug_doses


class TwinEngineServiceTests(TestCase):
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("svcuser", password="pass1234")
        self.patient = Patient.objects.create(
            mrn="MM-SVC-001",
            owner=self.user,
            first_name="Synthetic",
            last_name="Service",
            birth_date=date(1965, 5, 5),
            sex="M",
            diagnosis_date=date(2024, 2, 1),
            notes="Private note should never reach artifacts",
        )
        self.assessment = Assessment.objects.create(
            patient=self.patient,
            date=date(2025, 1, 1),
            m_protein_g_dl=1.0,
            flc_ratio=2.0,
            hemoglobin_g_dl=12.0,
            r_iss="II",
            beta2m_mg_l=3.0,
            ldH_u_l=230,
        )
        self.regimen = Regimen.objects.create(
            name="VRd",
            line="frontline",
            components="Lenalidomide, Bortezomib",
        )
        self.therapy = PatientTherapy.objects.create(
            patient=self.patient,
            regimen=self.regimen,
            start_date=date(2025, 1, 1),
            end_date=date(2025, 1, 28),
            doses={"lenalidomide": {"dose": 10, "unit": "mg", "schedule": "alternate day"}},
            cycle_length_days=2,
            days_on=[1],
            source_quality=PatientTherapy.SOURCE_QUALITY_CLINICAL_RECORD,
        )
        self.followup_assessment = Assessment.objects.create(
            patient=self.patient,
            date=date(2025, 1, 15),
            m_protein_g_dl=0.8,
            flc_ratio=1.7,
            hemoglobin_g_dl=12.2,
            r_iss="II",
            beta2m_mg_l=2.8,
            ldH_u_l=210,
        )

    def test_initialize_from_assessment_reuses_build_patient_twin(self) -> None:
        with patch(
            "twin_engine.state_model.build_patient_twin",
            return_value={
                "risk_score": 0.5,
                "tumor_growth_rate": 0.03,
                "healthy_growth_rate": 0.015,
                "carrying_capacity_tumor": 1.0e11,
                "carrying_capacity_healthy": 5.0e11,
                "immune_compromise_index": 1.1,
            },
        ) as mock_build:
            state = initialize_from_assessment(self.assessment, user=self.user)
        mock_build.assert_called_once_with(self.assessment)
        self.assertEqual(state.parameters["tumor_growth_rate"], 0.03)
        self.assertEqual(state.parameters["immune_compromise_index"], 1.1)

    def test_observation_model_returns_expected_predicted_values_for_deterministic_inputs(self) -> None:
        predicted = predict_biomarkers(
            {"tumor_cells": 1.0e9, "healthy_cells": 5.0e11},
            self.assessment,
            0,
            {
                "alpha_M": 2.0e-9,
                "beta_M": 0.1,
                "alpha_F": 4.0,
                "gamma_F": 1.0,
                "T_ref": 1.0e9,
                "Hb_baseline": 12.0,
                "H_ref": 5.0e11,
                "eta_H": 1.0,
            },
        )
        self.assertAlmostEqual(predicted["m_protein_g_dl"], 2.1)
        self.assertAlmostEqual(predicted["flc_ratio"], 4.0)
        self.assertAlmostEqual(predicted["hemoglobin_g_dl"], 12.0)

    def test_residuals_computed_correctly(self) -> None:
        payload = compute_residuals(
            predicted={"m_protein_g_dl": 1.0, "flc_ratio": 2.0},
            observed={"m_protein_g_dl": 1.5, "flc_ratio": 1.0},
            weights={"m_protein_g_dl": 1.0, "flc_ratio": 1.0},
        )
        self.assertAlmostEqual(payload["residuals"]["m_protein_g_dl"], 0.5)
        self.assertAlmostEqual(payload["residuals"]["flc_ratio"], -1.0)
        self.assertIsNotNone(payload["rmse"])

    def test_therapy_schedule_builder_refuses_to_silently_infer_unknown_doses(self) -> None:
        schedule = build_therapy_schedule(self.patient, date(2025, 1, 1), date(2025, 1, 28))
        with self.assertRaises(ValueError):
            convert_patient_therapies_to_drug_doses(schedule, strict=True)

    def test_simulation_bridge_builds_numeric_solver_inputs(self) -> None:
        state = initialize_from_assessment(self.assessment, user=self.user)
        solver_inputs = build_solver_inputs_from_twin_state(
            state,
            therapy_schedule=None,
            horizon_days=30,
            overrides={"drug_doses": {"lenalidomide": 25.0}},
        )
        self.assertIsInstance(solver_inputs["baseline_tumor_cells"], float)
        self.assertIsInstance(solver_inputs["baseline_healthy_cells"], float)
        self.assertTrue(all(isinstance(value, float) for value in solver_inputs["drug_doses"].values()))

    def test_run_patient_simulation_summary_emits_predicted_biomarkers(self) -> None:
        state = initialize_from_assessment(self.assessment, user=self.user)
        result = run_patient_simulation(state, therapy_schedule=None, horizon_days=30)
        self.assertIn("predicted_biomarkers", result["summary"])
        self.assertIn("m_protein_g_dl", result["summary"]["predicted_biomarkers"])
        self.assertIn("predicted_biomarkers_milestones", result["summary"])

    def test_calibration_persists_pre_and_post_residuals_and_diagnostics(self) -> None:
        state = initialize_from_assessment(self.assessment, user=self.user)

        def fake_optimizer(**kwargs):
            baseline = kwargs["baseline_diagnostics"]
            improved = copy.deepcopy(baseline)
            improved["objective"] = float(baseline["objective"]) * 0.5
            improved["rmse"] = max(float(baseline["rmse"] or 1.0) * 0.8, 0.0)
            improved["mae"] = max(float(baseline["mae"] or 1.0) * 0.8, 0.0)
            return {
                "optimizer": "least_squares",
                "success": True,
                "message": "synthetic improvement",
                "parameters": copy.deepcopy(state.parameters),
                "diagnostics": improved,
            }

        with patch("twin_engine.calibration._run_optimizer_sequence", side_effect=fake_optimizer):
            result = calibrate_patient_parameters(
                self.patient,
                state,
                [self.assessment, self.followup_assessment],
                self.patient.therapies.all(),
            )

        self.assertGreater(
            ObservationResidual.objects.filter(twin_state=state, stage=ObservationResidual.STAGE_PRE_CALIBRATION).count(),
            0,
        )
        self.assertGreater(
            ObservationResidual.objects.filter(
                twin_state=result["state"],
                stage=ObservationResidual.STAGE_POST_CALIBRATION,
            ).count(),
            0,
        )
        diagnostics = result["state"].parameter_uncertainty
        self.assertIn("objective_before", diagnostics)
        self.assertIn("objective_after", diagnostics)
        self.assertIn("rmse_before", diagnostics)
        self.assertIn("rmse_after", diagnostics)

    def test_failed_optimizer_marks_status_unreliable(self) -> None:
        state = initialize_from_assessment(self.assessment, user=self.user)

        def fake_unreliable_optimizer(**kwargs):
            baseline = kwargs["baseline_diagnostics"]
            unreliable = copy.deepcopy(baseline)
            unreliable["objective"] = max(float(baseline["objective"]) - 1.0, 0.0)
            unreliable["rmse"] = float(baseline["rmse"] or 1.0) + 0.2
            unreliable["mae"] = float(baseline["mae"] or 1.0) + 0.2
            unreliable["n_observations"] = 1
            return {
                "optimizer": "least_squares",
                "success": False,
                "message": "synthetic unreliable optimizer",
                "parameters": copy.deepcopy(state.parameters),
                "diagnostics": unreliable,
            }

        with patch("twin_engine.calibration._run_optimizer_sequence", side_effect=fake_unreliable_optimizer):
            result = calibrate_patient_parameters(
                self.patient,
                state,
                [self.assessment, self.followup_assessment],
                self.patient.therapies.all(),
            )

        self.assertEqual(result["state"].parameter_uncertainty["calibration_status"], CALIBRATION_STATUS_FAILED_UNRELIABLE)

    def test_counterfactual_rich_intervention_schema_emits_utility_and_sanitized_artifact(self) -> None:
        state = initialize_from_assessment(self.assessment, user=self.user)
        TherapyInterruption.objects.create(
            patient=self.patient,
            patient_therapy=self.therapy,
            start_date=date(2025, 1, 10),
            end_date=date(2025, 1, 12),
            drug="lenalidomide",
            reason=TherapyInterruption.REASON_HYPERTRANSAMINASEMIA,
        )
        intervention = {
            "label": "LEN_TEST_DAILY_30D",
            "classification": "mechanistic_simulation_only",
            "intervention": {
                "drug": "lenalidomide",
                "dose_mg": 10.0,
                "schedule": {"type": "daily"},
                "duration_days": 30,
                "start_day": 0,
            },
            "comparison": {"baseline": "current_recorded_therapy"},
            "outcomes": ["tumor_reduction", "healthy_loss", "durability_index"],
            "random_seed": 42,
            "notes": "Research simulation only.",
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            with override_settings(MEDIA_ROOT=tmpdir, MEDIA_URL="/media/"):
                run = run_counterfactual(self.patient, state, intervention, 30, user=self.user)
                report_path = Path(tmpdir) / run.report_artifact.replace("/media/", "")
                content = report_path.read_text(encoding="utf-8")

        self.assertIn("research_utility", run.comparison_metrics)
        self.assertIn("toxicity_constraint_penalty", run.comparison_metrics)
        self.assertIn("predicted_biomarkers", run.simulation_summary)
        self.assertNotIn('"notes"', content)
        self.assertNotIn(self.patient.notes, content)

    def test_counterfactual_run_persists_status_and_summary(self) -> None:
        state = initialize_from_assessment(self.assessment, user=self.user)
        intervention = {
            "drug_doses": {"lenalidomide": 25.0, "bortezomib": 1.3},
            "start_day": 1,
            "duration_days": 14,
            "schedule": {},
            "parameter_overrides": {},
            "random_seed": 7,
        }
        with tempfile.TemporaryDirectory() as tmpdir:
            with override_settings(MEDIA_ROOT=tmpdir, MEDIA_URL="/media/"):
                run = run_counterfactual(self.patient, state, intervention, 14, user=self.user)
                self.assertEqual(run.status, run.STATUS_COMPLETED)
                self.assertEqual(run.simulation_summary["label"], "research simulation")
                self.assertIn("predicted_biomarkers", run.simulation_summary)
                self.assertTrue(run.report_artifact)
                self.assertTrue(Path(tmpdir, "research_reports").exists())

    def test_provenance_hashes_are_stable(self) -> None:
        payload = {"a": 1, "b": [2, 3]}
        self.assertEqual(hash_json(payload), hash_json(payload))
