from __future__ import annotations

from datetime import date
from io import StringIO

from django.core.management import call_command
from django.test import TestCase

from clinic.models import Assessment, Patient
from twin_engine.models import AdverseEvent, CounterfactualRun, LongitudinalLabResult, SimulationRunMetadata, TherapyInterruption
from twin_engine.state_model import initialize_from_assessment



class DiagnosticCommandsTests(TestCase):
    def setUp(self) -> None:
        self.patient = Patient.objects.create(
            mrn="MM-CMD-001",
            first_name="Research",
            last_name="CommandPatient",
            birth_date=date(1965, 5, 5),
            sex="F",
            diagnosis_date=date(2024, 1, 1),
        )
        self.assessments = [
            Assessment.objects.create(patient=self.patient, date=date(2025, 1, 1), m_protein_g_dl=1.2, flc_ratio=2.4, hemoglobin_g_dl=11.8, r_iss="II", beta2m_mg_l=3.0, ldH_u_l=220),
            Assessment.objects.create(patient=self.patient, date=date(2025, 1, 15), m_protein_g_dl=1.1, flc_ratio=2.2, hemoglobin_g_dl=11.9, r_iss="II", beta2m_mg_l=2.9, ldH_u_l=215),
            Assessment.objects.create(patient=self.patient, date=date(2025, 1, 29), m_protein_g_dl=1.0, flc_ratio=2.0, hemoglobin_g_dl=12.0, r_iss="II", beta2m_mg_l=2.8, ldH_u_l=210),
        ]
        self.state = initialize_from_assessment(self.assessments[0])
        CounterfactualRun.objects.create(
            patient=self.patient,
            base_twin_state=self.state,
            intervention_definition={
                "label": "LEN_A",
                "classification": "mechanistic_simulation_only",
                "intervention": {"drug": "lenalidomide", "dose_mg": 5.0, "schedule": {"type": "daily"}, "duration_days": 14, "start_day": 0},
            },
            comparison_metrics={"research_utility_v2": 1.0},
            status=CounterfactualRun.STATUS_COMPLETED,
        )
        CounterfactualRun.objects.create(
            patient=self.patient,
            base_twin_state=self.state,
            intervention_definition={
                "label": "LEN_B",
                "classification": "mechanistic_simulation_only",
                "intervention": {"drug": "lenalidomide", "dose_mg": 10.0, "schedule": {"type": "interval", "every_days": 2}, "duration_days": 14, "start_day": 0},
            },
            comparison_metrics={"research_utility_v2": 0.9},
            status=CounterfactualRun.STATUS_COMPLETED,
        )

    def test_backtest_command_persists_summary_without_phi(self) -> None:
        output = StringIO()
        call_command("run_patient_backtest", patient_id=self.patient.id, stdout=output)
        rendered = output.getvalue()
        self.assertIn('"backtest"', rendered)
        self.assertNotIn(self.patient.first_name, rendered)
        self.assertNotIn(self.patient.last_name, rendered)
        self.assertTrue(SimulationRunMetadata.objects.filter(solver_name="rolling_origin_backtest").exists())

    def test_uncertainty_sensitivity_and_robustness_commands_store_diagnostics(self) -> None:
        uncertainty_output = StringIO()
        call_command("run_patient_uncertainty", patient_id=self.patient.id, horizon_days=14, samples=3, seed=5, stdout=uncertainty_output)
        sensitivity_output = StringIO()
        call_command("run_patient_sensitivity", patient_id=self.patient.id, horizon_days=14, stdout=sensitivity_output)
        robustness_output = StringIO()
        call_command("run_patient_robustness", patient_id=self.patient.id, stdout=robustness_output)

        for rendered in (uncertainty_output.getvalue(), sensitivity_output.getvalue(), robustness_output.getvalue()):
            self.assertNotIn(self.patient.first_name, rendered)
            self.assertNotIn(self.patient.last_name, rendered)

        run = CounterfactualRun.objects.filter(patient=self.patient, intervention_definition__label="LEN_A").first()
        self.assertEqual(run.comparison_metrics, {"research_utility_v2": 1.0})
        uncertainty = run.metadata_records.filter(solver_name="counterfactual_uncertainty").get()
        sensitivity = run.metadata_records.filter(solver_name="counterfactual_sensitivity").get()
        robustness = SimulationRunMetadata.objects.filter(solver_name="robust_scenario_ranking").get()
        self.assertEqual(uncertainty.solver_parameters["diagnostic_summary"]["status"], "completed")
        self.assertEqual(sensitivity.solver_parameters["diagnostic_summary"]["status"], "completed")
        self.assertEqual(robustness.solver_parameters["diagnostic_summary"]["status"], "completed")
        self.assertIsNotNone(uncertainty.manifest_id)
        self.assertIsNotNone(sensitivity.manifest_id)
        self.assertIsNotNone(robustness.manifest_id)
