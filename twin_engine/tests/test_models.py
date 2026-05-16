from __future__ import annotations

from datetime import date

from django.contrib.auth import get_user_model
from django.test import TestCase

from clinic.models import Assessment, Patient, Regimen
from twin_engine.causal import classify_estimand
from twin_engine.models import (
    AdverseEvent,
    CausalAssumptionSet,
    CounterfactualRun,
    LongitudinalLabResult,
    ObservationResidual,
    TherapyInterruption,
)
from twin_engine.provenance import CURRENT_MODEL_VERSION, record_simulation_metadata
from twin_engine.state_model import initialize_from_assessment


class TwinEngineModelTests(TestCase):
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("researcher", password="pass1234")
        self.patient = Patient.objects.create(
            mrn="MM-RES-001",
            owner=self.user,
            first_name="Synthetic",
            last_name="Patient",
            birth_date=date(1968, 1, 1),
            sex="F",
            diagnosis_date=date(2024, 1, 1),
        )
        self.assessment = Assessment.objects.create(
            patient=self.patient,
            date=date(2025, 1, 15),
            m_protein_g_dl=1.2,
            flc_ratio=3.4,
            hemoglobin_g_dl=11.8,
            r_iss="II",
            beta2m_mg_l=3.1,
            ldH_u_l=240,
        )
        self.regimen = Regimen.objects.create(
            name="VRd",
            line="frontline",
            components="Lenalidomide, Bortezomib",
        )

    def test_patient_twin_state_creation(self) -> None:
        state = initialize_from_assessment(self.assessment, user=self.user)
        self.assertEqual(state.patient, self.patient)
        self.assertEqual(state.assessment, self.assessment)
        self.assertTrue(state.is_current)
        self.assertEqual(state.model_version, CURRENT_MODEL_VERSION)

    def test_one_current_state_behavior(self) -> None:
        state_one = initialize_from_assessment(self.assessment, user=self.user)
        second_assessment = Assessment.objects.create(
            patient=self.patient,
            date=date(2025, 2, 15),
            m_protein_g_dl=0.9,
            flc_ratio=2.6,
            hemoglobin_g_dl=12.1,
            r_iss="II",
            beta2m_mg_l=2.9,
            ldH_u_l=220,
        )
        state_two = initialize_from_assessment(second_assessment, user=self.user)
        state_one.refresh_from_db()
        state_two.refresh_from_db()
        self.assertFalse(state_one.is_current)
        self.assertTrue(state_two.is_current)

    def test_observation_residual_metrics_persistence(self) -> None:
        state = initialize_from_assessment(self.assessment, user=self.user)
        residual = ObservationResidual.objects.create(
            patient=self.patient,
            twin_state=state,
            assessment=self.assessment,
            predicted_values={"m_protein_g_dl": 1.0},
            observed_values={"m_protein_g_dl": 1.2},
            residuals={"m_protein_g_dl": 0.2},
            normalized_residuals={"m_protein_g_dl": 0.1667},
            rmse=0.2,
            mae=0.2,
            biomarker_weights={"m_protein_g_dl": 1.0},
            stage=ObservationResidual.STAGE_PRE_CALIBRATION,
        )
        self.assertEqual(residual.patient, self.patient)
        self.assertAlmostEqual(residual.rmse, 0.2)
        self.assertAlmostEqual(residual.mae, 0.2)
        self.assertEqual(residual.stage, ObservationResidual.STAGE_PRE_CALIBRATION)

    def test_longitudinal_lab_result_creation(self) -> None:
        result = LongitudinalLabResult.objects.create(
            patient=self.patient,
            assessment=self.assessment,
            date=self.assessment.date,
            analyte=LongitudinalLabResult.ANALYTE_AST,
            value=25.0,
            unit="U/L",
            source_quality=LongitudinalLabResult.SOURCE_QUALITY_MANUALLY_CURATED,
        )
        self.assertEqual(result.analyte, LongitudinalLabResult.ANALYTE_AST)
        self.assertEqual(result.unit, "U/L")

    def test_therapy_interruption_creation(self) -> None:
        interruption = TherapyInterruption.objects.create(
            patient=self.patient,
            start_date=date(2025, 1, 20),
            end_date=date(2025, 1, 28),
            drug="lenalidomide",
            reason=TherapyInterruption.REASON_NEUTROPENIA,
            evidence={"absolute_neutropenia": True},
        )
        self.assertEqual(interruption.reason, TherapyInterruption.REASON_NEUTROPENIA)

    def test_adverse_event_creation(self) -> None:
        event = AdverseEvent.objects.create(
            patient=self.patient,
            date=date(2025, 5, 7),
            event_type=AdverseEvent.TYPE_HYPERTRANSAMINASEMIA,
            suspected_drug="lenalidomide",
            observed_values={"AST": 406, "ALT": 360},
        )
        self.assertEqual(event.event_type, AdverseEvent.TYPE_HYPERTRANSAMINASEMIA)

    def test_counterfactual_run_persistence(self) -> None:
        state = initialize_from_assessment(self.assessment, user=self.user)
        run = CounterfactualRun.objects.create(
            patient=self.patient,
            base_twin_state=state,
            alternative_regimen=self.regimen,
            intervention_definition={"drug_doses": {"lenalidomide": 25.0}},
            simulation_summary={"label": "research simulation"},
            status=CounterfactualRun.STATUS_COMPLETED,
            created_by=self.user,
        )
        self.assertEqual(run.patient, self.patient)
        self.assertEqual(run.status, CounterfactualRun.STATUS_COMPLETED)
        self.assertEqual(run.simulation_summary["label"], "research simulation")

    def test_simulation_run_metadata_persistence(self) -> None:
        state = initialize_from_assessment(self.assessment, user=self.user)
        metadata = record_simulation_metadata(
            twin_state=state,
            model_version=CURRENT_MODEL_VERSION,
            solver_name="MathematicalModel",
            input_payload={"state_id": state.id},
            solver_parameters={"tumor_growth_rate": 0.023},
            output_payload={"summary": {"tumor_reduction": 0.4}},
        )
        self.assertEqual(metadata.twin_state, state)
        self.assertEqual(metadata.model_version, CURRENT_MODEL_VERSION)
        self.assertTrue(metadata.input_hash)

    def test_causal_assumption_set_classification(self) -> None:
        assumption_set = CausalAssumptionSet.objects.create(
            patient=self.patient,
            name="Mechanistic only",
            graph_definition={"nodes": ["therapy", "response"], "edges": [["therapy", "response"]]},
            variables=["therapy", "response"],
            intervention={"therapy": "VRd"},
            outcome={"response": "tumor_reduction"},
            identification_status=CausalAssumptionSet.IDENT_MECHANISTIC_ONLY,
            created_by=self.user,
        )
        classification = classify_estimand(
            graph_definition=assumption_set.graph_definition,
            intervention=assumption_set.intervention,
            outcome=assumption_set.outcome,
            adjustment_set=assumption_set.adjustment_set,
            identification_status=assumption_set.identification_status,
        )
        self.assertEqual(classification, "mechanistic_simulation_only")
