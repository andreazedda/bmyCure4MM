from __future__ import annotations

from datetime import date, timedelta

from django.contrib.auth import get_user_model
from django.test import TestCase

from clinic.models import Assessment, Patient
from twin_engine.exposure_bridge import build_exposure_profile
from twin_engine.models import AdverseEvent, LongitudinalLabResult
from twin_engine.toxicity_dynamics import compute_toxicity_dynamics


class ToxicityDynamicsTests(TestCase):
    def setUp(self) -> None:
        user = get_user_model().objects.create_user("toxicity-user", password="pass1234")
        self.patient = Patient.objects.create(
            mrn="MM-TOX-001",
            owner=user,
            first_name="Synthetic",
            last_name="Toxicity",
            birth_date=date(1964, 1, 1),
            sex="F",
            diagnosis_date=date(2024, 1, 1),
        )
        self.assessment = Assessment.objects.create(
            patient=self.patient,
            date=date(2024, 1, 1),
            m_protein_g_dl=1.2,
            flc_ratio=2.1,
            hemoglobin_g_dl=11.6,
        )

    def test_liver_signal_increases_with_exposure(self) -> None:
        low = build_exposure_profile(
            drug="lenalidomide",
            dose_mg=5,
            horizon_days=10,
            schedule={"type": "daily"},
        ).to_dict()
        high = build_exposure_profile(
            drug="lenalidomide",
            dose_mg=10,
            horizon_days=10,
            schedule={"type": "daily"},
        ).to_dict()
        low_signal = compute_toxicity_dynamics(self.patient, low)
        high_signal = compute_toxicity_dynamics(self.patient, high)
        self.assertGreater(high_signal["liver_toxicity_signal_0_1"], low_signal["liver_toxicity_signal_0_1"])

    def test_steatosis_increases_liver_signal(self) -> None:
        profile = build_exposure_profile(
            drug="lenalidomide",
            dose_mg=8,
            horizon_days=10,
            schedule={"type": "daily"},
        ).to_dict()
        baseline = compute_toxicity_dynamics(self.patient, profile)
        AdverseEvent.objects.create(
            patient=self.patient,
            date=date(2024, 1, 5),
            event_type=AdverseEvent.TYPE_HEPATIC_STEATOSIS,
        )
        steatosis = compute_toxicity_dynamics(self.patient, profile)
        self.assertGreater(steatosis["liver_toxicity_signal_0_1"], baseline["liver_toxicity_signal_0_1"])

    def test_prior_peak_increases_liver_signal(self) -> None:
        profile = build_exposure_profile(
            drug="lenalidomide",
            dose_mg=8,
            horizon_days=10,
            schedule={"type": "daily"},
        ).to_dict()
        baseline = compute_toxicity_dynamics(self.patient, profile)
        for offset, value in enumerate([20.0, 80.0, 140.0]):
            LongitudinalLabResult.objects.create(
                patient=self.patient,
                assessment=self.assessment,
                date=date(2024, 1, 1) + timedelta(days=offset),
                analyte=LongitudinalLabResult.ANALYTE_AST,
                value=value,
                unit="U/L",
            )
        peaked = compute_toxicity_dynamics(self.patient, profile)
        self.assertGreater(peaked["predicted_ast_signal"], baseline["predicted_ast_signal"])

    def test_neutropenia_signal_increases_with_cumulative_exposure(self) -> None:
        short = build_exposure_profile(
            drug="lenalidomide",
            dose_mg=5,
            horizon_days=5,
            schedule={"type": "daily"},
        ).to_dict()
        long = build_exposure_profile(
            drug="lenalidomide",
            dose_mg=5,
            horizon_days=20,
            schedule={"type": "daily"},
        ).to_dict()
        short_signal = compute_toxicity_dynamics(self.patient, short)
        long_signal = compute_toxicity_dynamics(self.patient, long)
        self.assertGreater(long_signal["neutropenia_signal_0_1"], short_signal["neutropenia_signal_0_1"])

    def test_insufficient_data_marks_default(self) -> None:
        profile = build_exposure_profile(
            drug="lenalidomide",
            dose_mg=5,
            horizon_days=10,
            schedule={"type": "daily"},
        ).to_dict()
        result = compute_toxicity_dynamics(self.patient, profile)
        self.assertEqual(result["coefficient_diagnostics"]["fitted_or_default"], "insufficient_data_default")
        self.assertEqual(result["coefficient_diagnostics"]["parameter_status"]["ast"], "insufficient_data_default")

    def test_enough_data_marks_fitted_or_estimable(self) -> None:
        for analyte, values in [
            (LongitudinalLabResult.ANALYTE_AST, [20.0, 40.0, 60.0]),
            (LongitudinalLabResult.ANALYTE_ALT, [25.0, 35.0, 55.0]),
            (LongitudinalLabResult.ANALYTE_NEU, [2.2, 1.6, 1.1]),
        ]:
            for offset, value in enumerate(values):
                LongitudinalLabResult.objects.create(
                    patient=self.patient,
                    assessment=self.assessment,
                    date=date(2024, 1, 1) + timedelta(days=offset),
                    analyte=analyte,
                    value=value,
                    unit="U/L" if analyte in {LongitudinalLabResult.ANALYTE_AST, LongitudinalLabResult.ANALYTE_ALT} else "10^9/L",
                )
        profile = build_exposure_profile(
            drug="lenalidomide",
            dose_mg=10,
            horizon_days=10,
            schedule={"type": "daily"},
        ).to_dict()
        result = compute_toxicity_dynamics(self.patient, profile)
        self.assertEqual(result["coefficient_diagnostics"]["fitted_or_default"], "fitted_estimate")
        self.assertEqual(result["coefficient_diagnostics"]["parameter_status"]["ast"], "fitted_estimate")
        self.assertEqual(result["coefficient_diagnostics"]["parameter_status"]["neu"], "fitted_estimate")

    def test_no_fake_ast_alt_prediction_if_units_unjustified(self) -> None:
        profile = build_exposure_profile(
            drug="lenalidomide",
            dose_mg=10,
            horizon_days=10,
            schedule={"type": "daily"},
        ).to_dict()
        result = compute_toxicity_dynamics(self.patient, profile)
        self.assertEqual(result["signal_units"], "normalized_risk_signal_0_1")
        self.assertIn("prototype risk signals", result["limitation"])
