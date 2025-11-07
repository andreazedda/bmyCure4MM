from __future__ import annotations

import datetime
import json
from pathlib import Path

from django.conf import settings
from django.contrib.auth import get_user_model
from django.test import TestCase, override_settings

from clinic.models import Assessment, Patient
from simulator import models
from simulator.twin import build_patient_twin


class PatientTwinTests(TestCase):
    def setUp(self) -> None:
        self.patient = Patient.objects.create(
            mrn="TST-001",
            first_name="Test",
            last_name="Twin",
            birth_date=datetime.date(1970, 1, 1),
            sex="M",
            diagnosis_date=datetime.date(2020, 1, 1),
        )
        self.assessment_low = Assessment.objects.create(
            patient=self.patient,
            date=datetime.date(2024, 1, 1),
            r_iss="I",
            ldH_u_l=220,
            beta2m_mg_l=3.0,
            flc_ratio=0.9,
        )
        self.assessment_high = Assessment.objects.create(
            patient=self.patient,
            date=datetime.date(2024, 2, 1),
            r_iss="III",
            ldH_u_l=620,
            beta2m_mg_l=12.0,
            flc_ratio=12.0,
        )
        self.user = get_user_model().objects.create_user("predlab", password="securepass")
        self.scenario = models.Scenario.objects.create(
            title="Twin validation",
            clinical_stage="relapsed_refractory",
            summary="Testing twin impact",
            risk_stratification="Lab-based",
        )

    def test_risk_score_monotonic_with_stage(self) -> None:
        low = build_patient_twin(self.assessment_low)
        high = build_patient_twin(self.assessment_high)
        self.assertLess(low["tumor_growth_rate"], high["tumor_growth_rate"])
        self.assertGreater(low["carrying_capacity_healthy"], high["carrying_capacity_healthy"])

    @override_settings(PREDLAB_V2=True)
    def test_twin_serialization_written_to_results(self) -> None:
        attempt = models.SimulationAttempt.objects.create(
            scenario=self.scenario,
            user=self.user,
            parameters={
                "baseline_tumor_cells": 1.0e9,
                "baseline_healthy_cells": 5.0e11,
                "lenalidomide_dose": 25.0,
                "bortezomib_dose": 1.3,
                "daratumumab_dose": 16.0,
                "time_horizon": 30.0,
                "interaction_strength": 0.05,
                "twin_assessment_id": self.assessment_high.pk,
                "use_twin": True,
            },
        )
        attempt.run_model()
        attempt.refresh_from_db()
        self.assertIn("twin_params.json", attempt.results)
        twin_url = attempt.results["twin_params.json"]
        filename = Path(twin_url).name
        twin_path = Path(settings.MEDIA_ROOT) / "sim_data" / filename
        self.assertTrue(twin_path.exists())
        with twin_path.open("r", encoding="utf-8") as handle:
            payload = json.load(handle)
        self.assertIn("risk_score", payload)
