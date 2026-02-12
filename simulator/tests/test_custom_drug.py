from __future__ import annotations

from django.contrib.auth import get_user_model
from django.test import TestCase

from simulator import models


class CustomDrugSimulationTests(TestCase):
    def setUp(self) -> None:
        User = get_user_model()
        self.user = User.objects.create_user("custom", password="pass1234", is_staff=True)
        self.scenario = models.Scenario.objects.create(
            title="Custom drug scenario",
            clinical_stage="newly_diagnosed",
            summary="Scenario for custom drug coverage",
            risk_stratification="",
        )

    def test_custom_drug_adds_auc_entry(self) -> None:
        params = {
            "baseline_tumor_cells": 1.0e9,
            "baseline_healthy_cells": 5.0e11,
            "lenalidomide_dose": 25.0,
            "bortezomib_dose": 1.3,
            "daratumumab_dose": 16.0,
            "carfilzomib_dose": 0.0,
            "time_horizon": 60.0,
            "tumor_growth_rate": 0.023,
            "healthy_growth_rate": 0.015,
            "interaction_strength": 0.05,
            "custom_drug_enabled": True,
            "custom_drug_name": "NeverGivenX",
            "custom_drug_dose": 100.0,
            "custom_pk_half_life": 48.0,
            "custom_pk_vd": 40.0,
            "custom_pd_emax": 0.6,
            "custom_pd_ec50": 2.0,
            "custom_drug_key": "custom_nevergivenx",
        }
        attempt = models.SimulationAttempt.objects.create(
            scenario=self.scenario,
            user=self.user,
            parameters=params,
        )
        summary = attempt.run_model()

        auc = summary.get("auc") or {}
        self.assertIn(
            "custom_nevergivenx",
            auc,
            f"Expected custom drug AUC key present. Got: {list(auc.keys())}",
        )
        self.assertGreater(float(auc["custom_nevergivenx"]), 0.0)
