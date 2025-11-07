from __future__ import annotations

from django.contrib.auth import get_user_model
from django.test import TestCase, override_settings

from simulator import models


class PKPDRegistryTests(TestCase):
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("pkpd", password="demo1234")
        self.scenario = models.Scenario.objects.create(
            title="PK/PD validation",
            clinical_stage="newly_diagnosed",
            summary="Test scenario",
            risk_stratification="",
        )

    def _parameters(self) -> dict:
        return {
            "baseline_tumor_cells": 1.1e9,
            "baseline_healthy_cells": 5.2e11,
            "lenalidomide_dose": 25.0,
            "bortezomib_dose": 1.3,
            "daratumumab_dose": 16.0,
            "time_horizon": 84.0,
            "tumor_growth_rate": 0.023,
            "healthy_growth_rate": 0.015,
            "interaction_strength": 0.08,
        }

    @override_settings(PREDLAB_V2=True)
    def test_yaml_registry_alters_auc(self) -> None:
        attempt_flag = models.SimulationAttempt.objects.create(
            scenario=self.scenario,
            user=self.user,
            parameters=self._parameters(),
        )
        summary_flag = attempt_flag.run_model()

        with self.settings(PREDLAB_V2=False):
            attempt_legacy = models.SimulationAttempt.objects.create(
                scenario=self.scenario,
                user=self.user,
                parameters=self._parameters(),
            )
            summary_legacy = attempt_legacy.run_model()

        self.assertNotAlmostEqual(
            summary_flag["auc"]["lenalidomide"],
            summary_legacy["auc"]["lenalidomide"],
        )
