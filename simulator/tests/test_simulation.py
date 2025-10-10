from __future__ import annotations

from django.contrib.auth import get_user_model
from django.test import TestCase
from django.urls import reverse

from simulator import models


class SimulationModelTests(TestCase):
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("simuser", password="pass1234")
        self.scenario = models.Scenario.objects.create(
            title="Benchmark scenario",
            clinical_stage="newly_diagnosed",
            summary="Baseline scenario for testing",
            risk_stratification="Standard",
            guideline_notes="",
        )

    def _parameters(self, lenalidomide: float) -> dict:
        return {
            "baseline_tumor_cells": 1.0e9,
            "baseline_healthy_cells": 5.0e11,
            "lenalidomide_dose": lenalidomide,
            "bortezomib_dose": 1.3,
            "daratumumab_dose": 16.0,
            "time_horizon": 60.0,
            "tumor_growth_rate": 0.023,
            "healthy_growth_rate": 0.015,
            "interaction_strength": 0.05,
        }

    def test_higher_dose_reduces_tumor_more(self) -> None:
        low_attempt = models.SimulationAttempt.objects.create(
            scenario=self.scenario,
            user=self.user,
            parameters=self._parameters(lenalidomide=10.0),
        )
        high_attempt = models.SimulationAttempt.objects.create(
            scenario=self.scenario,
            user=self.user,
            parameters=self._parameters(lenalidomide=40.0),
        )
        low_summary = low_attempt.run_model()
        high_summary = high_attempt.run_model()
        self.assertGreater(high_summary["tumor_reduction"], low_summary["tumor_reduction"])

    def test_results_persist_to_database(self) -> None:
        attempt = models.SimulationAttempt.objects.create(
            scenario=self.scenario,
            user=self.user,
            parameters=self._parameters(lenalidomide=25.0),
        )
        attempt.run_model()
        attempt.refresh_from_db()
        self.assertIn("tumor_reduction", attempt.results_summary)
        self.assertIn("csv", attempt.results)

    def test_simulation_view_returns_partial(self) -> None:
        self.client.force_login(self.user)
        url = reverse("simulator:scenario_simulate", args=[self.scenario.pk])
        payload = self._parameters(lenalidomide=25.0)
        response = self.client.post(url, payload, HTTP_HX_REQUEST="true")
        self.assertEqual(response.status_code, 200)
        self.assertIn("Simulation Results", response.content.decode())
