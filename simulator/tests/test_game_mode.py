from __future__ import annotations

import tempfile

from django.contrib.auth import get_user_model
from django.test import Client, TestCase, override_settings
from django.urls import reverse

from simulator.models import Scenario


class GameModeTests(TestCase):
    def setUp(self) -> None:
        self.client = Client()
        User = get_user_model()
        self.user = User.objects.create_user(username="player", password="pass")
        self.scenario = Scenario.objects.create(
            title="Game scenario",
            summary="A minimal scenario for game mode testing.",
            active=True,
        )

    def test_simulate_scenario_renders_game_panel_when_enabled(self) -> None:
        self.client.login(username="player", password="pass")

        with tempfile.TemporaryDirectory() as tmpdir:
            with override_settings(MEDIA_ROOT=tmpdir, MEDIA_URL="/media/"):
                url = reverse("simulator:scenario_simulate", args=[self.scenario.pk])
                payload = {
                    "baseline_tumor_cells": 1.0e9,
                    "baseline_healthy_cells": 5.0e11,
                    "lenalidomide_dose": 25.0,
                    "bortezomib_dose": 1.3,
                    "daratumumab_dose": 16.0,
                    "time_horizon": 30,
                    "tumor_growth_rate": 0.023,
                    "healthy_growth_rate": 0.015,
                    "interaction_strength": 0.05,
                    "preset": "VRd",
                    "creatinine_clearance": 90.0,
                    "neuropathy_grade": 0,
                    "anc": 2.0,
                    "platelets": 150,
                    "cohort_size": 1,
                    "use_twin": "",
                    "twin_biology_mode": "auto",
                    "game_mode": "1",
                }
                response = self.client.post(url, data=payload)
                self.assertEqual(response.status_code, 200)
                self.assertIn("Game Mode", response.content.decode("utf-8"))
