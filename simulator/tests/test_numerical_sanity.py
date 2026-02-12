from __future__ import annotations

import math
import tempfile
from pathlib import Path

import pandas as pd
from django.conf import settings
from django.contrib.auth import get_user_model
from django.test import TestCase, override_settings

from simulator import models


class SimulationNumericalSanityTests(TestCase):
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("simuser", password="pass1234")
        self.scenario = models.Scenario.objects.create(
            title="Numerical sanity",
            clinical_stage="newly_diagnosed",
            summary="Sanity checks for solver outputs",
            risk_stratification="Standard",
            guideline_notes="",
        )

    def _parameters(self) -> dict:
        return {
            "baseline_tumor_cells": 1.0e9,
            "baseline_healthy_cells": 5.0e11,
            "lenalidomide_dose": 25.0,
            "bortezomib_dose": 1.3,
            "daratumumab_dose": 16.0,
            "time_horizon": 60.0,
            "tumor_growth_rate": 0.023,
            "healthy_growth_rate": 0.015,
            "interaction_strength": 0.05,
            "immune_compromise_index": 1.0,
        }

    def test_outputs_are_finite_and_well_formed(self) -> None:
        with tempfile.TemporaryDirectory(prefix="bmycure4mm-media-") as tmp:
            media_root = Path(tmp)
            with override_settings(MEDIA_ROOT=media_root, MEDIA_URL="/media/"):
                attempt = models.SimulationAttempt.objects.create(
                    scenario=self.scenario,
                    user=self.user,
                    parameters=self._parameters(),
                )
                summary = attempt.run_model()

                # Summary sanity
                for key in ("tumor_reduction", "healthy_loss", "durability_index"):
                    self.assertIn(key, summary)
                    self.assertTrue(math.isfinite(float(summary[key])))

                durability = float(summary["durability_index"])
                self.assertGreaterEqual(durability, 0.0)
                self.assertLessEqual(durability, 1.0)

                tumor_reduction = float(summary["tumor_reduction"])
                healthy_loss = float(summary["healthy_loss"])

                # Allow some growth (negative reduction / negative loss), but avoid absurd values.
                self.assertGreater(tumor_reduction, -5.0)
                self.assertLess(tumor_reduction, 1.5)
                self.assertGreater(healthy_loss, -5.0)
                self.assertLess(healthy_loss, 1.5)

                # Artifacts written
                attempt.refresh_from_db()
                self.assertIn("csv", attempt.results)
                self.assertIn("plot", attempt.results)

                csv_name = Path(attempt.results["csv"]).name
                plot_name = Path(attempt.results["plot"]).name

                csv_path = Path(settings.MEDIA_ROOT) / "sim_data" / csv_name
                plot_path = Path(settings.MEDIA_ROOT) / "sim_plots" / plot_name

                self.assertTrue(csv_path.exists())
                self.assertTrue(plot_path.exists())
                self.assertGreater(csv_path.stat().st_size, 100)
                self.assertGreater(plot_path.stat().st_size, 1000)

                # CSV trajectory sanity
                df = pd.read_csv(csv_path)
                expected_cols = {
                    "time",
                    "tumor_cells",
                    "healthy_cells",
                    "lenalidomide_concentration",
                    "bortezomib_concentration",
                    "daratumumab_concentration",
                }
                self.assertTrue(expected_cols.issubset(set(df.columns)))

                numeric_df = df[list(expected_cols)].astype(float)
                self.assertFalse(numeric_df.isna().any().any())
                self.assertTrue((numeric_df.abs() < 1e30).all().all())
                self.assertTrue((numeric_df[["tumor_cells", "healthy_cells"]] >= 0.0).all().all())

                times = numeric_df["time"].to_list()
                self.assertGreaterEqual(min(times), 0.0)
                self.assertLessEqual(max(times), float(self._parameters()["time_horizon"]) + 1e-6)
                self.assertEqual(times, sorted(times))
