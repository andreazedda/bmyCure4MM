from __future__ import annotations

from django.test import TestCase

from simulator import forms


class SimulationParameterValidationTests(TestCase):
    def _base_payload(self) -> dict[str, float]:
        return {
            "baseline_tumor_cells": 1.0e9,
            "baseline_healthy_cells": 5.0e11,
            "lenalidomide_dose": 25.0,
            "bortezomib_dose": 1.3,
            "daratumumab_dose": 16.0,
            "time_horizon": 180,
            "tumor_growth_rate": 0.023,
            "healthy_growth_rate": 0.015,
            "interaction_strength": 0.05,
        }

    def test_form_rejects_extreme_values(self) -> None:
        payload = self._base_payload()
        payload["lenalidomide_dose"] = 60.0
        form = forms.SimulationParameterForm(data=payload)
        self.assertFalse(form.is_valid())
        self.assertIn("lenalidomide_dose", form.errors)

    def test_combined_high_doses_raise_error(self) -> None:
        payload = self._base_payload()
        payload["lenalidomide_dose"] = 45.0
        payload["bortezomib_dose"] = 1.6
        form = forms.SimulationParameterForm(data=payload)
        self.assertFalse(form.is_valid())
        self.assertIn("Combined high doses", form.non_field_errors()[0])

    def test_warning_emitted_for_borderline_values(self) -> None:
        payload = self._base_payload()
        payload["lenalidomide_dose"] = 38.0
        form = forms.SimulationParameterForm(data=payload)
        self.assertTrue(form.is_valid())
        self.assertTrue(form.warnings)
        self.assertTrue(any("Lenalidomide dose" in warning for warning in form.warnings))
