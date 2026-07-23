from __future__ import annotations

from django.contrib.auth import get_user_model
from django.test import TestCase
from django.urls import reverse

from simulator import models
from simulator.views import (
    _build_scenario_input_classification_rows,
    _build_scenario_page_purpose_rows,
)


class ScenarioDetailWorkspaceTests(TestCase):
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("workspace", password=None)
        self.scenario = models.Scenario.objects.create(
            title="Workspace scenario",
            clinical_stage="newly_diagnosed",
            summary="Scenario summary for exploratory workspace testing.",
            active=True,
        )

    def _payload(self) -> dict:
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
            "preset": "VRd",
            "creatinine_clearance": 90.0,
            "neuropathy_grade": 0,
            "anc": 2.0,
            "platelets": 150,
            "cohort_size": 1,
            "use_twin": "",
            "twin_biology_mode": "auto",
        }

    def test_scenario_detail_uses_exploratory_workspace_language(self) -> None:
        self.client.force_login(self.user)
        response = self.client.get(reverse("simulator:scenario_detail", args=[self.scenario.pk]))

        self.assertEqual(response.status_code, 200)
        content = response.content.decode("utf-8").lower()

        for phrase in [
            "exploratory simulation workspace",
            "not a prescribing tool",
            "mechanistic simulations",
            "what can be concluded",
            "what cannot be concluded",
            "input classification",
            "simulated outputs",
            "next actions",
        ]:
            self.assertIn(phrase, content)

        for phrase in [
            "best treatment",
            "recommended therapy",
            "clinically superior",
            "should have been treated",
            "proven outcome",
            "strong recommendation",
            "moderate recommendation",
            "treatment recommendation",
            "good balance",
            "effective?",
            "action plan",
        ]:
            self.assertNotIn(phrase, content)

    def test_simulation_partial_uses_safe_output_labels(self) -> None:
        self.client.force_login(self.user)
        response = self.client.post(
            reverse("simulator:scenario_simulate", args=[self.scenario.pk]),
            self._payload(),
            HTTP_HX_REQUEST="true",
        )

        self.assertEqual(response.status_code, 200)
        content = response.content.decode("utf-8").lower()
        compact = content.replace("\n", "")

        for phrase in [
            "simulated outputs",
            "simulated efficacy signal",
            "simulated healthy-cell impact",
            "exploratory trade-off",
        ]:
            self.assertIn(phrase, content)

        for phrase in [
            "good balance",
            "effective?",
            "action plan",
        ]:
            self.assertNotIn(phrase, content)

        self.assertNotIn(">toxic<", compact)


    def test_input_classification_separates_virtual_and_clinic_sources(self) -> None:
        inactive_rows = _build_scenario_input_classification_rows(
            twin_context_active=False
        )
        inactive_by_label = {
            row["input_label"]: row
            for row in inactive_rows
        }

        self.assertEqual(
            inactive_by_label["Scenario laboratory snapshot"]["classification"],
            "SCENARIO-DEFINED",
        )
        self.assertEqual(
            inactive_by_label["Clinic-linked Patient Twin assessment"]["classification"],
            "UNKNOWN",
        )

        active_rows = _build_scenario_input_classification_rows(
            twin_context_active=True
        )
        active_by_label = {
            row["input_label"]: row
            for row in active_rows
        }

        self.assertEqual(
            active_by_label["Clinic-linked Patient Twin assessment"]["classification"],
            "CLINIC-LINKED STRUCTURED",
        )

    def test_page_purpose_summary_does_not_repeat_direct_identifiers(self) -> None:
        rows = _build_scenario_page_purpose_rows(
            scenario=self.scenario,
            twin_context_active=True,
            latest_simulation=None,
            simulation_runs_count=0,
            decision_logs_count=0,
        )
        rows_by_label = {
            row["label"]: row["value"]
            for row in rows
        }

        self.assertEqual(
            rows_by_label["Clinic-linked Patient Twin"],
            "Active accessible clinic-linked assessment",
        )
        self.assertNotIn("MRN", rows_by_label["Clinic-linked Patient Twin"])
