from __future__ import annotations

from django.contrib.auth import get_user_model
from django.test import TestCase
from django.urls import reverse

from simulator import models
from simulator.models_help import HelpArticle


class HelpApiTests(TestCase):
    def setUp(self) -> None:
        self.article = HelpArticle.objects.create(
            slug="test_help",
            title_en="Test Help",
            body_en="<p>English body</p>",
            title_it="Aiuto Test",
            body_it="<p>Corpo italiano</p>",
        )

    def test_help_item_language_and_type(self) -> None:
        url = reverse("api_help_item", kwargs={"slug": self.article.slug})
        resp_en = self.client.get(url, {"lang": "en"})
        self.assertEqual(resp_en.status_code, 200)
        data_en = resp_en.json()
        self.assertEqual(data_en["title"], "Test Help")
        self.assertEqual(data_en["type"], "article")

        resp_it = self.client.get(url, {"lang": "it"})
        self.assertEqual(resp_it.status_code, 200)
        data_it = resp_it.json()
        self.assertEqual(data_it["title"], "Aiuto Test")

    def test_help_search_matches_slug_with_type(self) -> None:
        url = reverse("api_help_search")
        resp = self.client.get(url, {"q": "test", "lang": "en"})
        self.assertEqual(resp.status_code, 200)
        items = resp.json()["results"]
        self.assertIn(
            self.article.slug,
            {item["slug"] for item in items},
        )
        self.assertTrue(any(item["type"] == "article" for item in items))

    def test_help_search_finds_field_and_kpi_types(self) -> None:
        url = reverse("api_help_search")
        resp = self.client.get(url, {"q": "lenalidomide", "lang": "en"})
        self.assertEqual(resp.status_code, 200)
        payload = resp.json()["results"]
        self.assertIn("lenalidomide_dose", {item["slug"] for item in payload})
        self.assertTrue(any(item["type"] == "field" and item["slug"] == "lenalidomide_dose" for item in payload))

        resp_kpi = self.client.get(url, {"q": "kpi_auc", "lang": "en"})
        self.assertEqual(resp_kpi.status_code, 200)
        self.assertTrue(any(item["type"] == "kpi" for item in resp_kpi.json()["results"]))

    def test_preset_search_returns_preset_type(self) -> None:
        url = reverse("api_help_search")
        resp = self.client.get(url, {"q": "VRd", "lang": "en"})
        self.assertEqual(resp.status_code, 200)
        self.assertTrue(any(item["type"] == "preset" for item in resp.json()["results"]))


class HelpTemplateIntegrationTests(TestCase):
    def setUp(self) -> None:
        self.user = get_user_model().objects.create_user("demo", password="demo1234")
        self.scenario = models.Scenario.objects.create(
            title="Scenario",
            clinical_stage="newly_diagnosed",
            summary="Summary",
        )
        params = {
            "baseline_tumor_cells": 1.0e9,
            "baseline_healthy_cells": 5.0e11,
            "lenalidomide_dose": 25.0,
            "bortezomib_dose": 1.3,
            "daratumumab_dose": 16.0,
            "time_horizon": 30.0,
            "tumor_growth_rate": 0.023,
            "healthy_growth_rate": 0.015,
            "interaction_strength": 0.05,
        }
        self.attempt = models.SimulationAttempt.objects.create(
            scenario=self.scenario,
            user=self.user,
            parameters=params,
            seed=123,
        )
        self.attempt.run_model()

    def test_data_help_present_in_results(self) -> None:
        self.client.force_login(self.user)
        resp = self.client.get(reverse("simulator:scenario_detail", args=[self.scenario.pk]))
        self.assertContains(resp, 'data-help="kpi_tumor_reduction"')

    def test_seed_persisted_in_artifacts(self) -> None:
        self.assertEqual(self.attempt.artifacts.get("seed"), 123)
