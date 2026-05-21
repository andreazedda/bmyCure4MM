from __future__ import annotations

from django.contrib.auth import get_user_model
from django.test import Client, TestCase
from django.urls import reverse


class AlgorithmTransparencyViewTests(TestCase):
    def setUp(self) -> None:
        self.client = Client()
        User = get_user_model()
        self.user = User.objects.create_user(username="algo-user", password="pass")
        self.client.force_login(self.user)

    def test_algorithm_transparency_uses_exploratory_language(self) -> None:
        response = self.client.get(reverse("simulator:algorithm_transparency"))

        self.assertEqual(response.status_code, 200)

        content = response.content.decode("utf-8").lower()
        self.assertIn("algorithm logic transparency", content)
        self.assertIn("heuristic", content)
        self.assertIn("literature-informed", content)
        self.assertIn("not a prescribing tool", content)
        self.assertIn("what can be concluded", content)
        self.assertIn("what cannot be concluded", content)
        self.assertIn("heuristic", content)
        self.assertIn("literature-informed", content)
        self.assertIn("threshold-based", content)
        self.assertIn("back to simulator scenarios", content)
        self.assertIn("open simple research view", content)
        self.assertIn("open scientific cockpit", content)

        forbidden_phrases = [
            "best treatment",
            "recommended therapy",
            "clinically superior",
            "should have been treated",
            "proven outcome",
            "strong recommendation",
            "moderate recommendation",
            "treatment recommendation",
            "recommendations are made",
            "generates treatment recommendations",
        ]
        for phrase in forbidden_phrases:
            self.assertNotIn(phrase, content)