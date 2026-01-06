from __future__ import annotations

from django.test import TestCase
from django.urls import reverse


class SimulatorAliasRoutingTests(TestCase):
    def test_reverse_canonical_twin_preview_path(self) -> None:
        url = reverse("simulator:twin_preview")
        self.assertTrue(url.startswith("/sim/"), url)
        self.assertEqual(url, "/sim/api/twin/preview/")

    def test_reverse_alias_twin_preview_path(self) -> None:
        url = reverse("simulator_alias:twin_preview")
        self.assertTrue(url.startswith("/simulator/"), url)
        self.assertEqual(url, "/simulator/api/twin/preview/")
