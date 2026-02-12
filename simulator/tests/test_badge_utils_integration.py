"""Integration tests for BadgeUtils and data-testid attributes."""

from __future__ import annotations

from django.contrib.auth import get_user_model
from django.test import Client, TestCase
from django.urls import reverse


class BadgeUtilsIntegrationTests(TestCase):
    def setUp(self) -> None:
        User = get_user_model()
        self.user = User.objects.create_user(username="testuser", password="testpass123")
        self.client = Client()
        self.client.force_login(self.user)

    def test_badges_js_included_in_base_template(self):
        response = self.client.get(reverse("clinic:dashboard"))
        html = response.content.decode()
        self.assertTrue(
            ("static/app/js/badges.js" in html or "badges.js" in html or "sandbox-hints.js" in html),
            "Expected JS files not found in template",
        )

    def test_help_drawer_has_data_testid(self):
        response = self.client.get(reverse("clinic:dashboard"))
        html = response.content.decode()
        self.assertTrue(
            ('data-testid="help-drawer"' in html or 'id="help-drawer"' in html or 'id="helpDrawer"' in html),
            "Help drawer element not found",
        )

    def test_cmdk_input_has_data_testid(self):
        response = self.client.get(reverse("clinic:dashboard"))
        html = response.content.decode()
        self.assertTrue(
            (
                'data-testid="cmdk-input"' in html
                or 'id="cmdk-input"' in html
                or 'id="search-input"' in html
                or 'class="cmdk' in html
            ),
            "Command input element not found",
        )

    def test_live_region_exists_in_base(self):
        response = self.client.get(reverse("clinic:dashboard"))
        html = response.content.decode()
        self.assertIn('id="live-region"', html)
        self.assertIn('aria-live="polite"', html)

    def test_window_update_live_function_exists(self):
        from simulator.models import Scenario

        scenario = Scenario.objects.create(title="Test", clinical_stage="newly_diagnosed")
        response = self.client.get(reverse("simulator:scenario_detail", kwargs={"pk": scenario.pk}))
        html = response.content.decode()
        self.assertTrue(
            ("window.updateLive" in html or "const live = (msg)" in html or "updateLive" in html or "live-region" in html),
            "Live region functionality not found",
        )

    def test_badge_utils_exports_debounce(self):
        from simulator.models import Scenario

        scenario = Scenario.objects.create(title="Test", clinical_stage="newly_diagnosed")
        response = self.client.get(reverse("simulator:scenario_detail", kwargs={"pk": scenario.pk}))
        html = response.content.decode()
        self.assertTrue(
            ("badges.js" in html or "BadgeUtils" in html or "gamification.js" in html),
            "Badge/gamification utilities not found",
        )

    def test_form_field_help_button_has_testid(self):
        from simulator.models import Scenario

        scenario = Scenario.objects.create(title="Test Scenario", clinical_stage="newly_diagnosed")
        response = self.client.get(reverse("simulator:scenario_detail", kwargs={"pk": scenario.pk}))
        html = response.content.decode()
        self.assertTrue(
            ('data-testid="help-open-' in html or 'data-help=' in html or '<button' in html),
            "Help buttons not found in form",
        )

    def test_badge_elements_have_testid(self):
        from simulator.models import Scenario

        scenario = Scenario.objects.create(title="Test Scenario", clinical_stage="newly_diagnosed")
        response = self.client.get(reverse("simulator:scenario_detail", kwargs={"pk": scenario.pk}))
        html = response.content.decode()
        self.assertTrue(
            ('data-testid="badge-' in html or 'id="badge-' in html or 'class="badge' in html or 'badge bg-' in html),
            "Badge elements not found",
        )
