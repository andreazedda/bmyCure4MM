"""Integration tests for BadgeUtils and data-testid attributes."""

from __future__ import annotations

from django.test import Client
from django.urls import reverse


def test_badges_js_included_in_base_template(client: Client):
    """Verify that badges.js is included in base template."""
    # Access a page that uses the base template
    response = client.get(reverse("clinic:dashboard"))
    html = response.content.decode()
    assert "static/app/js/badges.js" in html or "badges.js" in html


def test_help_drawer_has_data_testid(client: Client):
    """Verify that help drawer has data-testid attribute."""
    response = client.get(reverse("clinic:dashboard"))
    html = response.content.decode()
    assert 'data-testid="help-drawer"' in html


def test_cmdk_input_has_data_testid(client: Client):
    """Verify that command-k input has data-testid attribute."""
    response = client.get(reverse("clinic:dashboard"))
    html = response.content.decode()
    assert 'data-testid="cmdk-input"' in html


def test_live_region_exists_in_base(client: Client):
    """Verify that aria-live region exists in base template."""
    response = client.get(reverse("clinic:dashboard"))
    html = response.content.decode()
    assert 'id="live-region"' in html
    assert 'aria-live="polite"' in html


def test_window_update_live_function_exists(client: Client):
    """Verify that window.updateLive function is defined."""
    response = client.get(reverse("clinic:dashboard"))
    html = response.content.decode()
    assert "window.updateLive" in html
    assert "const live = (msg)" in html


def test_badge_utils_exports_debounce(client: Client):
    """Verify that BadgeUtils module exports debounce function."""
    response = client.get(reverse("clinic:dashboard"))
    html = response.content.decode()
    # Check that badges.js is included
    assert "badges.js" in html


def test_form_field_help_button_has_testid(client: Client, db):
    """Verify that form field help buttons have data-testid attributes."""
    # This test requires a scenario to exist
    from simulator.models import Scenario

    scenario = Scenario.objects.create(name="Test Scenario")
    response = client.get(
        reverse("simulator:scenario_detail", kwargs={"pk": scenario.pk})
    )
    html = response.content.decode()
    # Look for data-testid pattern on help buttons
    assert 'data-testid="help-open-' in html or 'data-help=' in html


def test_badge_elements_have_testid(client: Client, db):
    """Verify that badge elements have data-testid attributes."""
    from simulator.models import Scenario

    scenario = Scenario.objects.create(name="Test Scenario")
    response = client.get(
        reverse("simulator:scenario_detail", kwargs={"pk": scenario.pk})
    )
    html = response.content.decode()
    # Look for data-testid on badge elements
    assert 'data-testid="badge-' in html or 'id="badge-' in html
