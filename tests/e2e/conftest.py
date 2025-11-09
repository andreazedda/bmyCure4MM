"""
Pytest configuration for E2E tests with Playwright.
"""
import pytest
from django.contrib.auth import get_user_model
from simulator.models import Scenario

User = get_user_model()


@pytest.fixture(scope="function")
def test_user(db):
    """Create a test user for authentication."""
    return User.objects.create_user(
        username="testuser",
        password="testpass123",
        email="test@example.com"
    )


@pytest.fixture(scope="function")
def test_scenario(db):
    """Create a test scenario for simulation tests."""
    return Scenario.objects.create(
        title="E2E Test Scenario",
        clinical_stage="newly_diagnosed",
        description="Test scenario for E2E testing"
    )


@pytest.fixture(scope="function")
def authenticated_page(page, live_server, test_user):
    """Return an authenticated Playwright page."""
    # Navigate to login page
    page.goto(f"{live_server.url}/accounts/login/")
    
    # Fill login form
    page.fill('input[name="username"]', test_user.username)
    page.fill('input[name="password"]', "testpass123")
    page.click('button[type="submit"]')
    
    # Wait for navigation after login
    page.wait_for_url("**")
    
    return page
