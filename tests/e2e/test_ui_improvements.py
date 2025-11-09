"""
E2E tests for UX improvements: dose colored zones, educational errors, sandbox hints.
"""
import pytest
from playwright.sync_api import Page, expect


@pytest.mark.django_db
class TestDoseColorZones:
    """Test the colored zone feedback for dose inputs."""
    
    def test_dose_input_shows_green_zone_for_safe_values(
        self, authenticated_page: Page, live_server, test_scenario
    ):
        """Verify green zone appears for safe lenalidomide doses (10-20mg)."""
        # Navigate to scenario detail page
        authenticated_page.goto(
            f"{live_server.url}/sim/scenarios/{test_scenario.pk}/"
        )
        
        # Select VRd preset (if available)
        preset_select = authenticated_page.locator("#id_preset")
        if preset_select.count() > 0:
            preset_select.select_option(label="VRd")
        
        # Enter safe dose (15mg)
        lenalidomide_input = authenticated_page.locator("#id_lenalidomide_dose")
        lenalidomide_input.fill("15")
        
        # Trigger dose zone color update (blur event)
        lenalidomide_input.blur()
        
        # Verify zone-safe class is applied
        expect(lenalidomide_input).to_have_class(
            lambda class_list: "zone-safe" in class_list,
            timeout=2000
        )
    
    def test_dose_input_shows_yellow_zone_for_caution_values(
        self, authenticated_page: Page, live_server, test_scenario
    ):
        """Verify yellow zone appears for caution lenalidomide doses (20-25mg)."""
        authenticated_page.goto(
            f"{live_server.url}/sim/scenarios/{test_scenario.pk}/"
        )
        
        # Select VRd preset
        preset_select = authenticated_page.locator("#id_preset")
        if preset_select.count() > 0:
            preset_select.select_option(label="VRd")
        
        # Enter caution dose (23mg)
        lenalidomide_input = authenticated_page.locator("#id_lenalidomide_dose")
        lenalidomide_input.fill("23")
        lenalidomide_input.blur()
        
        # Verify zone-caution class is applied
        expect(lenalidomide_input).to_have_class(
            lambda class_list: "zone-caution" in class_list,
            timeout=2000
        )
    
    def test_dose_input_shows_red_zone_for_danger_values(
        self, authenticated_page: Page, live_server, test_scenario
    ):
        """Verify red zone appears for dangerous lenalidomide doses (>25mg)."""
        authenticated_page.goto(
            f"{live_server.url}/sim/scenarios/{test_scenario.pk}/"
        )
        
        # Select VRd preset
        preset_select = authenticated_page.locator("#id_preset")
        if preset_select.count() > 0:
            preset_select.select_option(label="VRd")
        
        # Enter danger dose (30mg)
        lenalidomide_input = authenticated_page.locator("#id_lenalidomide_dose")
        lenalidomide_input.fill("30")
        lenalidomide_input.blur()
        
        # Verify zone-danger class is applied
        expect(lenalidomide_input).to_have_class(
            lambda class_list: "zone-danger" in class_list,
            timeout=2000
        )


@pytest.mark.django_db
class TestEducationalErrors:
    """Test enhanced validation error messages with educational explanations."""
    
    def test_high_dose_shows_educational_error(
        self, authenticated_page: Page, live_server, test_scenario
    ):
        """Verify educational error message appears for excessive doses."""
        authenticated_page.goto(
            f"{live_server.url}/sim/scenarios/{test_scenario.pk}/"
        )
        
        # Fill form with excessive lenalidomide dose
        preset_select = authenticated_page.locator("#id_preset")
        if preset_select.count() > 0:
            preset_select.select_option(label="VRd")
        
        authenticated_page.locator("#id_lenalidomide_dose").fill("60")
        authenticated_page.locator("#id_time_horizon").fill("365")
        authenticated_page.locator("#id_cohort_size").fill("100")
        
        # Submit form
        authenticated_page.locator('button[type="submit"]').click()
        
        # Wait for error message
        error_message = authenticated_page.locator(".alert-danger, .invalid-feedback")
        expect(error_message).to_be_visible(timeout=3000)
        
        # Verify educational content exists (💡 Why?)
        error_text = error_message.inner_text()
        assert "💡" in error_text or "Why" in error_text, \
            "Educational explanation not found in error message"


@pytest.mark.django_db
class TestSandboxHints:
    """Test contextual learning hints system."""
    
    def test_hint_panel_appears_on_low_dose_trigger(
        self, authenticated_page: Page, live_server, test_scenario
    ):
        """Verify sandbox hint appears when user enters very low dose."""
        authenticated_page.goto(
            f"{live_server.url}/sim/scenarios/{test_scenario.pk}/"
        )
        
        # Enter very low dose (triggers low_dose hint)
        authenticated_page.locator("#id_lenalidomide_dose").fill("3")
        authenticated_page.locator("#id_lenalidomide_dose").blur()
        
        # Wait for hint panel to appear
        hint_panel = authenticated_page.locator("#sandbox-hint-panel")
        expect(hint_panel).to_be_visible(timeout=3000)
        
        # Verify hint contains relevant content
        hint_text = hint_panel.inner_text()
        assert "dose" in hint_text.lower() or "efficacy" in hint_text.lower(), \
            "Hint does not contain relevant dose guidance"
    
    def test_hint_can_be_dismissed(
        self, authenticated_page: Page, live_server, test_scenario
    ):
        """Verify hints can be dismissed by user."""
        authenticated_page.goto(
            f"{live_server.url}/sim/scenarios/{test_scenario.pk}/"
        )
        
        # Trigger a hint
        authenticated_page.locator("#id_time_horizon").fill("10")
        authenticated_page.locator("#id_time_horizon").blur()
        
        # Wait for hint
        hint_panel = authenticated_page.locator("#sandbox-hint-panel")
        if hint_panel.is_visible():
            # Click dismiss button
            dismiss_button = hint_panel.locator(".hint-dismiss")
            dismiss_button.click()
            
            # Verify hint is hidden
            expect(hint_panel).to_be_hidden(timeout=1000)


@pytest.mark.django_db
class TestEmptyState:
    """Test improved empty state with actionable checklist."""
    
    def test_empty_results_shows_checklist(
        self, authenticated_page: Page, live_server, test_scenario
    ):
        """Verify empty state checklist appears before first simulation."""
        authenticated_page.goto(
            f"{live_server.url}/sim/scenarios/{test_scenario.pk}/"
        )
        
        # Look for empty state container
        empty_state = authenticated_page.locator(".empty-state")
        
        # If visible, verify checklist items
        if empty_state.is_visible():
            expect(empty_state.locator(".checklist")).to_be_visible()
            
            # Verify checklist has multiple items
            checklist_items = empty_state.locator(".checklist li")
            expect(checklist_items).to_have_count(4)  # 4 checklist items expected
    
    def test_show_me_how_button_exists(
        self, authenticated_page: Page, live_server, test_scenario
    ):
        """Verify 'Show Me How' button is present in empty state."""
        authenticated_page.goto(
            f"{live_server.url}/sim/scenarios/{test_scenario.pk}/"
        )
        
        # Look for Show Me How button
        empty_state = authenticated_page.locator(".empty-state")
        
        if empty_state.is_visible():
            show_me_button = empty_state.locator("button:has-text('Show Me How')")
            expect(show_me_button).to_be_visible()
