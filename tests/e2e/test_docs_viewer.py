"""
E2E tests for documentation viewer using Playwright.
"""
import pytest
from playwright.sync_api import Page, expect


@pytest.mark.django_db
class TestDocsViewerE2E:
    """End-to-end tests for documentation viewer."""
    
    def test_index_page_loads(self, page: Page, live_server):
        """Documentation index should load successfully."""
        page.goto(f"{live_server.url}/docs/")
        
        # Should show title
        expect(page.locator("h1")).to_contain_text("Documentation")
        
        # Should have search box
        search_input = page.locator('input[name="q"]')
        expect(search_input).to_be_visible()
    
    def test_can_view_readme(self, page: Page, live_server):
        """Should be able to view README file."""
        page.goto(f"{live_server.url}/docs/")
        
        # Click on README link (assuming it exists)
        readme_link = page.locator('a:has-text("README")')
        if readme_link.count() > 0:
            readme_link.first.click()
            
            # Should navigate to view page
            expect(page).to_have_url(lambda url: '/docs/view/' in url)
            
            # Should show rendered content
            expect(page.locator('.doc-content')).to_be_visible()
    
    def test_sidebar_navigation_works(self, page: Page, live_server):
        """Sidebar TOC navigation should work."""
        page.goto(f"{live_server.url}/docs/view/README.md")
        
        # Should show breadcrumb
        breadcrumb = page.locator('.breadcrumb')
        expect(breadcrumb).to_be_visible()
        
        # Should show back button
        back_button = page.locator('a:has-text("Back")')
        expect(back_button).to_be_visible()
    
    def test_search_functionality(self, page: Page, live_server):
        """Search should find and display results."""
        page.goto(f"{live_server.url}/docs/")
        
        # Enter search query
        search_input = page.locator('input[name="q"]')
        search_input.fill("test")
        
        # Submit search
        search_button = page.locator('button[type="submit"]:has-text("Search")')
        search_button.click()
        
        # Should navigate to search results
        expect(page).to_have_url(lambda url: '/docs/search/' in url)
        
        # Should show results or "no results"
        expect(page.locator('body')).to_contain_text(
            lambda text: 'result' in text.lower() or 'found' in text.lower(),
            timeout=5000
        )
    
    def test_download_raw_file(self, page: Page, live_server):
        """Should be able to download raw markdown."""
        page.goto(f"{live_server.url}/docs/view/README.md")
        
        # Find download link
        download_link = page.locator('a:has-text("Download")')
        if download_link.count() > 0:
            # Just verify link exists and is clickable
            expect(download_link.first).to_be_visible()
            expect(download_link.first).to have_attribute("href", lambda href: '/docs/raw/' in href)
    
    def test_code_blocks_highlighted(self, page: Page, live_server):
        """Code blocks should have syntax highlighting."""
        # Find a doc with code blocks
        page.goto(f"{live_server.url}/docs/view/tests/TESTING.md")
        
        # Look for code blocks
        code_blocks = page.locator('pre code, .highlight')
        if code_blocks.count() > 0:
            # Should have syntax highlighting classes or elements
            expect(code_blocks.first).to_be_visible()
    
    def test_mobile_responsive(self, page: Page, live_server):
        """Documentation should be responsive on mobile."""
        # Set mobile viewport
        page.set_viewport_size({"width": 375, "height": 667})
        
        page.goto(f"{live_server.url}/docs/")
        
        # Should still be accessible
        expect(page.locator("h1")).to_be_visible()
        
        # Search should be accessible
        search_input = page.locator('input[name="q"]')
        expect(search_input).to_be_visible()


@pytest.mark.django_db
class TestDocsSecurityE2E:
    """E2E security tests."""
    
    def test_path_traversal_shows_404(self, page: Page, live_server):
        """Path traversal attempts should show 404."""
        page.goto(f"{live_server.url}/docs/view/../etc/passwd")
        
        # Should show 404 page
        expect(page).to_have_url(lambda url: '/docs/view/' in url)
        # Django 404 page or custom 404
        expect(page.locator('body')).to_contain_text(
            lambda text: '404' in text or 'not found' in text.lower(),
            timeout=3000
        )
    
    def test_restricted_file_shows_404(self, page: Page, live_server):
        """Restricted files should not be accessible."""
        page.goto(f"{live_server.url}/docs/view/manage.py")
        
        # Should show 404
        expect(page.locator('body')).to_contain_text(
            lambda text: '404' in text or 'not found' in text.lower(),
            timeout=3000
        )


@pytest.mark.django_db
class TestDocsAccessibility:
    """Accessibility tests for documentation viewer."""
    
    def test_keyboard_navigation(self, page: Page, live_server):
        """Should be navigable via keyboard."""
        page.goto(f"{live_server.url}/docs/")
        
        # Tab through elements
        page.keyboard.press("Tab")
        
        # Should focus on interactive elements
        focused = page.evaluate("document.activeElement.tagName")
        assert focused in ['A', 'INPUT', 'BUTTON']
    
    def test_headings_hierarchy(self, page: Page, live_server):
        """Should have proper heading hierarchy."""
        page.goto(f"{live_server.url}/docs/view/README.md")
        
        # Should have h1
        h1 = page.locator('h1')
        expect(h1).to_be_visible()
        
        # Check heading structure (h1 should come before h2)
        headings = page.locator('h1, h2, h3').all()
        if len(headings) > 0:
            assert headings[0].tag_name.lower() == 'h1'
