"""Optional real-browser checks for public documentation and contextual help."""
import os
import re
import json
from pathlib import Path

from django.contrib.staticfiles.testing import StaticLiveServerTestCase
from django.test import override_settings
from playwright.sync_api import expect, sync_playwright

from simulator.api_help import GOVERNED_STATIC_HELP
from simulator.models_help import HelpArticle


@override_settings(DEBUG=False, SECURE_SSL_REDIRECT=False)
class DocumentationAccessibilityJourney(StaticLiveServerTestCase):
    host = '127.0.0.1'

    def setUp(self):
        super().setUp()
        # LiveServer tests flush migration-seeded rows between classes. Declare
        # the actual search index record so this journey is independent of order.
        english = GOVERNED_STATIC_HELP['quickstart']['en']
        italian = GOVERNED_STATIC_HELP['quickstart']['it']
        HelpArticle.objects.update_or_create(
            slug='quickstart',
            defaults={
                'title_en': english[0], 'body_en': english[1],
                'title_it': italian[0], 'body_it': italian[1],
            },
        )

    def test_keyboard_reflow_help_recovery_and_contrast(self):
        errors = []
        metrics = {}
        captures = os.environ.get('BMY_AUDIT_CAPTURE_DIR')

        def capture(page, name):
            if captures:
                root = Path(captures)
                root.mkdir(parents=True, exist_ok=True)
                page.screenshot(path=str(root / (name + '.png')), full_page=False)

        def contrast(rgb, background):
            if len(rgb) == 4:
                rgb = [rgb[i] * rgb[3] + background[i] * (1 - rgb[3]) for i in range(3)]
            def luminance(color):
                values = [v / 255 for v in color[:3]]
                return sum(w * (v / 12.92 if v <= .04045 else ((v + .055) / 1.055) ** 2.4)
                           for w, v in zip((.2126, .7152, .0722), values))
            a, b = sorted((luminance(rgb), luminance(background)))
            return (b + .05) / (a + .05)

        with sync_playwright() as playwright:
            browser = playwright.chromium.launch()
            context = browser.new_context(viewport={'width': 1280, 'height': 900}, reduced_motion='reduce')
            page = context.new_page()
            page.on('pageerror', lambda error: errors.append(str(error)))
            page.goto(self.live_server_url + '/docs/')
            page.keyboard.press('Tab')
            expect(page.locator('.skip-link')).to_be_focused()
            page.keyboard.press('Enter')
            expect(page.locator('#main-content')).to_be_focused()
            colors = page.locator('.navbar .nav-link').first.evaluate("el => [getComputedStyle(el).color.match(/[\\d.]+/g).map(Number), getComputedStyle(document.querySelector('.navbar')).backgroundColor.match(/[\\d.]+/g).map(Number)]")
            metrics.setdefault('navbar_contrast', round(contrast(*colors), 3))
            self.assertGreaterEqual(contrast(*colors), 4.5)
            nav_link = page.locator('.navbar .nav-link').first
            nav_link.hover()
            colors = nav_link.evaluate("el => [getComputedStyle(el).color.match(/[\\d.]+/g).map(Number), getComputedStyle(el).backgroundColor.match(/[\\d.]+/g).map(Number)]")
            metrics['navbar_hover_contrast'] = round(contrast(*colors), 3)
            self.assertGreaterEqual(contrast(*colors), 4.5)
            page.locator('#doc-search').focus()
            page.keyboard.press('Control+k')
            expect(page.locator('#cmdk')).to_be_focused()
            page.locator('#cmdk').fill('quickstart')
            help_result = page.locator('#cmdkResults [data-help="quickstart"]')
            expect(help_result).to_be_visible()
            help_result.press('Enter')
            expect(page.locator('#helpDrawer')).to_have_class(re.compile(r'\bshow\b'))
            expect(page.locator('#helpDrawerContent')).to_contain_text('synthetic scenario')
            page.keyboard.press('Escape')
            expect(page.locator('#helpDrawer')).not_to_have_class(re.compile(r'\bshow\b'))
            expect(page.locator('#doc-search')).to_be_focused()
            page.keyboard.press('Control+k')
            page.locator('#cmdk').fill('zz-no-help-result-9911')
            expect(page.locator('#cmdkResults')).to_contain_text('No help found')
            page.route('**/api/help/search/**', lambda route: route.fulfill(status=503, body='Unavailable'))
            page.locator('#cmdk').fill('temporary failure')
            expect(page.locator('#cmdkResults')).to_contain_text('temporarily unavailable')
            page.unroute('**/api/help/search/**')
            page.locator('#cmdk').fill('quickstart')
            expect(help_result).to_be_visible()
            page.set_viewport_size({'width': 320, 'height': 700})
            self.assertLessEqual(page.locator('#cmdk').bounding_box()['width'], 320)
            capture(page, 'viewer-01-help-search-mobile')
            page.keyboard.press('Escape')
            page.goto(self.live_server_url + '/docs/view/docs/development/DEVELOPMENT.md')
            expect(page.locator('.doc-navigation')).not_to_have_attribute('open', '')
            expect(page.locator('.doc-metadata')).to_contain_text('Authored status:')
            capture(page, 'viewer-02-handbook-320px')
            metrics['heading_top_320px'] = page.locator('.doc-content h1').first.bounding_box()['y']
            self.assertLess(metrics['heading_top_320px'], 700)
            self.assertFalse(page.evaluate('document.documentElement.scrollWidth > innerWidth + 1'))
            page.locator('.doc-navigation summary').focus()
            page.keyboard.press('Enter')
            expect(page.locator('.doc-navigation')).to_have_attribute('open', '')
            page.keyboard.press('Enter')
            expect(page.locator('.doc-navigation')).not_to_have_attribute('open', '')
            table = page.locator('.doc-content .doc-table-scroll').first
            expect(table).to_have_attribute('tabindex', '0')
            table.focus()
            page.keyboard.press('ArrowRight')
            page.wait_for_function("document.querySelector('.doc-content .doc-table-scroll').scrollLeft > 0")
            self.assertGreater(table.evaluate('el => el.scrollLeft'), 0)
            capture(page, 'viewer-04-table-keyboard-scroll')
            page.set_viewport_size({'width': 1280, 'height': 900})
            expect(page.locator('.doc-navigation')).to_have_attribute('open', '')
            colors = page.locator('.doc-content .highlight pre').first.evaluate("""pre => {
                const sample = document.createElement('span'); sample.className = 'c1'; sample.textContent = '# contrast test'; pre.append(sample);
                const result = [getComputedStyle(sample).color.match(/[\\d.]+/g).map(Number), getComputedStyle(pre).backgroundColor.match(/[\\d.]+/g).map(Number)];
                sample.remove(); return result;
            }""")
            metrics['code_comment_contrast'] = round(contrast(*colors), 3)
            self.assertGreaterEqual(contrast(*colors), 4.5)
            page.set_viewport_size({'width': 640, 'height': 900})
            expect(page.locator('.doc-navigation')).not_to_have_attribute('open', '')
            self.assertFalse(page.evaluate('document.documentElement.scrollWidth > innerWidth + 1'))
            page.locator('#lang-it').click()
            expect(page.locator('#lang-it')).to_have_attribute('aria-pressed', 'true')
            expect(page.locator('#lang-en')).to_have_attribute('aria-pressed', 'false')
            page.locator('.doc-copy').first.click()
            expect(page.locator('#live-region')).to_contain_text(re.compile('Copiato|Codice selezionato'))
            capture(page, 'viewer-03-code-copy-status')
            self.assertEqual(errors, [])
            print('DOCS_BROWSER_METRICS=' + json.dumps(metrics, sort_keys=True))
            browser.close()
