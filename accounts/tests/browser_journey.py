"""Real forms, sessions and SMTP in a disposable database; no production identities."""

import io
import os
import re
import secrets
import socket
import traceback
import uuid
from concurrent.futures import ThreadPoolExecutor
from email import message_from_bytes
from pathlib import Path

import pyotp
import zxingcpp
from aiosmtpd.controller import Controller
from django.contrib.staticfiles.testing import StaticLiveServerTestCase
from django.db import connection
from django.test import override_settings
from PIL import Image
from playwright.sync_api import expect, sync_playwright

from accounts import services
from accounts.models import AccessProfile, Notification, SupportTicket
from accounts.notifications import dispatch_notification
from simulator.models import Scenario, SimulationAttempt

from .helpers import TEST_SETTINGS, open_policy


class MailSink:
    def __init__(self):
        self.messages = []

    async def handle_DATA(self, server, session, envelope):
        self.messages.append((envelope.rcpt_tos, message_from_bytes(envelope.content)))
        return "250 Message accepted by isolated test sink"

    def link(self, email, path):
        for recipients, message in reversed(self.messages):
            if email in recipients:
                plain = next(p for p in message.walk() if p.get_content_type() == "text/plain")
                body = plain.get_payload(decode=True).decode(plain.get_content_charset() or "utf-8")
                match = re.search(r"http://127\.0\.0\.1:\d+" + path + r"#token=[A-Za-z0-9_-]+", body)
                if match:
                    return match.group()
        raise AssertionError("Expected message was not received by the isolated sink")


@override_settings(**TEST_SETTINGS)
class RegistrationBrowserJourney(StaticLiveServerTestCase):
    host = "127.0.0.1"

    def test_desktop_and_mobile_existing_site(self):
        with socket.socket() as sock:
            sock.bind(("127.0.0.1", 0))
            smtp_port = sock.getsockname()[1]
        sink = MailSink()
        smtp = Controller(sink, hostname="127.0.0.1", port=smtp_port)
        smtp.start()
        password = secrets.token_urlsafe(24)
        errors, bad_assets = [], []
        stage = "account setup"
        capture_root = os.environ.get("BMY_AUDIT_CAPTURE_DIR")

        def capture(page, name):
            if capture_root:
                target = Path(capture_root)
                target.mkdir(parents=True, exist_ok=True)
                page.screenshot(path=str(target / (name + ".png")), full_page=True, animations="disabled")

        def check_accessibility(page, width):
            self.assertLessEqual(page.evaluate("document.documentElement.scrollWidth"), width)
            self.assertTrue(
                page.evaluate(
                    """() => [...document.querySelectorAll('input:not([type=hidden]),select,textarea')].filter(e=>e.getClientRects().length).every(e=>e.labels.length || e.getAttribute('aria-label'))"""
                ),
                "Visible field has no accessible label",
            )
            self.assertTrue(
                page.evaluate(
                    """() => [...document.querySelectorAll('[aria-describedby]')].every(e=>e.getAttribute('aria-describedby').split(/\s+/).every(id=>document.getElementById(id)))"""
                ),
                "Broken field description reference",
            )

        def db_call(fn):
            def operation():
                try:
                    return fn()
                finally:
                    connection.close()

            return executor.submit(operation).result(timeout=20)

        def deliver():
            while db_call(dispatch_notification):
                pass

        def login(page, email):
            page.goto(self.live_server_url + "/accounts/login/")
            page.get_by_label("Email", exact=True).fill(email)
            page.get_by_label("Password", exact=True).fill(password)
            page.get_by_role("button", name="Sign in", exact=True).click()

        def logout(page):
            toggle = page.get_by_role("button", name="Toggle navigation")
            if toggle.is_visible():
                toggle.click()
            page.locator("#userMenu").click()
            page.get_by_role("button", name="Log out", exact=True).click()

        def set_password(page, link, capture_name=None):
            page.goto(link)
            expect(page.locator('input[name="token"]')).not_to_have_value("")
            self.assertNotIn("#", page.url)
            if capture_name:
                expect(page.get_by_label("Password", exact=True)).to_have_value("")
                capture(page, capture_name)
            page.get_by_label("Password", exact=True).fill(password)
            page.get_by_label("Password confirm", exact=True).fill(password)
            page.get_by_role("button", name="Save password").click()
            expect(page.get_by_role("heading", name="Sign in", exact=True)).to_be_visible()

        try:
            with (
                self.settings(
                    DEBUG=True,
                    PUBLIC_BASE_URL=self.live_server_url,
                    EMAIL_BACKEND="django.core.mail.backends.smtp.EmailBackend",
                    EMAIL_HOST="127.0.0.1",
                    EMAIL_PORT=smtp_port,
                    EMAIL_USE_TLS=False,
                    EMAIL_USE_SSL=False,
                    EMAIL_HOST_USER="",
                    EMAIL_HOST_PASSWORD="",
                ),
                ThreadPoolExecutor(max_workers=1) as executor,
            ):
                open_policy()
                services.bootstrap_owner()
                scenario_id = db_call(
                    lambda: Scenario.objects.create(
                        title="Synthetic browser learning example",
                        clinical_stage="newly_diagnosed",
                        summary="Independent synthetic software verification",
                        risk_stratification="Standard",
                    ).pk
                )
                with sync_playwright() as playwright:
                    browser = playwright.chromium.launch()
                    owner_context = browser.new_context(viewport={"width": 1280, "height": 900})
                    owner_page = owner_context.new_page()
                    deliver()
                    set_password(
                        owner_page, sink.link(TEST_SETTINGS["PLATFORM_OWNER_EMAIL"], "/owner/activate/")
                    )
                    login(owner_page, TEST_SETTINGS["PLATFORM_OWNER_EMAIL"])
                    owner_page.goto(self.live_server_url + "/manage/access-requests/")
                    expect(owner_page.get_by_role("heading", name="Owner security")).to_be_visible()
                    decoded_uris = []
                    for size in ({"width": 1280, "height": 900}, {"width": 390, "height": 844}):
                        owner_page.set_viewport_size(size)
                        qr = owner_page.get_by_alt_text("Authenticator setup QR code")
                        expect(qr).to_be_visible()
                        decoded = zxingcpp.read_barcodes(Image.open(io.BytesIO(qr.screenshot())))
                        self.assertEqual(len(decoded), 1)
                        decoded_uris.append(decoded[0].text)
                    self.assertTrue(
                        decoded_uris[0] == decoded_uris[1], "Setup identity changed between viewports"
                    )
                    authenticator = pyotp.parse_uri(decoded_uris[0])
                    self.assertEqual(authenticator.issuer, "bmyCure4MM")
                    owner_page.get_by_label("Code", exact=True).fill(authenticator.now())
                    owner_page.get_by_role("button", name="Confirm authenticator code").click()
                    owner_page.set_viewport_size({"width": 1280, "height": 900})
                    expect(
                        owner_page.get_by_role("heading", name="Access requests", exact=True)
                    ).to_be_visible()

                    capture(owner_page, "01-owner-console-desktop")
                    first_run = None
                    for index, viewport in enumerate(
                        ({"width": 1440, "height": 1000}, {"width": 390, "height": 844})
                    ):
                        email = f"browser{index}@example.test"
                        context = browser.new_context(viewport=viewport)
                        page = context.new_page()
                        page.on("pageerror", lambda error: errors.append(str(error)))
                        page.on(
                            "response",
                            lambda response: bad_assets.append(response.status)
                            if "/static/" in response.url and response.status >= 400
                            else None,
                        )
                        page.goto(self.live_server_url + "/access/request/")
                        page.keyboard.press("Tab")
                        expect(page.locator(".skip-link")).to_be_focused()
                        page.keyboard.press("Enter")
                        expect(page.locator("#main-content")).to_be_focused()
                        check_accessibility(page, viewport["width"])
                        capture(page, f"02-request-{index}")
                        page.get_by_label("Email", exact=True).fill(email)
                        page.get_by_label("Nickname", exact=True).fill(f"Browser researcher {index}")
                        page.get_by_label("Reason", exact=True).fill(
                            "Synthetic browser verification of the existing platform"
                        )
                        page.locator('input[name="acknowledgement"]').check()
                        page.get_by_role("button", name="Send verification email").click()
                        expect(page.get_by_role("heading", name="Check your email")).to_be_visible()
                        capture(page, f"03-submitted-{index}")
                        deliver()
                        verification = sink.link(email, "/access/verify/")
                        page.goto(verification)
                        self.assertEqual(
                            db_call(lambda: AccessProfile.objects.get(email=email).request_state),
                            "PENDING_EMAIL_VERIFICATION",
                        )
                        page.get_by_role("button", name="Confirm email").click()
                        expect(page.get_by_text("Waiting for review", exact=True)).to_be_visible()
                        capture(page, f"04-awaiting-review-{index}")
                        deliver()
                        owner_page.goto(self.live_server_url + "/manage/access-requests/")
                        card = owner_page.locator("article[data-access-request]").filter(has_text=email)
                        capture(owner_page, f"05-owner-review-{index}")
                        card.get_by_role("button", name="Approve", exact=True).click()
                        deliver()
                        set_password(
                            page, sink.link(email, "/access/activate/"), f"06a-activation-empty-{index}"
                        )
                        page.get_by_label("Email", exact=True).fill(email)
                        page.get_by_label("Password", exact=True).fill("synthetic-wrong-password")
                        page.get_by_role("button", name="Sign in", exact=True).click()
                        expect(page.locator("[data-error-summary]")).to_be_focused()
                        expect(page.get_by_label("Password", exact=True)).to_have_value("")
                        capture(page, f"11-sign-in-error-{index}")
                        login(page, email)
                        expect(page).to_have_url(self.live_server_url + "/sim/getting-started/")
                        self.assertFalse(
                            db_call(lambda: AccessProfile.objects.get(email=email).user.is_staff)
                        )
                        expect(
                            page.get_by_text("Explore the learning scenarios", exact=False)
                        ).to_be_visible()
                        capture(page, f"06-getting-started-{index}")
                        stage = "real simulation and authorized artifacts"
                        page.goto(self.live_server_url + "/learn/")
                        page.get_by_role(
                            "link", name="Synthetic browser learning example", exact=True
                        ).click()
                        for dose in (20, 25):
                            if not page.locator("#simulateScenarioPanel").is_visible():
                                page.get_by_role("link", name="Open simulation form", exact=True).click()
                            check_accessibility(page, viewport["width"])
                            expect(page.locator('[name="twin_assessment_id"]')).to_have_count(0)
                            expect(page.locator('[name="use_twin"]')).to_have_count(0)
                            expect(page.get_by_text("Load demo data (DEMO)", exact=False)).to_have_count(0)
                            page.locator('[name="lenalidomide_dose"]').fill(str(dose))
                            page.locator('[name="time_horizon"]').fill("30")
                            with page.expect_response(
                                lambda r: "/simulate/" in r.url and r.request.method == "POST"
                            ) as computed:
                                page.get_by_role("button", name="Run simulation", exact=True).click()
                            self.assertEqual(computed.value.status, 200)
                            expect(page.get_by_role("link", name="Export study", exact=True)).to_be_visible()
                            page.reload()
                        stage = "input classification disclosure and keyboard scrolling"
                        classification = page.locator("#input-classification")
                        summary = classification.locator("summary")
                        self.assertFalse(classification.evaluate("e => e.open"))
                        summary.focus()
                        page.keyboard.press("Enter")
                        self.assertTrue(classification.evaluate("e => e.open"))
                        expect(classification.locator("tbody tr")).to_have_count(8)
                        scroller = classification.get_by_role("region")
                        scroller.focus()
                        if viewport["width"] == 390:
                            self.assertTrue(scroller.evaluate("e => e.scrollWidth > e.clientWidth"))
                            page.keyboard.press("ArrowRight")
                            page.wait_for_function(
                                "document.querySelector('#input-classification [role=region]').scrollLeft > 0"
                            )
                        check_accessibility(page, viewport["width"])
                        capture(page, f"07a-input-classification-{index}")
                        summary.focus()
                        page.keyboard.press("Enter")
                        self.assertFalse(classification.evaluate("e => e.open"))
                        capture(page, f"07-results-{index}")
                        runs = db_call(
                            lambda: list(
                                SimulationAttempt.objects.filter(user__email=email, scenario_id=scenario_id)
                                .order_by("pk")
                                .values("pk", "parameters", "results_summary")
                            )
                        )
                        self.assertEqual(len(runs), 2)
                        self.assertEqual([r["parameters"]["lenalidomide_dose"] for r in runs], [20, 25])
                        self.assertNotEqual(runs[0]["results_summary"], runs[1]["results_summary"])
                        self.assertEqual(
                            db_call(lambda: Scenario.objects.get(pk=scenario_id).title),
                            "Synthetic browser learning example",
                        )
                        stage = "authorized plot and CSV export"
                        plot_url = page.get_by_role("link", name="View plot", exact=True).get_attribute(
                            "href"
                        )
                        plot_page = context.new_page()
                        self.assertEqual(plot_page.goto(self.live_server_url + plot_url).status, 200)
                        expect(plot_page.locator(".plotly-graph-div")).to_be_visible()
                        plot_page.close()
                        export_url = page.get_by_role("link", name="Export study", exact=True).get_attribute(
                            "href"
                        )
                        self.assertEqual(page.request.get(self.live_server_url + export_url).status, 200)
                        with page.expect_download() as download:
                            page.get_by_role("link", name="Download CSV", exact=True).click()
                        self.assertTrue(download.value.suggested_filename.endswith(".csv"))
                        if index:
                            self.assertIsNotNone(first_run)
                            forbidden = self.live_server_url + f"/sim/attempts/{first_run}/export/"
                            self.assertEqual(page.request.get(forbidden).status, 403)
                        else:
                            first_run = runs[0]["pk"]
                        page.goto(self.live_server_url + "/feedback/")
                        capture(page, f"08-feedback-{index}")
                        page.get_by_label("Category", exact=True).select_option("usability")
                        message = f"Synthetic suggestion {index}: make the first steps clearer."
                        page.get_by_label("Message", exact=True).fill(message)
                        page.get_by_role("button", name="Send ticket").click()
                        expect(page.get_by_text(message, exact=True)).to_be_visible()
                        if index:
                            expect(page.get_by_text("Synthetic suggestion 0:", exact=False)).to_have_count(0)
                        owner_page.goto(self.live_server_url + "/manage/feedback/")
                        ticket = owner_page.locator("article").filter(has_text=message)
                        check_accessibility(owner_page, 1280)
                        capture(owner_page, f"09-owner-ticket-{index}")
                        ticket.get_by_label("Reply to author").fill(
                            "Thank you. Your suggestion has been received."
                        )
                        ticket.get_by_label("Internal note").fill("Private operator note")
                        ticket.get_by_label("Status", exact=True).select_option("answered")
                        ticket.get_by_role("button", name="Save reply and status").click()
                        page.reload()
                        expect(
                            page.get_by_text("Thank you. Your suggestion has been received.", exact=True)
                        ).to_be_visible()
                        expect(page.get_by_text("Private operator note", exact=True)).to_have_count(0)
                        capture(page, f"15-user-ticket-reply-{index}")
                        self.assertLessEqual(
                            page.evaluate("document.documentElement.scrollWidth"), viewport["width"]
                        )
                        logout(page)
                        login(page, email)
                        page.goto(self.live_server_url + "/feedback/")
                        expect(page.get_by_text(message, exact=True)).to_be_visible()
                        self.assertEqual(page.goto(self.live_server_url + "/patients/").status, 403)
                        logout(page)
                        capture(page, f"16-logged-out-{index}")
                        context.close()
                    owner_page.set_viewport_size({"width": 390, "height": 844})
                    owner_page.goto(self.live_server_url + "/manage/feedback/")
                    expect(
                        owner_page.get_by_role("heading", name="User tickets and suggestions")
                    ).to_be_visible()
                    self.assertLessEqual(owner_page.evaluate("document.documentElement.scrollWidth"), 390)
                    capture(owner_page, "10-owner-tickets-mobile")
                    owner_page.get_by_role("button", name="IT", exact=True).click()
                    expect(
                        owner_page.get_by_role("heading", name="Ticket e suggerimenti degli utenti")
                    ).to_be_visible()
                    expect(owner_page.get_by_label("Risposta all’autore").first).to_be_visible()
                    expect(owner_page.locator("select[name=status] option:checked").first).to_have_text(
                        "Con risposta"
                    )
                    check_accessibility(owner_page, 390)
                    capture(owner_page, "12-owner-tickets-mobile-italian")
                    owner_page.get_by_role("button", name="EN", exact=True).click()
                    stage = "paginated ticket history"
                    db_call(
                        lambda: SupportTicket.objects.bulk_create(
                            [
                                SupportTicket(
                                    identity=AccessProfile.objects.get(email="browser0@example.test"),
                                    category="bug",
                                    message=f"Synthetic pagination ticket {i}",
                                    submission_key=uuid.uuid4(),
                                )
                                for i in range(21)
                            ]
                        )
                    )
                    owner_page.goto(self.live_server_url + "/manage/feedback/")
                    owner_page.get_by_role("link", name="Next", exact=True).click()
                    expect(owner_page.get_by_text("Synthetic suggestion 0:", exact=False)).to_be_visible()
                    check_accessibility(owner_page, 390)
                    capture(owner_page, "14-owner-pagination-mobile")
                    stage = "documentation and localized forms"
                    for width in (1280, 390):
                        docs_context = browser.new_context(viewport={"width": width, "height": 844})
                        docs_page = docs_context.new_page()
                        docs_page.goto(self.live_server_url + "/docs/")
                        docs_page.get_by_role("link", name="Open the handbook", exact=True).click()
                        expect(docs_page.locator(".mermaid svg").first).to_be_visible(timeout=20000)
                        expect(docs_page.locator("mjx-container").first).to_be_visible(timeout=20000)
                        expect(docs_page.get_by_role("button", name="Copy code").first).to_be_visible()
                        if width == 390:
                            self.assertFalse(docs_page.locator(".doc-navigation").evaluate("e=>e.open"))
                            self.assertLess(docs_page.locator(".doc-content h1").bounding_box()["y"], 900)
                        docs_page.get_by_role("button", name="Copy code").first.click()
                        expect(docs_page.locator("#live-region")).to_contain_text(
                            re.compile("Copied|Code selected")
                        )
                        self.assertTrue(
                            docs_page.evaluate(
                                "[...document.querySelectorAll('.doc-content img')].every(i=>i.complete && i.naturalWidth>0)"
                            )
                        )
                        self.assertLessEqual(
                            docs_page.evaluate("document.documentElement.scrollWidth"), width
                        )
                        docs_page.goto(self.live_server_url + "/accounts/login/")
                        docs_page.get_by_role("button", name="IT", exact=True).click()
                        expect(docs_page.locator("#lang-it")).to_have_attribute("aria-pressed", "true")
                        expect(docs_page.get_by_role("heading", name="Accedi", exact=True)).to_be_visible()
                        expect(docs_page.get_by_label("Email", exact=True)).to_be_visible()
                        docs_page.get_by_role("link", name="Richiedi accesso", exact=True).click()
                        expect(docs_page.get_by_label("Motivo della richiesta", exact=True)).to_be_visible()
                        self.assertEqual(docs_page.locator("html").get_attribute("lang"), "it")
                        expect(
                            docs_page.locator('select[name="professional_role"] option[value="student"]')
                        ).to_have_text("Studente")
                        docs_page.locator("#lang-it").focus()
                        docs_page.keyboard.press("Control+k")
                        expect(docs_page.locator("#cmdk")).to_be_focused()
                        docs_page.locator("#cmdk").fill("no-such-synthetic-help-key")
                        expect(docs_page.locator("#cmdkResults")).to_contain_text("Nessuna voce trovata")
                        docs_page.route(
                            "**/api/help/search/**",
                            lambda route: route.fulfill(status=503, body="Unavailable"),
                        )
                        docs_page.locator("#cmdk").fill("test service failure")
                        expect(docs_page.locator("#cmdkResults")).to_contain_text(
                            "Aiuto temporaneamente non disponibile"
                        )
                        docs_page.unroute("**/api/help/search/**")
                        docs_page.keyboard.press("Escape")
                        expect(docs_page.locator("#lang-it")).to_be_focused()
                        check_accessibility(docs_page, width)
                        capture(docs_page, f"13-request-italian-{width}")
                        self.assertLessEqual(
                            docs_page.evaluate("document.documentElement.scrollWidth"), width
                        )
                        docs_context.close()
                    self.assertFalse(errors, "Browser JavaScript errors")
                    self.assertFalse(bad_assets, "Missing local CSS or JavaScript")
                    self.assertTrue(db_call(lambda: Notification.objects.filter(state="sent").exists()))
                    browser.close()
        except Exception as exc:
            frame = next(
                f for f in reversed(traceback.extract_tb(exc.__traceback__)) if f.filename == __file__
            )
            # Playwright diagnostics can include capability URLs or entered credentials.
            raise AssertionError(
                f"Isolated browser journey failed during {stage}: {type(exc).__name__} at line {frame.lineno}"
            ) from None
        finally:
            smtp.stop()
