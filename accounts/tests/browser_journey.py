"""Real forms, sessions and SMTP in a disposable database; no production identities."""

import re
import secrets
import socket
from concurrent.futures import ThreadPoolExecutor
from email import message_from_bytes

import pyotp
from aiosmtpd.controller import Controller
from django.contrib.staticfiles.testing import StaticLiveServerTestCase
from django.db import connection
from django.test import override_settings
from playwright.sync_api import expect, sync_playwright

from accounts import services
from accounts.models import AccessProfile, Notification
from accounts.notifications import dispatch_notification

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
                body = message.get_payload(decode=True).decode()
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

        def set_password(page, link):
            page.goto(link)
            expect(page.locator('input[name="token"]')).not_to_have_value("")
            self.assertNotIn("#", page.url)
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
                    key = owner_page.locator("code").inner_text()
                    owner_page.get_by_label("Code", exact=True).fill(pyotp.TOTP(key).now())
                    owner_page.get_by_role("button", name="Confirm authenticator code").click()
                    expect(
                        owner_page.get_by_role("heading", name="Access requests", exact=True)
                    ).to_be_visible()

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
                        page.get_by_label("Email", exact=True).fill(email)
                        page.get_by_label("Nickname", exact=True).fill(f"Browser researcher {index}")
                        page.get_by_label("Reason", exact=True).fill(
                            "Synthetic browser verification of the existing platform"
                        )
                        page.locator('input[name="acknowledgement"]').check()
                        page.get_by_role("button", name="Send verification email").click()
                        expect(page.get_by_role("heading", name="Check your email")).to_be_visible()
                        deliver()
                        verification = sink.link(email, "/access/verify/")
                        page.goto(verification)
                        self.assertEqual(
                            db_call(lambda: AccessProfile.objects.get(email=email).request_state),
                            "PENDING_EMAIL_VERIFICATION",
                        )
                        page.get_by_role("button", name="Confirm email").click()
                        expect(page.get_by_text("Waiting for review", exact=True)).to_be_visible()
                        deliver()
                        owner_page.goto(self.live_server_url + "/manage/access-requests/")
                        card = owner_page.locator("article[data-access-request]").filter(has_text=email)
                        card.get_by_role("button", name="Approve", exact=True).click()
                        deliver()
                        set_password(page, sink.link(email, "/access/activate/"))
                        login(page, email)
                        expect(page).to_have_url(self.live_server_url + "/sim/getting-started/")
                        self.assertFalse(
                            db_call(lambda: AccessProfile.objects.get(email=email).user.is_staff)
                        )
                        expect(
                            page.get_by_text("Explore the learning scenarios", exact=False)
                        ).to_be_visible()
                        page.goto(self.live_server_url + "/feedback/")
                        page.get_by_label("Category", exact=True).select_option("usability")
                        message = f"Synthetic suggestion {index}: make the first steps clearer."
                        page.get_by_label("Message", exact=True).fill(message)
                        page.get_by_role("button", name="Send ticket").click()
                        expect(page.get_by_text(message, exact=True)).to_be_visible()
                        if index:
                            expect(page.get_by_text("Synthetic suggestion 0:", exact=False)).to_have_count(0)
                        owner_page.goto(self.live_server_url + "/manage/feedback/")
                        ticket = owner_page.locator("article").filter(has_text=message)
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
                        self.assertLessEqual(
                            page.evaluate("document.documentElement.scrollWidth"), viewport["width"]
                        )
                        logout(page)
                        login(page, email)
                        page.goto(self.live_server_url + "/feedback/")
                        expect(page.get_by_text(message, exact=True)).to_be_visible()
                        self.assertEqual(page.goto(self.live_server_url + "/patients/").status, 403)
                        logout(page)
                        context.close()
                    owner_page.set_viewport_size({"width": 390, "height": 844})
                    owner_page.goto(self.live_server_url + "/manage/feedback/")
                    expect(
                        owner_page.get_by_role("heading", name="User tickets and suggestions")
                    ).to_be_visible()
                    self.assertLessEqual(owner_page.evaluate("document.documentElement.scrollWidth"), 390)
                    self.assertFalse(errors, "Browser JavaScript errors")
                    self.assertFalse(bad_assets, "Missing local CSS or JavaScript")
                    self.assertTrue(db_call(lambda: Notification.objects.filter(state="sent").exists()))
                    browser.close()
        finally:
            smtp.stop()
