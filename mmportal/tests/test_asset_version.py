"""Local asset URLs follow the image revision without trusting arbitrary env text."""

import os
from html.parser import HTMLParser
from unittest.mock import patch

from django.contrib.auth.models import AnonymousUser
from django.template.loader import render_to_string
from django.test import RequestFactory, SimpleTestCase


class AssetReferences(HTMLParser):
    def __init__(self):
        super().__init__()
        self.urls = []

    def handle_starttag(self, tag, attrs):
        if tag in {"link", "script"}:
            self.urls.extend(value for name, value in attrs if name in {"href", "src"})


class AssetVersionTests(SimpleTestCase):
    revision = "0123456789abcdef0123456789abcdef01234567"

    def render_urls(self, template):
        request = RequestFactory().get("/docs/")
        request.user = AnonymousUser()
        parser = AssetReferences()
        parser.feed(render_to_string(template, {"path": "README.md"}, request=request))
        return parser.urls

    @patch.dict(os.environ, {"BMY_SOURCE_SHA": revision})
    def test_shared_docs_and_account_assets_use_the_image_revision(self):
        templates = {
            "base.html": [],
            "docs_viewer/view.html": ["docs_viewer/reference.css", "docs_viewer/reference.js"],
            "accounts/form.html": ["accounts/access.js"],
            "accounts/security.html": ["accounts/access.js"],
            "accounts/manage.html": ["accounts/access.js"],
            "accounts/feedback.html": ["accounts/access.js"],
            "accounts/manage_feedback.html": ["accounts/access.js"],
        }
        for template, assets in templates.items():
            with self.subTest(template=template):
                urls = self.render_urls(template)
                for asset in ["app/ui-guardrails.css", *assets]:
                    matching = [url for url in urls if url.startswith(f"/static/{asset}")]
                    self.assertEqual(matching, [f"/static/{asset}?v={self.revision}"])
                external = [url for url in urls if url.startswith("https://")]
                self.assertTrue(external)
                self.assertFalse(any("?v=" in url for url in external))

    def test_missing_or_untrusted_revision_uses_stable_development_version(self):
        for revision in (
            None,
            "",
            "abc123",
            "f" * 41,
            "f" * 40 + "\n",
            "<script>alert(1)</script>",
            "../../private?secret=value",
            "F" * 40,
        ):
            with self.subTest(revision=revision), patch.dict(os.environ):
                if revision is None:
                    os.environ.pop("BMY_SOURCE_SHA", None)
                else:
                    os.environ["BMY_SOURCE_SHA"] = revision
                urls = self.render_urls("base.html")
                self.assertIn("/static/app/ui-guardrails.css?v=dev", urls)
                self.assertFalse(any("secret=value" in url or "<script>" in url for url in urls))

    def test_next_image_changes_the_asset_url(self):
        with patch.dict(os.environ, {"BMY_SOURCE_SHA": self.revision}):
            first = self.render_urls("base.html")
        with patch.dict(os.environ, {"BMY_SOURCE_SHA": "f" * 40}):
            second = self.render_urls("base.html")
        self.assertNotEqual(first, second)
        self.assertIn("/static/app/ui-guardrails.css?v=" + "f" * 40, second)
