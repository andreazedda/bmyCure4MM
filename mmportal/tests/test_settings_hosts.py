from __future__ import annotations

import os
from unittest.mock import patch

from django.test import SimpleTestCase

from mmportal import settings as settings_module


class HostSettingsParsingTests(SimpleTestCase):
    def test_parse_csv_env_strips_ignores_empty_and_deduplicates(self) -> None:
        with patch.dict(os.environ, {"DJANGO_ALLOWED_HOSTS": " localhost, ,127.0.0.1, localhost ,100.64.0.1 "}, clear=False):
            self.assertEqual(
                settings_module._parse_csv_env("DJANGO_ALLOWED_HOSTS"),
                ["localhost", "127.0.0.1", "100.64.0.1"],
            )

    def test_allowed_hosts_include_env_hosts_and_defaults(self) -> None:
        with patch.dict(
            os.environ,
            {
                "DJANGO_ALLOWED_HOSTS": "localhost,127.0.0.1,0.0.0.0,100.64.0.1",
                "ALLOWED_HOSTS": "",
            },
            clear=False,
        ):
            hosts = settings_module._build_allowed_hosts()

        self.assertIn("localhost", hosts)
        self.assertIn("127.0.0.1", hosts)
        self.assertIn("0.0.0.0", hosts)
        self.assertIn("100.64.0.1", hosts)

    def test_csrf_trusted_origins_parse_env_values(self) -> None:
        with patch.dict(
            os.environ,
            {
                "DJANGO_CSRF_TRUSTED_ORIGINS": "http://localhost:8001, ,http://127.0.0.1:8001,http://100.64.0.1:8001",
                "CSRF_TRUSTED_ORIGINS": "",
            },
            clear=False,
        ):
            origins = settings_module._build_csrf_trusted_origins()

        self.assertEqual(
            origins,
            ["http://localhost:8001", "http://127.0.0.1:8001", "http://100.64.0.1:8001"],
        )

    def test_default_allowed_hosts_keep_localhost(self) -> None:
        with patch.dict(os.environ, {"DJANGO_ALLOWED_HOSTS": "", "ALLOWED_HOSTS": ""}, clear=False):
            hosts = settings_module._build_allowed_hosts()

        self.assertIn("localhost", hosts)
        self.assertIn("127.0.0.1", hosts)

    def test_no_wildcard_is_added_by_default(self) -> None:
        with patch.dict(os.environ, {"DJANGO_ALLOWED_HOSTS": "", "ALLOWED_HOSTS": ""}, clear=False):
            hosts = settings_module._build_allowed_hosts()

        self.assertNotIn("*", hosts)
