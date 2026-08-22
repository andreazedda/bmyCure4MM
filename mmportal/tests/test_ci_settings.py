from __future__ import annotations

import os
import subprocess
import sys
import tempfile
from pathlib import Path

from django.test import SimpleTestCase


class CISettingsTests(SimpleTestCase):
    def test_ci_settings_replace_external_secrets_and_isolate_artifacts(self) -> None:
        root = Path(__file__).resolve().parents[2]
        with tempfile.TemporaryDirectory(prefix="bmycure4mm-ci-settings-") as temp_dir:
            env = os.environ.copy()
            env.update(
                {
                    "BMYCURE4MM_CI_TEMP_ROOT": temp_dir,
                    "DJANGO_SECRET_KEY": "external-production-secret-must-not-be-imported",
                    "DJANGO_SETTINGS_MODULE": "mmportal.settings_ci",
                }
            )
            script = """
from pathlib import Path
from mmportal import settings_ci, settings_ci_deploy

assert settings_ci.SECRET_KEY == settings_ci.CI_SYNTHETIC_SECRET_KEY
assert settings_ci.SECRET_KEY != "external-production-secret-must-not-be-imported"
assert settings_ci_deploy.SECRET_KEY == settings_ci.CI_SYNTHETIC_SECRET_KEY
assert settings_ci.DEBUG is False
assert settings_ci_deploy.DEBUG is False
assert settings_ci_deploy.SECURE_SSL_REDIRECT is True
assert Path(settings_ci.DATABASES["default"]["NAME"]).is_relative_to(settings_ci.CI_TEMP_ROOT)
assert Path(settings_ci.MEDIA_ROOT).is_relative_to(settings_ci.CI_TEMP_ROOT)
assert Path(settings_ci.LOGS_DIR).is_relative_to(settings_ci.CI_TEMP_ROOT)
"""
            result = subprocess.run(
                [sys.executable, "-c", script],
                cwd=root,
                env=env,
                check=False,
                capture_output=True,
                text=True,
            )
        self.assertEqual(result.returncode, 0, result.stdout + result.stderr)
