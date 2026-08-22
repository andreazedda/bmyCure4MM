"""Synthetic CI-only settings for integrity and deployment checks.

These settings isolate runtime artifacts and exercise Django's deployment
checks. They are not a production architecture and must never consume a
production secret or external service credential.
"""

from __future__ import annotations

import os
import tempfile
from pathlib import Path

CI_SYNTHETIC_SECRET_KEY = (
    "ci-only-synthetic-key-for-deployment-checks-9Fj3Kp7Qm2Vx8Nz4Rw6Ty1Ua5Bc0De"
)
CI_TEMP_ROOT = Path(
    os.environ.get(
        "BMYCURE4MM_CI_TEMP_ROOT",
        str(Path(tempfile.gettempdir()) / "bmycure4mm-ci"),
    )
)
CI_TEMP_ROOT.mkdir(parents=True, exist_ok=True)

# Override, rather than inherit, secret-bearing environment variables before
# importing the base module. This is intentionally deterministic and synthetic.
os.environ["DJANGO_DEBUG"] = "0"
os.environ["DJANGO_SECRET_KEY"] = CI_SYNTHETIC_SECRET_KEY
os.environ["DJANGO_SQLITE_PATH"] = str(CI_TEMP_ROOT / "db" / "ci.sqlite3")
os.environ["DJANGO_MEDIA_ROOT"] = str(CI_TEMP_ROOT / "media")
os.environ["DJANGO_LOGS_ROOT"] = str(CI_TEMP_ROOT / "logs")

from .settings import *  # noqa: E402,F403

DEBUG = False
SECRET_KEY = CI_SYNTHETIC_SECRET_KEY
ALLOWED_HOSTS = ["localhost", "127.0.0.1", "testserver"]
CSRF_TRUSTED_ORIGINS = ["https://testserver"]

STATIC_ROOT = CI_TEMP_ROOT / "staticfiles"
MEDIA_ROOT = CI_TEMP_ROOT / "media"

EMAIL_BACKEND = "django.core.mail.backends.locmem.EmailBackend"
CELERY_BROKER_URL = "memory://"
CELERY_RESULT_BACKEND = "cache+memory://"
CELERY_TASK_ALWAYS_EAGER = True
CELERY_TASK_EAGER_PROPAGATES = True
CACHES = {
    "default": {
        "BACKEND": "django.core.cache.backends.locmem.LocMemCache",
        "LOCATION": "bmycure4mm-ci-only",
    }
}
PASSWORD_HASHERS = ["django.contrib.auth.hashers.MD5PasswordHasher"]

# CI writes no application log files; the runner retains console output only.
LOGGING = {
    "version": 1,
    "disable_existing_loggers": False,
    "handlers": {"console": {"class": "logging.StreamHandler"}},
    "root": {"handlers": ["console"], "level": "WARNING"},
}
