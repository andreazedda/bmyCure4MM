"""Disposable synthetic verification; never opens the workspace database or media."""

import os
import tempfile
from pathlib import Path

from cryptography.fernet import Fernet

TEST_ROOT = Path(os.environ.get("BMY_TEST_ROOT", tempfile.mkdtemp(prefix="bmy-test-")))
os.environ["DJANGO_SQLITE_PATH"] = os.environ.get("BMY_TEST_DATABASE_NAME", str(TEST_ROOT / "db.sqlite3"))
os.environ["DJANGO_MEDIA_ROOT"] = str(TEST_ROOT / "media")
os.environ["DJANGO_LOGS_ROOT"] = str(TEST_ROOT / "logs")
os.environ.setdefault("ACCOUNT_TOKEN_ENCRYPTION_KEY", Fernet.generate_key().decode())

from .settings import *  # noqa: E402,F403

DATABASES["default"]["TEST"] = {"NAME": str(TEST_ROOT / "test.sqlite3")}
STATIC_ROOT = TEST_ROOT / "static"
LOGGING = {
    "version": 1,
    "disable_existing_loggers": True,
    "handlers": {"null": {"class": "logging.NullHandler"}},
    "root": {"handlers": ["null"]},
}
EMAIL_BACKEND = "django.core.mail.backends.locmem.EmailBackend"
PASSWORD_HASHERS = ["django.contrib.auth.hashers.MD5PasswordHasher"]
CELERY_TASK_ALWAYS_EAGER = True
CELERY_TASK_EAGER_PROPAGATES = True
ACCOUNT_EMAIL_TASKS_ENABLED = False
PUBLIC_BASE_URL = "http://127.0.0.1:8765"
ACCOUNT_TOKEN_ENCRYPTION_KEY = os.environ["ACCOUNT_TOKEN_ENCRYPTION_KEY"]
DEFAULT_FROM_EMAIL = "service@example.test"
