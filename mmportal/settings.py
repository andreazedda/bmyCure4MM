"""Django settings for mmportal project."""

from __future__ import annotations

import os
import sys
from pathlib import Path
from django.core.exceptions import ImproperlyConfigured
from django.contrib.messages import constants as django_messages

# Build paths inside the project like this: BASE_DIR / "subdir".
BASE_DIR = Path(__file__).resolve().parent.parent
LOGS_DIR = BASE_DIR / "logs"
LOGS_DIR.mkdir(exist_ok=True)

# SECURITY WARNING: don't run with debug turned on in production!
DEBUG = os.environ.get("DJANGO_DEBUG", "1") == "1"

# SECURITY WARNING: keep the secret key used in production secret!
# In production (DEBUG=False), SECRET_KEY MUST be set via DJANGO_SECRET_KEY environment variable
# In development (DEBUG=True), a default key is provided for convenience
SECRET_KEY = os.environ.get("DJANGO_SECRET_KEY")

if not DEBUG:  # Production mode - strict validation
    if not SECRET_KEY:
        raise ImproperlyConfigured(
            "DJANGO_SECRET_KEY environment variable is required in production. "
            "Generate a secure key with: "
            "python -c 'from django.core.management.utils import get_random_secret_key; print(get_random_secret_key())'"
        )
    
    if "insecure" in SECRET_KEY.lower() or SECRET_KEY in (
        "your-secret-key-here-change-me-in-production",
        "django-insecure-please-change-me",
        "change-me",
        "changeme",
        "dev-key-not-for-production",
    ):
        raise ImproperlyConfigured(
            "DJANGO_SECRET_KEY contains an insecure or default value. "
            "Please set a strong, random secret key in your environment. "
            "Generate one with: "
            "python -c 'from django.core.management.utils import get_random_secret_key; print(get_random_secret_key())'"
        )
else:  # Development mode - allow default for convenience
    if not SECRET_KEY:
        SECRET_KEY = "dev-key-not-for-production"
        print("⚠️  Using default SECRET_KEY for development. Set DJANGO_SECRET_KEY in production!")

_hosts = os.environ.get("ALLOWED_HOSTS")
ALLOWED_HOSTS: list[str] = [host.strip() for host in _hosts.split(",") if host.strip()] if _hosts else []

_trusted = os.environ.get("CSRF_TRUSTED_ORIGINS")
CSRF_TRUSTED_ORIGINS = [origin.strip() for origin in _trusted.split(",") if origin.strip()] if _trusted else []

if DEBUG:
    for host in ("localhost", "127.0.0.1", "testserver"):
        if host not in ALLOWED_HOSTS:
            ALLOWED_HOSTS.append(host)

PREDLAB_V2 = os.environ.get("PREDLAB_V2", "0") == "1"


# Application definition

INSTALLED_APPS = [
    "django.contrib.admin",
    "django.contrib.auth",
    "django.contrib.contenttypes",
    "django.contrib.sessions",
    "django.contrib.messages",
    "django.contrib.staticfiles",
    "django_filters",
    "clinic",
    "chemtools",
    "simulator",
    "docs_viewer",
    "rest_framework",
]

MIDDLEWARE = [
    "django.middleware.security.SecurityMiddleware",
    "django.contrib.sessions.middleware.SessionMiddleware",
    "django.middleware.common.CommonMiddleware",
    "django.middleware.csrf.CsrfViewMiddleware",
    "django.contrib.auth.middleware.AuthenticationMiddleware",
    "mmportal.middleware.EmbedDebugMiddleware",
    "mmportal.middleware.ActivityLoggingMiddleware",
    "django.contrib.messages.middleware.MessageMiddleware",
    "django.middleware.clickjacking.XFrameOptionsMiddleware",
]

# The Clinic embeds locally-generated Plotly HTML reports (e.g., /media/sim_plots/*.html)
# inside iframes on same-origin pages.
X_FRAME_OPTIONS = "SAMEORIGIN"

ROOT_URLCONF = "mmportal.urls"

TEMPLATES = [
    {
        "BACKEND": "django.template.backends.django.DjangoTemplates",
        "DIRS": [BASE_DIR / "templates"],
        "APP_DIRS": True,
        "OPTIONS": {
            "context_processors": [
                "django.template.context_processors.debug",
                "django.template.context_processors.request",
                "django.contrib.auth.context_processors.auth",
                "django.contrib.messages.context_processors.messages",
            ],
        },
    },
]

WSGI_APPLICATION = "mmportal.wsgi.application"


# Database
# https://docs.djangoproject.com/en/5.0/ref/settings/#databases

DATABASES = {
    "default": {
        "ENGINE": "django.db.backends.sqlite3",
        "NAME": BASE_DIR / "db.sqlite3",
    }
}


# Password validation
# https://docs.djangoproject.com/en/5.0/ref/settings/#auth-password-validators

AUTH_PASSWORD_VALIDATORS = [
    {
        "NAME": "django.contrib.auth.password_validation.UserAttributeSimilarityValidator",
    },
    {
        "NAME": "django.contrib.auth.password_validation.MinimumLengthValidator",
    },
    {
        "NAME": "django.contrib.auth.password_validation.CommonPasswordValidator",
    },
    {
        "NAME": "django.contrib.auth.password_validation.NumericPasswordValidator",
    },
]


# Internationalization
# https://docs.djangoproject.com/en/5.0/topics/i18n/

LANGUAGE_CODE = "en-us"

TIME_ZONE = "UTC"

USE_I18N = True

USE_TZ = True


# Static files (CSS, JavaScript, Images)
# https://docs.djangoproject.com/en/5.0/howto/static-files/

STATIC_URL = "/static/"
STATIC_ROOT = BASE_DIR / "staticfiles"

STATICFILES_DIRS: list[Path] = [BASE_DIR / "mmportal" / "static"]

MEDIA_URL = "/media/"
MEDIA_ROOT = BASE_DIR / "media"

LOGIN_URL = "/admin/login/"
LOGIN_REDIRECT_URL = "/"

DEFAULT_AUTO_FIELD = "django.db.models.BigAutoField"

LOGGING = {
    "version": 1,
    "disable_existing_loggers": False,
    "formatters": {
        "activity": {
            "format": "%(asctime)s [%(levelname)s] %(message)s",
        },
        "verbose": {
            "format": "{levelname} {asctime} {module} {message}",
            "style": "{",
        },
        "simple": {
            "format": "{levelname} {message}",
            "style": "{",
        },
    },
    "handlers": {
        "activity_file": {
            "level": "INFO",
            "class": "logging.FileHandler",
            "filename": str(LOGS_DIR / "activity.log"),
            "formatter": "activity",
        },
        "embed_debug_file": {
            "level": "INFO",
            "class": "logging.FileHandler",
            "filename": str(LOGS_DIR / "embed_debug.log"),
            "formatter": "activity",
        },
        "console": {
            "class": "logging.StreamHandler",
            "formatter": "verbose",
        },
        "file": {
            "class": "logging.FileHandler",
            "filename": os.path.join(BASE_DIR, "logs", "django.log"),
            "formatter": "verbose",
        },
        "celery_file": {
            "class": "logging.FileHandler",
            "filename": os.path.join(BASE_DIR, "logs", "celery_tasks.log"),
            "formatter": "verbose",
        },
    },
    "loggers": {
        "activity": {
            "handlers": ["activity_file"],
            "level": "INFO",
            "propagate": False,
        },
        "embed_debug": {
            "handlers": ["embed_debug_file"],
            "level": "INFO",
            "propagate": False,
        },
        "ux.education": {
            "handlers": ["activity_file"],
            "level": "INFO",
            "propagate": False,
        },
        "django": {
            "handlers": ["console", "file"],
            "level": "INFO",
        },
        "chemtools.tasks": {
            "handlers": ["console", "celery_file"],
            "level": "INFO",
            "propagate": False,
        },
        "celery": {
            "handlers": ["console", "celery_file"],
            "level": "INFO",
            "propagate": False,
        },
    },
}

MESSAGE_STORAGE = "django.contrib.messages.storage.session.SessionStorage"
MESSAGE_TAGS = {
    django_messages.ERROR: "danger",
}

REST_FRAMEWORK = {
    "DEFAULT_PERMISSION_CLASSES": [
        "rest_framework.permissions.IsAuthenticated",
    ],
    "DEFAULT_AUTHENTICATION_CLASSES": [
        "rest_framework.authentication.SessionAuthentication",
        "rest_framework.authentication.BasicAuthentication",
    ],
}

CELERY_BROKER_URL = os.environ.get("CELERY_BROKER_URL", "redis://localhost:6379/0")
CELERY_RESULT_BACKEND = os.environ.get("CELERY_RESULT_BACKEND", CELERY_BROKER_URL)

# For development: run tasks synchronously without Redis/Celery worker
# Set CELERY_TASK_ALWAYS_EAGER=1 in environment or uncomment below:
# CELERY_TASK_ALWAYS_EAGER = True
# CELERY_TASK_EAGER_PROPAGATES = True


