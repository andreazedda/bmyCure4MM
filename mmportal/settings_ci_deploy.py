"""Synthetic deployment-check settings used only by CI.

This module exercises Django deployment checks without importing production
credentials or claiming that the repository has a complete production runtime.
"""

from .settings import *  # noqa: F403

DEBUG = False
ALLOWED_HOSTS = ["ci.invalid"]
CSRF_TRUSTED_ORIGINS = ["https://ci.invalid"]

SECURE_SSL_REDIRECT = True
SESSION_COOKIE_SECURE = True
CSRF_COOKIE_SECURE = True
SECURE_HSTS_SECONDS = 31536000
SECURE_HSTS_INCLUDE_SUBDOMAINS = True
SECURE_HSTS_PRELOAD = True
SECURE_CONTENT_TYPE_NOSNIFF = True
SECURE_REFERRER_POLICY = "same-origin"

# SAMEORIGIN is deliberate in the current product because same-origin Plotly
# report artifacts are embedded. The broader CSP/frame policy remains tracked
# under GitHub issue #9; this CI module does not redefine that product decision.
X_FRAME_OPTIONS = "SAMEORIGIN"
