from __future__ import annotations

import logging
import time
from typing import Callable

from django.http import HttpRequest, HttpResponse

activity_logger = logging.getLogger("activity")
embed_debug_logger = logging.getLogger("embed_debug")


class ActivityLoggingMiddleware:
    """Middleware to record every user interaction with the portal."""

    def __init__(self, get_response: Callable[[HttpRequest], HttpResponse]) -> None:
        self.get_response = get_response

    def __call__(self, request: HttpRequest) -> HttpResponse:
        start = time.perf_counter()
        try:
            response = self.get_response(request)
        except Exception as exc:
            duration_ms = (time.perf_counter() - start) * 1000
            self._log_request(request, status_code=500, duration_ms=duration_ms, error=repr(exc))
            raise
        duration_ms = (time.perf_counter() - start) * 1000
        self._log_request(request, status_code=response.status_code, duration_ms=duration_ms)
        return response

    def _log_request(
        self,
        request: HttpRequest,
        *,
        status_code: int,
        duration_ms: float,
        error: str | None = None,
    ) -> None:
        user = getattr(request, "user", None)
        username = user.get_username() if getattr(user, "is_authenticated", False) else "anonymous"
        activity_logger.info(
            "user=%s method=%s path=%s status=%s duration_ms=%.2f ip=%s%s",
            username,
            request.method,
            request.get_full_path(),
            status_code,
            duration_ms,
            request.META.get("REMOTE_ADDR", "unknown"),
            f" error={error}" if error else "",
        )


class EmbedDebugMiddleware:
    """Targeted logging for iframe embedding issues.

    Writes a compact, greppable record to logs/embed_debug.log for requests that
    are likely involved in plot embedding (e.g., /media/sim_plots/*.html).
    """

    def __init__(self, get_response: Callable[[HttpRequest], HttpResponse]) -> None:
        self.get_response = get_response

    def __call__(self, request: HttpRequest) -> HttpResponse:
        start = time.perf_counter()

        try:
            response = self.get_response(request)
        except Exception as exc:
            duration_ms = (time.perf_counter() - start) * 1000
            if self._should_log(request):
                self._log(request, response=None, duration_ms=duration_ms, error=repr(exc))
            raise

        duration_ms = (time.perf_counter() - start) * 1000
        if self._should_log(request):
            self._log(request, response=response, duration_ms=duration_ms)
        return response

    def _should_log(self, request: HttpRequest) -> bool:
        path = request.path or ""
        if path.startswith("/media/sim_plots/"):
            return True
        if path.startswith("/sim/attempts/") and "/plot/embed" in path:
            return True
        if path.startswith("/patients/"):
            # Helps correlate the page that *contains* the iframe.
            return True
        return False

    def _log(
        self,
        request: HttpRequest,
        *,
        response: HttpResponse | None,
        duration_ms: float,
        error: str | None = None,
    ) -> None:
        try:
            user = getattr(request, "user", None)
            username = user.get_username() if getattr(user, "is_authenticated", False) else "anonymous"

            # Request headers that are useful for iframe debugging.
            origin = request.META.get("HTTP_ORIGIN", "")
            referer = request.META.get("HTTP_REFERER", "")
            sec_fetch_dest = request.META.get("HTTP_SEC_FETCH_DEST", "")
            sec_fetch_site = request.META.get("HTTP_SEC_FETCH_SITE", "")
            sec_fetch_mode = request.META.get("HTTP_SEC_FETCH_MODE", "")
            ua = request.META.get("HTTP_USER_AGENT", "")

            # Response headers that commonly break embeds.
            status = getattr(response, "status_code", 500)
            xfo = response.get("X-Frame-Options", "") if response is not None else ""
            csp = response.get("Content-Security-Policy", "") if response is not None else ""
            ct = response.get("Content-Type", "") if response is not None else ""

            embed_debug_logger.info(
                "user=%s method=%s path=%s status=%s duration_ms=%.2f xfo=%r csp=%r ct=%r origin=%r referer=%r sec_dest=%r sec_site=%r sec_mode=%r ua=%r%s",
                username,
                request.method,
                request.get_full_path(),
                status,
                duration_ms,
                xfo,
                csp,
                ct,
                origin,
                referer,
                sec_fetch_dest,
                sec_fetch_site,
                sec_fetch_mode,
                ua,
                f" error={error}" if error else "",
            )
        except Exception:
            # Never break requests due to logging.
            return

