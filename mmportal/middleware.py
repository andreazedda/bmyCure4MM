from __future__ import annotations

import logging
import time
from typing import Callable

from django.http import HttpRequest, HttpResponse

activity_logger = logging.getLogger("activity")


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

