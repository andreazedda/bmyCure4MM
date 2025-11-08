from __future__ import annotations

import json
import logging

from django.http import HttpResponseBadRequest, JsonResponse
from django.views.decorators.csrf import csrf_exempt
from django.views.decorators.http import require_POST

logger = logging.getLogger("ux.education")

MAX_BODY_BYTES = 2048


@csrf_exempt
@require_POST
def audit(request):
    body = request.body[:MAX_BODY_BYTES]
    try:
        data = json.loads(body.decode("utf-8"))
    except (json.JSONDecodeError, UnicodeDecodeError):
        return HttpResponseBadRequest("invalid-json")
    event = str(data.get("event", ""))[:64]
    detail = data.get("detail") or {}
    timestamp = data.get("ts")
    remote_addr = request.META.get("REMOTE_ADDR")
    logger.info(
        "ux_event=%s detail=%s ts=%s ip=%s", event, detail, timestamp, remote_addr
    )
    return JsonResponse({"ok": True})
