from __future__ import annotations

import os
import re

from .governance import template_governance_context


def governance(request):
    del request
    context = template_governance_context()
    # Keep local asset URLs aligned with the deployed image, including behind
    # caches that retain the same /static/ path across releases. Never expose
    # an arbitrary environment value in an HTML URL.
    revision = os.environ.get("BMY_SOURCE_SHA", "")
    context["asset_version"] = revision if re.fullmatch(r"[0-9a-f]{40}", revision) else "dev"
    return context
