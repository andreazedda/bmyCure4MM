from __future__ import annotations

from .governance import template_governance_context


def governance(request):
    del request
    return template_governance_context()
