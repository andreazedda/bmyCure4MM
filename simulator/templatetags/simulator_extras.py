from __future__ import annotations

from django import template

from simulator.permissions import is_editor

register = template.Library()


@register.simple_tag
def simulator_is_editor(user) -> bool:
    """Return True when the given user can manage simulator content."""
    return is_editor(user)


@register.filter
def dict_get(mapping, key):
    """Safe dictionary lookup in templates."""
    if isinstance(mapping, dict):
        return mapping.get(key, "")
    return ""
