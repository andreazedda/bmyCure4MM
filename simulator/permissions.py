from __future__ import annotations

from django.contrib.auth.models import User


def is_editor(user: User) -> bool:
    """Return True when user can manage simulator content."""
    if not user.is_authenticated:
        return False
    if user.is_staff:
        return True
    return user.groups.filter(name="Simulator Editors").exists()
