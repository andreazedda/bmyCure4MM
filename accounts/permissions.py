from rest_framework.permissions import BasePermission

from .models import AccessProfile


class ExistingDataAccess(BasePermission):
    """Registration alone never grants legacy clinical API or Basic-auth access."""

    def has_permission(self, request, view):
        if not request.user or not request.user.is_authenticated:
            return False
        profile = AccessProfile.objects.filter(user=request.user).first()
        return profile is None or (profile.is_owner and request.user.is_staff)
