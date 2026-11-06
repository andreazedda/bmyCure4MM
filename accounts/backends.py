from django.contrib.auth.backends import ModelBackend
from django.contrib.auth.hashers import check_password, make_password

from .models import AccessProfile
from .services import normalize_email, usable


class VerifiedEmailBackend(ModelBackend):
    def authenticate(self, request, email=None, password=None, **kwargs):
        if email is None or password is None:
            return None
        try:
            email = normalize_email(email)
        except Exception:
            return None
        identity = AccessProfile.objects.select_related("user").filter(email=email).first()
        if not identity or not identity.user_id:
            # Match the expensive password operation for unknown identities.
            make_password(password)
            return None
        valid = check_password(password, identity.user.password)
        if valid and usable(identity) and identity.user.is_active:
            return identity.user
        return None

    def get_user(self, user_id):
        user = super().get_user(user_id)
        if user is None:
            return None
        identity = AccessProfile.objects.select_related("user").filter(user=user).first()
        return user if usable(identity) else None
