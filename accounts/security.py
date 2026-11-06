import base64
import io
import time

import pyotp
import qrcode
from django.db import transaction

from .models import AccessProfile
from .services import audit, cipher, locked_policy


def mfa_current(request):
    profile = request.access_profile
    return bool(
        profile
        and profile.totp_confirmed
        and request.session.get("owner_mfa_version") == profile.session_version
        and request.session.get("owner_mfa_at", 0) > time.time() - 900
    )


def owner_setup_qr(profile, secret):
    """Render locally; authenticator secrets never leave this authenticated response."""
    uri = pyotp.TOTP(secret).provisioning_uri(name=profile.email, issuer_name="bmyCure4MM")
    output = io.BytesIO()
    qrcode.make(uri, box_size=6, border=4).save(output, format="PNG")
    return "data:image/png;base64," + base64.b64encode(output.getvalue()).decode("ascii")


@transaction.atomic
def owner_secret(profile):
    locked_policy()
    current = AccessProfile.objects.select_for_update().get(pk=profile.pk)
    if not current.totp_secret:
        current.totp_secret = cipher().encrypt(pyotp.random_base32().encode()).decode()
        current.save(update_fields=["totp_secret"])
    return cipher().decrypt(current.totp_secret.encode()).decode()


@transaction.atomic
def confirm_owner_code(profile, code):
    locked_policy()
    current = AccessProfile.objects.select_for_update().get(pk=profile.pk)
    if not current.totp_secret:
        return False
    totp = pyotp.TOTP(cipher().decrypt(current.totp_secret.encode()).decode())
    now_step = int(time.time()) // 30
    for offset in (0, -1, 1):
        step = now_step + offset
        if step > current.totp_last_step and totp.verify(code, for_time=step * 30):
            current.totp_confirmed, current.totp_last_step = True, step
            current.save(update_fields=["totp_confirmed", "totp_last_step"])
            audit("owner_mfa_verified", current, current.user)
            return True
    return False
