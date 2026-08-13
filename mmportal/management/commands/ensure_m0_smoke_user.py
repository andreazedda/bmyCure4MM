from __future__ import annotations

import stat
import sys
from pathlib import Path

from django.contrib.auth import get_user_model
from django.core.management.base import BaseCommand, CommandError


SMOKE_PREFIX = "m0-smoke-"
SMOKE_EMAIL_DOMAIN = "@m0-smoke.invalid"


def _read_password(*, from_stdin: bool, password_file: str) -> str:
    if from_stdin:
        value = sys.stdin.readline().rstrip("\r\n")
    else:
        path = Path(password_file)
        try:
            mode = stat.S_IMODE(path.stat().st_mode)
        except OSError as exc:
            raise CommandError("Smoke password file is not readable.") from exc
        if mode & 0o077:
            raise CommandError("Smoke password file must not be accessible by group or others.")
        try:
            value = path.read_text(encoding="utf-8").splitlines()[0]
        except (OSError, IndexError) as exc:
            raise CommandError("Smoke password file is empty or unreadable.") from exc
    if not value:
        raise CommandError("Smoke password must not be empty.")
    return value


class Command(BaseCommand):
    help = "Create, rotate, or disable a least-privilege temporary M0 smoke user."

    def add_arguments(self, parser):
        parser.add_argument("--username", required=True)
        action = parser.add_mutually_exclusive_group(required=True)
        action.add_argument("--password-stdin", action="store_true")
        action.add_argument("--password-file", default="")
        action.add_argument("--disable", action="store_true")

    def handle(self, *args, **options):
        username = str(options["username"]).strip()
        if not username.startswith(SMOKE_PREFIX):
            raise CommandError(f"Username must use the reserved {SMOKE_PREFIX!r} prefix.")

        user_model = get_user_model()
        user = user_model.objects.filter(username=username).first()
        marker_email = f"{username}{SMOKE_EMAIL_DOMAIN}"

        if user is not None and str(getattr(user, "email", "") or "") != marker_email:
            raise CommandError("Refusing to modify an existing account without the M0 smoke marker.")

        created = user is None
        if created:
            user = user_model(username=username)
            if hasattr(user, "email"):
                user.email = marker_email

        user.is_staff = False
        user.is_superuser = False
        if options["disable"]:
            user.is_active = False
            user.set_unusable_password()
            action = "disabled"
        else:
            password = _read_password(
                from_stdin=bool(options["password_stdin"]),
                password_file=str(options["password_file"] or ""),
            )
            user.is_active = True
            user.set_password(password)
            action = "created" if created else "rotated"

        user.save()
        if hasattr(user, "groups"):
            user.groups.clear()
        if hasattr(user, "user_permissions"):
            user.user_permissions.clear()

        self.stdout.write(
            "M0_SMOKE_USER "
            f"username={username} action={action} active={str(bool(user.is_active)).lower()} "
            "staff=false superuser=false groups=0 direct_permissions=0"
        )
