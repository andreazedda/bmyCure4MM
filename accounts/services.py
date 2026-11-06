"""Registration lifecycle, access control and encrypted notification outbox."""

import hashlib
import json
import secrets
import time
import uuid
from datetime import timedelta
from urllib.parse import urlsplit

from cryptography.fernet import Fernet
from django.conf import settings
from django.contrib.auth import get_user_model
from django.core.exceptions import PermissionDenied, ValidationError
from django.core.validators import validate_email
from django.db import transaction
from django.db.models import Q
from django.utils import timezone
from django.utils.crypto import salted_hmac

from .models import (
    AccessPolicy,
    AuditEvent,
    Capability,
    Notification,
    RateBucket,
    SupportTicket,
)
from .models import AccessProfile as Identity


def normalize_email(value):
    value = str(value).strip().lower()
    try:
        value.encode("ascii")
    except UnicodeEncodeError as exc:
        raise ValidationError("Use an ASCII email address.") from exc
    validate_email(value)
    if len(value) > 254:
        raise ValidationError("Email address is too long.")
    return value


def locked_policy():
    # Atomic no-op UPDATE also acquires the write lock on the deployed SQLite backend.
    AccessPolicy.objects.filter(pk=1).update(id=1)
    return AccessPolicy.objects.select_for_update().get(pk=1)


def audit(action, identity=None, actor=None, code=""):
    AuditEvent.objects.create(action=action, identity=identity, actor=actor, code=code)


def usable(identity, *, pending=False):
    return bool(
        identity
        and identity.user_id
        and identity.verified_at
        and identity.user.email.strip().lower() == identity.email
        and identity.request_state == Identity.RequestState.APPROVED
        and identity.access_state
        in (
            {Identity.AccessState.ACTIVE, Identity.AccessState.PENDING}
            if pending
            else {Identity.AccessState.ACTIVE}
        )
        and (identity.is_owner or (identity.expires_at and identity.expires_at > timezone.now()))
        and (identity.is_owner or (not identity.user.is_staff and not identity.user.is_superuser))
    )


def require_manager(actor):
    identity = Identity.objects.select_related("user").filter(user=actor).first()
    if not usable(identity) or not identity.user.is_active or not (identity.is_owner or identity.is_approver):
        raise PermissionDenied
    return identity


def require_community(identity):
    if not usable(identity) or identity.is_owner or identity.is_approver:
        raise PermissionDenied
    if not identity.user.is_active:
        raise PermissionDenied


def public_base():
    base = settings.PUBLIC_BASE_URL.rstrip("/")
    parts = urlsplit(base)
    if parts.username or parts.password or parts.query or parts.fragment or parts.path not in {"", "/"}:
        raise ValidationError("PUBLIC_BASE_URL must be a configured origin.")
    if parts.scheme != "https" and not (settings.DEBUG and parts.hostname in {"localhost", "127.0.0.1"}):
        raise ValidationError("PUBLIC_BASE_URL requires HTTPS.")
    if not parts.hostname:
        raise ValidationError("PUBLIC_BASE_URL is required.")
    return base


def cipher():
    return Fernet(settings.ACCOUNT_TOKEN_ENCRYPTION_KEY.encode("ascii"))


def notify(identity, kind, key, body):
    public_base()
    encrypted = cipher().encrypt(json.dumps({"body": body}).encode()).decode()
    notification, created = Notification.objects.get_or_create(
        deduplication_key=key,
        defaults={"identity": identity, "kind": kind, "encrypted_payload": encrypted},
    )
    if created:
        from .tasks import request_delivery

        transaction.on_commit(request_delivery)
    return notification


def issue_capability(identity, purpose, *, revoke=True):
    if revoke:
        previous = Capability.objects.filter(
            identity=identity, purpose=purpose, consumed_at=None, revoked_at=None
        )
        keys = [f"capability:{pk}" for pk in previous.values_list("pk", flat=True)]
        previous.update(revoked_at=timezone.now())
        Notification.objects.filter(deduplication_key__in=keys).exclude(state="sent").update(
            encrypted_payload="", state="failed", error_code="CAPABILITY_REVOKED", attempts=5
        )
    raw = secrets.token_urlsafe(32)
    capability = Capability.objects.create(
        identity=identity,
        purpose=purpose,
        digest=hashlib.sha256(raw.encode()).hexdigest(),
        expires_at=timezone.now() + timedelta(hours=24 if purpose == "verify" else 2),
    )
    paths = {
        "verify": "/access/verify/",
        "activate": "/access/activate/",
        "owner": "/owner/activate/",
        "reset": "/password/reset/confirm/",
    }
    # Fragments are not transmitted in HTTP requests or proxy access logs.
    link = f"{public_base()}{paths[purpose]}#token={raw}"
    notify(
        identity,
        purpose,
        f"capability:{capability.pk}",
        f"Confirm your {purpose} action and, when requested, set your password:\n{link}\nThis link expires and can only be used once. Do not forward it.",
    )
    return capability


@transaction.atomic
def request_access(data):
    policy = locked_policy()
    if not policy.requests_open:
        raise ValidationError("Access requests are temporarily paused.")
    email = normalize_email(data["email"])
    # Public submission cannot claim the configured owner or an existing identity.
    configured_owner = getattr(settings, "PLATFORM_OWNER_EMAIL", "").strip().lower()
    if email == configured_owner or get_user_model().objects.filter(email__iexact=email).exists():
        return None
    existing = Identity.objects.filter(email=email).first()
    if existing:
        # A forgotten/expired verification must not reserve an email forever.
        # Preserve valid links and all original request fields. Public retries
        # never grant a status-session capability merely from knowing an email.
        if (
            existing.request_state == Identity.RequestState.UNVERIFIED
            and existing.access_state == Identity.AccessState.PENDING
            and existing.verified_at is None
            and existing.user_id is None
            and not existing.is_owner
            and not existing.is_approver
            and not Capability.objects.filter(
                identity=existing,
                purpose="verify",
                consumed_at=None,
                revoked_at=None,
                expires_at__gt=timezone.now(),
            ).exists()
        ):
            issue_capability(existing, "verify")
            audit("verification_reissued", existing)
        return None
    identity = Identity.objects.create(
        email=email,
        nickname=data["nickname"],
        reason=data["reason"],
        affiliation=data.get("affiliation", ""),
        professional_role=data.get("professional_role", ""),
        acknowledged_at=timezone.now(),
    )
    issue_capability(identity, "verify")
    audit("request_created", identity)
    return identity


@transaction.atomic
def bootstrap_owner():
    locked_policy()
    email = normalize_email(settings.PLATFORM_OWNER_EMAIL)
    username = settings.PLATFORM_OWNER_USERNAME.strip()
    name = settings.PLATFORM_OWNER_DISPLAY_NAME.strip()
    if not username or not name:
        raise ValidationError("Owner configuration is incomplete.")
    existing = Identity.objects.filter(is_owner=True).select_related("user").first()
    if existing:
        if (
            existing.email != email
            or existing.user.username != username
            or existing.user.email.strip().lower() != email
        ):
            raise ValidationError("Owner identity conflict; use the documented recovery procedure.")
        return existing, False
    if Identity.objects.filter(email=email).exists():
        raise ValidationError("Owner identity conflict; no existing account was promoted.")
    matches = get_user_model().objects.filter(Q(username=username) | Q(email__iexact=email))
    if matches.exists():
        if matches.count() != 1:
            raise ValidationError("Owner identity conflict.")
        user = matches.get()
        if (
            user.username != username
            or user.email.strip().lower() != email
            or not user.is_active
            or not (user.is_staff and user.is_superuser)
        ):
            raise ValidationError("Only an exact, already privileged owner identity may be linked.")
        identity = Identity.objects.create(
            email=email,
            nickname=name[:60],
            user=user,
            is_owner=True,
            request_state=Identity.RequestState.APPROVED,
            access_state=Identity.AccessState.ACTIVE,
            verified_at=timezone.now(),
        )
        audit("existing_owner_linked", identity)
        return identity, True
    user = get_user_model()(
        username=username, email=email, first_name=name, is_active=False, is_staff=False, is_superuser=False
    )
    user.set_unusable_password()
    user.save()
    identity = Identity.objects.create(
        email=email,
        nickname=name[:60],
        user=user,
        is_owner=True,
        request_state=Identity.RequestState.APPROVED,
        verified_at=timezone.now(),
    )
    issue_capability(identity, "owner")
    audit("owner_bootstrapped", identity)
    return identity, True


@transaction.atomic
def confirm_capability(raw, purpose, password=None):
    policy = locked_policy()
    digest = hashlib.sha256(str(raw).encode()).hexdigest()
    token = (
        Capability.objects.select_for_update(of=("self",))
        .select_related("identity__user")
        .filter(digest=digest, purpose=purpose)
        .first()
    )
    if not token or token.consumed_at or token.revoked_at or token.expires_at <= timezone.now():
        raise ValidationError("This link is invalid, expired or already used.")
    identity = token.identity
    if purpose == "verify":
        if identity.request_state != Identity.RequestState.UNVERIFIED:
            raise ValidationError("This request cannot be verified.")
        identity.verified_at = timezone.now()
        identity.request_state = Identity.RequestState.REVIEW
        owner = Identity.objects.filter(is_owner=True).first()
        if owner:
            notify(
                owner,
                "verified_request",
                f"verified:{identity.pk}",
                "A verified access request is waiting in /manage/access-requests/.",
            )
    else:
        if not usable(identity, pending=True):
            raise ValidationError("Approval or access is no longer valid.")
        if purpose == "owner" and not identity.is_owner or purpose == "activate" and identity.is_owner:
            raise ValidationError("This link has a different purpose.")
        if purpose in {"activate", "owner"} and identity.access_state != Identity.AccessState.PENDING:
            raise ValidationError("This account has already been activated.")
        if purpose == "reset" and identity.access_state != Identity.AccessState.ACTIVE:
            raise ValidationError("Password reset cannot activate an account.")
        if not identity.is_owner:
            if not policy.activations_open:
                raise ValidationError("Community activations are temporarily paused.")
        from django.contrib.auth.password_validation import validate_password

        validate_password(password, identity.user)
        identity.user.set_password(password)
        identity.user.is_active = True
        identity.user.save(update_fields=["password", "is_active"])
        identity.access_state = Identity.AccessState.ACTIVE
        identity.session_version += 1
    identity.save()
    token.consumed_at = timezone.now()
    token.save(update_fields=["consumed_at"])
    audit(f"capability_{purpose}_confirmed", identity)
    return identity


@transaction.atomic
def decide_request(actor, identity_id, *, approve, days=14):
    policy = locked_policy()
    require_manager(actor)
    identity = Identity.objects.select_for_update().get(pk=identity_id)
    if identity.request_state in {Identity.RequestState.APPROVED, Identity.RequestState.REJECTED}:
        return identity
    if identity.request_state != Identity.RequestState.REVIEW or not identity.verified_at:
        raise ValidationError("Only a verified request can be reviewed.")
    if not approve:
        identity.request_state = Identity.RequestState.REJECTED
        identity.save(update_fields=["request_state"])
        notify(
            identity,
            "rejected",
            f"rejected:{identity.pk}",
            "Your community request was not approved. Contact the platform owner if you need clarification.",
        )
        audit("request_rejected", identity, actor)
        return identity
    if not policy.activations_open:
        raise ValidationError("Account activations are temporarily paused.")
    if not 1 <= days <= 14:
        raise ValidationError("Initial access is limited to 14 days.")
    if get_user_model().objects.filter(email__iexact=identity.email).exists():
        raise ValidationError("Email identity conflict; no accounts were merged.")
    user = get_user_model()(username=f"community-{uuid.uuid4().hex}", email=identity.email, is_active=False)
    user.set_unusable_password()
    user.save()
    identity.user = user
    identity.expires_at = timezone.now() + timedelta(days=days)
    identity.request_state = Identity.RequestState.APPROVED
    identity.save()
    issue_capability(identity, "activate")
    audit("request_approved", identity, actor)
    return identity


@transaction.atomic
def resend_activation(actor, identity_id):
    locked_policy()
    require_manager(actor)
    identity = Identity.objects.select_related("user").get(pk=identity_id)
    if not usable(identity, pending=True) or identity.access_state != Identity.AccessState.PENDING:
        raise ValidationError("Only an approved, unexpired pending account can receive an activation.")
    issue_capability(identity, "owner" if identity.is_owner else "activate")
    audit("activation_reissued", identity, actor)


@transaction.atomic
def change_access(actor, identity_id, state):
    locked_policy()
    require_manager(actor)
    identity = Identity.objects.select_for_update().get(pk=identity_id)
    if identity.is_owner or state not in {Identity.AccessState.SUSPENDED, Identity.AccessState.REVOKED}:
        raise ValidationError("This access transition is not permitted.")
    if identity.access_state == state:
        return
    identity.access_state = state
    identity.session_version += 1
    identity.save(update_fields=["access_state", "session_version"])
    if identity.user_id:
        identity.user.is_active = False
        identity.user.save(update_fields=["is_active"])
    Capability.objects.filter(identity=identity, consumed_at=None).update(revoked_at=timezone.now())
    notify(
        identity,
        state.lower(),
        f"state:{identity.pk}:{identity.session_version}",
        f"Your community access is now {state.lower()}.",
    )
    audit("access_changed", identity, actor, state)


@transaction.atomic
def renew_access(actor, identity_id, days=14):
    policy = locked_policy()
    require_manager(actor)
    identity = Identity.objects.select_related("user").get(pk=identity_id)
    if (
        not policy.activations_open
        or not 1 <= days <= 14
        or identity.is_owner
        or identity.is_approver
        or identity.request_state != Identity.RequestState.APPROVED
        or identity.access_state not in {"ACTIVE", "SUSPENDED", "EXPIRED"}
        or not identity.user_id
        or identity.user.is_staff
        or identity.user.is_superuser
        or identity.user.email.strip().lower() != identity.email
        or not identity.user.has_usable_password()
    ):
        raise ValidationError(
            "Only an approved activated account can be renewed; revoked access stays revoked."
        )
    identity.expires_at = timezone.now() + timedelta(days=days)
    identity.access_state = Identity.AccessState.ACTIVE
    identity.session_version += 1
    identity.save(update_fields=["expires_at", "access_state", "session_version"])
    identity.user.is_active = True
    identity.user.save(update_fields=["is_active"])
    audit("access_renewed", identity, actor)
    notify(
        identity,
        "renewed",
        f"renewed:{identity.pk}:{identity.session_version}",
        "Your access was renewed by the platform owner. Sign in again to continue.",
    )
    return identity


@transaction.atomic
def rate_limit(scope, request, identity_key="", limit=10, seconds=3600):
    locked_policy()
    now_window = int(time.time()) // seconds
    ip = request.META.get("REMOTE_ADDR", "unknown")
    keys = [(f"{scope}:global", limit * 20), (f"{scope}:ip:{ip}", limit)]
    if identity_key:
        keys.append((f"{scope}:identity:{identity_key}", limit))
    buckets = []
    for raw, maximum in keys:
        key = salted_hmac("bmy-community-rate", raw, algorithm="sha256").hexdigest()
        bucket, _ = RateBucket.objects.get_or_create(key=key, defaults={"window": now_window})
        if bucket.window != now_window:
            bucket.window, bucket.count = now_window, 0
        if bucket.count >= maximum:
            return False
        buckets.append(bucket)
    for bucket in buckets:
        bucket.count += 1
        bucket.save()
    return True


@transaction.atomic
def submit_feedback(identity, data, key):
    locked_policy()
    identity = Identity.objects.select_related("user").get(pk=identity.pk)
    require_community(identity)
    from mmportal.governance import APP_VERSION

    ticket, created = SupportTicket.objects.get_or_create(
        submission_key=key,
        defaults={
            "identity": identity,
            "category": data["category"],
            "message": data["message"],
            "app_version": APP_VERSION,
        },
    )
    if ticket.identity_id != identity.pk:
        raise PermissionDenied
    if created:
        audit("feedback_created", identity)
        notify(
            identity,
            "feedback_ack",
            f"feedback:{ticket.pk}",
            "Your feedback is saved. You can follow its status at /feedback/.",
        )
        owner = Identity.objects.filter(is_owner=True).first()
        if owner:
            notify(
                owner,
                "feedback",
                f"feedback-owner:{ticket.pk}",
                "New community feedback is waiting at /manage/feedback/.",
            )
    return ticket
