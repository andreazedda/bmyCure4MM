"""Encrypted email outbox, retries and account expiry housekeeping."""

import json
import time
from datetime import timedelta
from urllib.parse import urlsplit

from django.conf import settings
from django.core.mail import EmailMultiAlternatives
from django.db import transaction
from django.db.models import Q
from django.utils import timezone

from .models import AccessProfile as Identity
from .models import Capability, Notification, RateBucket
from .services import (
    audit,
    cipher,
    locked_policy,
    notify,
    public_base,
)


@transaction.atomic
def dispatch_notification():
    locked_policy()
    # The outbox is committed before dispatch. Database row locks prevent two
    # mail workers from dispatching the same item concurrently.
    item = (
        Notification.objects.select_for_update(skip_locked=True)
        .filter(
            Q(next_attempt_at=None) | Q(next_attempt_at__lte=timezone.now()),
            state__in=["queued", "failed"],
            attempts__lt=5,
        )
        .order_by("created_at")
        .first()
    )
    if not item:
        return False
    capability = None
    if item.deduplication_key.startswith("capability:"):
        capability = Capability.objects.filter(pk=int(item.deduplication_key.split(":", 1)[1])).first()
        if (
            not capability
            or capability.consumed_at
            or capability.revoked_at
            or capability.expires_at <= timezone.now()
        ):
            item.state, item.error_code, item.encrypted_payload = "failed", "CAPABILITY_NO_LONGER_VALID", ""
            item.attempts = 5
            item.save()
            return True
    if item.created_at < timezone.now() - timedelta(hours=24) or not item.encrypted_payload:
        item.state, item.error_code, item.encrypted_payload = "failed", "PAYLOAD_EXPIRED", ""
        item.attempts = 5
        item.save()
        return True
    from django.utils.crypto import salted_hmac

    window = int(time.time()) // 3600
    buckets = []
    for name, ceiling in (
        ("global", settings.ACCOUNT_EMAIL_GLOBAL_HOURLY),
        (str(item.identity_id), settings.ACCOUNT_EMAIL_IDENTITY_HOURLY),
    ):
        key = salted_hmac("bmy-community-mail", name, algorithm="sha256").hexdigest()
        RateBucket.objects.get_or_create(key=key, defaults={"window": window})
        bucket = RateBucket.objects.select_for_update().get(key=key)
        if bucket.window != window:
            bucket.window, bucket.count = window, 0
        if bucket.count >= ceiling:
            item.next_attempt_at = timezone.now() + timedelta(seconds=3600 - int(time.time()) % 3600)
            item.error_code = "MAIL_BUDGET_WAIT"
            item.save(update_fields=["next_attempt_at", "error_code"])
            return True
        buckets.append(bucket)
    for bucket in buckets:
        bucket.count += 1
        bucket.save()
    item.attempts += 1
    try:
        payload = json.loads(cipher().decrypt(item.encrypted_payload.encode(), ttl=86400))
        from .emails import presentation

        subject, plain, html, sender = presentation(item, payload, capability)
        message = EmailMultiAlternatives(
            subject=subject,
            body=plain,
            from_email=sender,
            to=[item.identity.email],
            reply_to=[sender],
            headers={"Message-ID": f"<{item.pk}@{urlsplit(public_base()).hostname}>"},
        )
        message.attach_alternative(html, "text/html")
        if message.send(fail_silently=False) != 1:
            raise RuntimeError("MAIL_NOT_ACCEPTED")
    except Exception:
        item.state, item.error_code = "failed", "MAIL_BACKEND_FAILED"
        item.next_attempt_at = timezone.now() + timedelta(seconds=min(3600, 30 * 2**item.attempts))
    else:
        item.state, item.error_code, item.sent_at = "sent", "", timezone.now()
        item.encrypted_payload = ""
    item.save()
    return True


@transaction.atomic
def housekeeping():
    locked_policy()
    now = timezone.now()
    for identity in Identity.objects.filter(
        is_owner=False,
        request_state=Identity.RequestState.APPROVED,
        access_state__in=[Identity.AccessState.ACTIVE, Identity.AccessState.PENDING],
        expires_at__lte=now,
    ).select_related("user"):
        identity.access_state = Identity.AccessState.EXPIRED
        identity.session_version += 1
        identity.save(update_fields=["access_state", "session_version"])
        if identity.user_id:
            identity.user.is_active = False
            identity.user.save(update_fields=["is_active"])
        Capability.objects.filter(identity=identity, consumed_at=None).update(revoked_at=now)
        audit("access_expired", identity)
    for identity in Identity.objects.filter(
        is_owner=False,
        access_state=Identity.AccessState.ACTIVE,
        expires_at__gt=now,
        expires_at__lte=now + timedelta(days=1),
    ):
        notify(
            identity,
            "expiry",
            f"expiry:{identity.pk}:{identity.expires_at.isoformat()}",
            "Your account access expires within 24 hours. Contact the owner through a ticket if you need continued access.",
        )
    Notification.objects.filter(created_at__lt=now - timedelta(hours=24)).exclude(
        encrypted_payload=""
    ).update(encrypted_payload="")
    Identity.objects.filter(
        is_owner=False,
        created_at__lt=now - timedelta(days=30),
        request_state__in=[Identity.RequestState.UNVERIFIED, Identity.RequestState.REJECTED],
    ).update(reason="", affiliation="", professional_role="")
