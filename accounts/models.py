"""Registration and ticket metadata extending the existing Django user."""

import uuid

from django.conf import settings
from django.db import models
from django.db.models import Q


class AccessPolicy(models.Model):
    # This seeded singleton is the common lock for admissions and lifecycle changes.
    id = models.PositiveSmallIntegerField(primary_key=True, default=1, editable=False)
    requests_open = models.BooleanField(default=False)
    activations_open = models.BooleanField(default=False)

    class Meta:
        constraints = [models.CheckConstraint(condition=Q(pk=1), name="account_policy_singleton")]


class AccessProfile(models.Model):
    class RequestState(models.TextChoices):
        UNVERIFIED = "PENDING_EMAIL_VERIFICATION", "Verify your email"
        REVIEW = "PENDING_REVIEW", "Waiting for review"
        APPROVED = "APPROVED", "Approved"
        REJECTED = "REJECTED", "Rejected"

    class AccessState(models.TextChoices):
        PENDING = "PENDING_ACTIVATION", "Set your password"
        ACTIVE = "ACTIVE", "Active"
        SUSPENDED = "SUSPENDED", "Suspended"
        EXPIRED = "EXPIRED", "Expired"
        REVOKED = "REVOKED", "Revoked"

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    email = models.EmailField(unique=True)
    nickname = models.CharField(max_length=60)
    user = models.OneToOneField(
        settings.AUTH_USER_MODEL, null=True, on_delete=models.PROTECT, related_name="community_identity"
    )
    request_state = models.CharField(
        max_length=32, choices=RequestState.choices, default=RequestState.UNVERIFIED
    )
    access_state = models.CharField(max_length=24, choices=AccessState.choices, default=AccessState.PENDING)
    is_owner = models.BooleanField(default=False)
    is_approver = models.BooleanField(default=False)
    reason = models.CharField(max_length=1000, blank=True)
    affiliation = models.CharField(max_length=120, blank=True)
    professional_role = models.CharField(max_length=60, blank=True)
    acknowledged_at = models.DateTimeField(null=True)
    verified_at = models.DateTimeField(null=True)
    created_at = models.DateTimeField(auto_now_add=True)
    expires_at = models.DateTimeField(null=True)
    session_version = models.PositiveIntegerField(default=1)
    totp_secret = models.TextField(blank=True)
    totp_confirmed = models.BooleanField(default=False)
    totp_last_step = models.BigIntegerField(default=-1)

    class Meta:
        constraints = [
            models.UniqueConstraint(
                fields=["is_owner"], condition=Q(is_owner=True), name="one_platform_owner"
            ),
            models.CheckConstraint(
                condition=Q(email=models.functions.Lower(models.functions.Trim("email"))),
                name="normalized_identity_email",
            ),
        ]


class Capability(models.Model):
    class Purpose(models.TextChoices):
        VERIFY = "verify"
        OWNER = "owner"
        ACTIVATE = "activate"
        RESET = "reset"

    identity = models.ForeignKey(AccessProfile, on_delete=models.CASCADE)
    purpose = models.CharField(max_length=12, choices=Purpose.choices)
    digest = models.CharField(max_length=64, unique=True)
    expires_at = models.DateTimeField()
    consumed_at = models.DateTimeField(null=True)
    revoked_at = models.DateTimeField(null=True)
    created_at = models.DateTimeField(auto_now_add=True)


class Notification(models.Model):
    class State(models.TextChoices):
        QUEUED = "queued"
        SENT = "sent", "Accepted by email backend (receipt unverified)"
        FAILED = "failed"

    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    identity = models.ForeignKey(AccessProfile, on_delete=models.CASCADE)
    deduplication_key = models.CharField(max_length=160, unique=True)
    kind = models.CharField(max_length=32)
    encrypted_payload = models.TextField()
    state = models.CharField(max_length=12, choices=State.choices, default=State.QUEUED)
    attempts = models.PositiveSmallIntegerField(default=0)
    error_code = models.CharField(max_length=40, blank=True)
    created_at = models.DateTimeField(auto_now_add=True)
    sent_at = models.DateTimeField(null=True)
    next_attempt_at = models.DateTimeField(null=True)


class AuditEvent(models.Model):
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    actor = models.ForeignKey(settings.AUTH_USER_MODEL, null=True, on_delete=models.PROTECT)
    identity = models.ForeignKey(AccessProfile, null=True, on_delete=models.PROTECT)
    action = models.CharField(max_length=48)
    occurred_at = models.DateTimeField(auto_now_add=True)
    # Codes and identifiers only; never request payloads or free text.
    code = models.CharField(max_length=48, blank=True)


class RateBucket(models.Model):
    key = models.CharField(max_length=64, primary_key=True)
    window = models.BigIntegerField()
    count = models.PositiveIntegerField(default=0)


class SupportTicket(models.Model):
    CATEGORIES = [
        (x, x.title()) for x in ("bug", "usability", "feature", "scientific_question", "interpretation")
    ]
    id = models.UUIDField(primary_key=True, default=uuid.uuid4, editable=False)
    identity = models.ForeignKey(AccessProfile, on_delete=models.PROTECT)
    category = models.CharField(max_length=24, choices=CATEGORIES)
    message = models.CharField(max_length=2000)
    status = models.CharField(
        max_length=12,
        choices=[("open", "Open"), ("answered", "Answered"), ("closed", "Closed")],
        default="open",
    )
    reply = models.CharField(max_length=2000, blank=True)
    internal_note = models.CharField(max_length=1000, blank=True)
    app_version = models.CharField(max_length=80)
    created_at = models.DateTimeField(auto_now_add=True)
    updated_at = models.DateTimeField(auto_now=True)
    submission_key = models.UUIDField(unique=True, default=uuid.uuid4)
