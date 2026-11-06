import uuid
from functools import wraps

from django.contrib import messages
from django.contrib.auth import login, logout
from django.core.exceptions import PermissionDenied, ValidationError
from django.core.paginator import Paginator
from django.db import DatabaseError, transaction
from django.http import HttpResponse
from django.shortcuts import get_object_or_404, redirect, render
from django.utils import timezone
from django.views.decorators.debug import sensitive_post_parameters
from django.views.decorators.http import require_POST

from . import forms, services
from .backends import VerifiedEmailBackend
from .models import AccessPolicy, Notification, SupportTicket
from .models import AccessProfile as Identity
from .ui import bilingual


def limited(request, scope, identity_key="", limit=10):
    try:
        allowed = services.rate_limit(scope, request, identity_key, limit=limit)
    except DatabaseError:
        return render(
            request,
            "accounts/message.html",
            {
                "title": "Please try again later",
                "text": "The service is temporarily unavailable. No action was taken.",
            },
            status=503,
        )
    if not allowed:
        response = render(
            request,
            "accounts/message.html",
            {"title": "Please try again later", "text": "Too many attempts. Please try again in an hour."},
            status=429,
        )
        response["Retry-After"] = "3600"
        return response


def member_required(view):
    @wraps(view)
    def wrapped(request, *args, **kwargs):
        if not request.user.is_authenticated:
            return redirect("accounts:login")
        services.require_community(request.access_profile)
        return view(request, *args, **kwargs)

    return wrapped


def manager_required(view):
    @wraps(view)
    def wrapped(request, *args, **kwargs):
        if not request.user.is_authenticated:
            return redirect("accounts:login")
        services.require_manager(request.user)
        from .security import mfa_current

        if not mfa_current(request):
            return redirect("accounts:security")
        return view(request, *args, **kwargs)

    return wrapped


@sensitive_post_parameters()
def request_access(request):
    if not AccessPolicy.objects.get(pk=1).requests_open:
        return render(
            request,
            "accounts/message.html",
            {
                "title": "Registration is temporarily paused",
                "text": "New access requests are currently paused. No verification email is sent while registration is paused. Please try again later. Existing accounts can still sign in.",
            },
        )
    form = forms.RequestForm(request.POST or None)
    if request.method == "POST":
        throttle = limited(request, "request", request.POST.get("email", "").strip().lower(), 5)
        if throttle is not None:
            return throttle
        if form.is_valid():
            try:
                identity = services.request_access(form.cleaned_data)
            except ValidationError as exc:
                form.add_error(None, exc)
            else:
                if identity:
                    request.session["access_request"] = str(identity.pk)
                return redirect("accounts:submitted")
    return render(
        request,
        "accounts/form.html",
        {
            "title": "Request community access",
            "title_it": "Richiedi accesso alla piattaforma",
            "form": form,
            "button": "Send verification email",
            "button_it": "Invia email di verifica",
            "privacy_notice": True,
            "notice": "Confirm your email, wait for Andrea’s approval, then choose your password. Use synthetic data only; never include patient information.",
            "notice_it": "Conferma l’email, attendi l’approvazione di Andrea, poi scegli la password. Usa solo dati sintetici; non inserire informazioni su pazienti.",
        },
    )


def submitted(request):
    return render(
        request,
        "accounts/message.html",
        {
            "title": "Check your email",
            "text": "If this request can be processed, a verification email will arrive. Delivery may take a few minutes. Verification does not approve access.",
            "status_link": True,
            "hide_login": True,
        },
    )


def request_status(request):
    identity = (
        Identity.objects.filter(pk=request.session.get("access_request")).first()
        if request.session.get("access_request")
        else None
    )
    return render(
        request,
        "accounts/message.html",
        {
            "title": "Your access request",
            "next_step": (
                "Your email is verified. Andrea will review your request and notify you by email."
                if identity.request_state == "PENDING_REVIEW"
                else "Open the verification email to send your request for review."
                if identity.request_state == "PENDING_EMAIL_VERIFICATION"
                else "This request was not approved. Reply to the decision email if you need clarification."
                if identity.request_state == "REJECTED"
                else "Your request was approved. Use the activation email to set your password; ask Andrea to reissue it if expired."
                if identity.access_state == "PENDING_ACTIVATION"
                else "Your access is active. Sign in to continue."
                if identity.access_state == "ACTIVE"
                else "Access is unavailable. Reply to an account email for help. Password reset does not renew access."
            )
            if identity
            else "",
            "hide_login": bool(
                identity and (identity.request_state != "APPROVED" or identity.access_state != "ACTIVE")
            ),
            "text": identity.get_request_state_display()
            if identity
            else "Request status is available only in the browser that submitted or verified the request.",
        },
    )


@sensitive_post_parameters()
def capability(request, purpose):
    form_type = forms.CapabilityForm if purpose == "verify" else forms.ActivationForm
    form = form_type(request.POST or None)
    if request.method == "POST":
        throttle = limited(request, "capability", limit=20)
        if throttle is not None:
            return throttle
        if form.is_valid():
            try:
                identity = services.confirm_capability(
                    form.cleaned_data["token"], purpose, form.cleaned_data.get("password")
                )
            except ValidationError as exc:
                form.add_error(None, exc)
            else:
                if purpose == "verify":
                    request.session["access_request"] = str(identity.pk)
                    return redirect("accounts:status")
                messages.success(
                    request, bilingual("Password saved. You can now sign in with your verified email.")
                )
                return redirect("accounts:login")
    titles = {
        "verify": "Confirm your email",
        "activate": "Activate your community access",
        "owner": "Set the platform owner password",
        "reset": "Choose a new password",
    }
    # Never echo posted tokens after an error; the user must reopen the email link.
    form.data = form.data.copy() if form.is_bound else form.data
    if form.is_bound:
        form.data["token"] = ""
    return render(
        request,
        "accounts/form.html",
        {
            "title": titles[purpose],
            "form": form,
            "button": "Confirm email" if purpose == "verify" else "Save password",
            "token_page": True,
            "notice": "Confirm to send your request for review. Verification does not approve access."
            if purpose == "verify"
            else "This action requires confirmation. Email previews do not activate access. If the form fails, reopen the email link. Use at least 12 characters for your password.",
        },
    )


@sensitive_post_parameters()
def sign_in(request):
    form = forms.LoginForm(request.POST or None)
    if request.method == "POST":
        throttle = limited(request, "login", request.POST.get("email", "").strip().lower(), 20)
        if throttle is not None:
            return throttle
        if form.is_valid():
            user = VerifiedEmailBackend().authenticate(request, **form.cleaned_data)
            identity = Identity.objects.filter(user=user).first() if user else None
            if user:
                login(request, user, backend="accounts.backends.VerifiedEmailBackend")
                request.session["access_version"] = identity.session_version
                services.audit("product_login", identity, user)
                return redirect(
                    "accounts:manage_requests"
                    if identity.is_owner or identity.is_approver
                    else "simulator:getting_started"
                )
            form.add_error(
                None, "Unable to sign in. Check your credentials, approval, expiry and access network."
            )
    return render(
        request,
        "accounts/form.html",
        {
            "title": "Sign in",
            "title_it": "Accedi",
            "form": form,
            "button": "Sign in",
            "button_it": "Accedi",
            "login_page": True,
        },
    )


@require_POST
def sign_out(request):
    logout(request)
    return redirect("docs_viewer:index")


@sensitive_post_parameters()
def reset_request(request):
    form = forms.EmailForm(request.POST or None)
    if request.method == "POST":
        throttle = limited(request, "reset", request.POST.get("email", "").strip().lower(), 5)
        if throttle is not None:
            return throttle
        if form.is_valid():
            with transaction.atomic():
                services.locked_policy()
                identity = (
                    Identity.objects.select_related("user")
                    .filter(email=form.cleaned_data["email"].strip().lower())
                    .first()
                )
                if services.usable(identity) and identity.user.is_active:
                    services.issue_capability(identity, "reset")
            return render(
                request,
                "accounts/message.html",
                {
                    "title": "Check your email",
                    "text": "If an eligible active account exists, password reset instructions will be sent. Reset does not change approval or access status.",
                },
            )
    return render(
        request,
        "accounts/form.html",
        {
            "title": "Reset your password",
            "title_it": "Reimposta la password",
            "form": form,
            "button": "Request reset",
            "button_it": "Richiedi nuova password",
        },
    )


@member_required
def feedback(request):
    form = forms.SupportTicketForm(request.POST or None, initial={"submission_key": uuid.uuid4()})
    if request.method == "POST":
        throttle = limited(request, "feedback", str(request.access_profile.pk), 10)
        if throttle is not None:
            return throttle
        if form.is_valid():
            services.submit_feedback(
                request.access_profile, form.cleaned_data, form.cleaned_data["submission_key"]
            )
            messages.success(request, bilingual("Your ticket has been saved. Replies will appear below."))
            return redirect("accounts:feedback")
    return render(
        request,
        "accounts/feedback.html",
        {
            "form": form,
            "tickets": Paginator(
                SupportTicket.objects.filter(identity=request.access_profile).order_by("-created_at", "-pk"),
                20,
            ).get_page(request.GET.get("page")),
        },
    )


@manager_required
def manage_requests(request):
    form = forms.DecisionForm(request.POST or None)
    if form.is_bound:
        try:
            form.auto_id = f"%s-{uuid.UUID(request.POST.get('identity_id', ''))}"
        except (ValueError, AttributeError):
            pass
    if request.method == "POST" and form.is_valid():
        throttle = limited(request, "management", str(request.user.pk), 100)
        if throttle is not None:
            return throttle
        data = form.cleaned_data
        if data["action"] == "revoke" and not data["confirm_revoke"]:
            identity = get_object_or_404(Identity, pk=data["identity_id"], is_owner=False, is_approver=False)
            return render(request, "accounts/revoke_confirm.html", {"identity": identity})
        try:
            if data["action"] in {"approve", "reject"}:
                services.decide_request(
                    request.user,
                    data["identity_id"],
                    approve=data["action"] == "approve",
                    days=data["days"] or 14,
                )
            elif data["action"] == "resend":
                services.resend_activation(request.user, data["identity_id"])
            elif data["action"] == "renew":
                services.renew_access(request.user, data["identity_id"], data["days"] or 14)
            else:
                services.change_access(
                    request.user,
                    data["identity_id"],
                    "SUSPENDED" if data["action"] == "suspend" else "REVOKED",
                )
        except (ValidationError, Identity.DoesNotExist) as exc:
            form.add_error(None, exc if isinstance(exc, ValidationError) else "Request unavailable.")
        else:
            messages.success(request, bilingual("Access decision saved."))
            return redirect("accounts:manage_requests")
    rows = Paginator(
        Identity.objects.filter(verified_at__isnull=False, is_owner=False).order_by("-created_at", "-pk"), 20
    ).get_page(request.GET.get("page"))
    return render(
        request,
        "accounts/manage.html",
        {
            "rows": rows,
            "submitted_identity_id": request.POST.get("identity_id"),
            "submitted_days": request.POST.get("days", "14"),
            "form": form,
            "policy": AccessPolicy.objects.get(pk=1),
            "failed": Paginator(
                Notification.objects.select_related("identity")
                .filter(state="failed")
                .order_by("-created_at", "-pk"),
                10,
            ).get_page(request.GET.get("failed_page")),
            "switch_form": forms.SwitchForm(
                initial={k: getattr(AccessPolicy.objects.get(pk=1), k) for k in forms.SwitchForm.base_fields}
            ),
        },
    )


@manager_required
def manage_users(request):
    rows = Paginator(
        Identity.objects.filter(user__isnull=False, is_owner=False).order_by("-created_at", "-pk"), 20
    ).get_page(request.GET.get("page"))
    return render(request, "accounts/users.html", {"rows": rows})


@require_POST
@manager_required
def switches(request):
    if not request.access_profile.is_owner:
        raise PermissionDenied
    form = forms.SwitchForm(request.POST)
    if not form.is_valid():
        return HttpResponse("Invalid control settings.", status=400)
    with transaction.atomic():
        policy = services.locked_policy()
        for key, value in form.cleaned_data.items():
            setattr(policy, key, value)
        policy.save()
        services.audit("pilot_controls_changed", request.access_profile, request.user)
    messages.success(request, bilingual("Registration controls saved."))
    return redirect("accounts:manage_requests")


@require_POST
@manager_required
def retry_notification(request):
    throttle = limited(request, "notification_retry", str(request.user.pk), 20)
    if throttle is not None:
        return throttle
    form = forms.RetryForm(request.POST)
    if not form.is_valid():
        return HttpResponse("Invalid retry request.", status=400)
    with transaction.atomic():
        notification = get_object_or_404(
            Notification.objects.select_for_update(), pk=form.cleaned_data["notification_id"]
        )
        if notification.state == "failed" and notification.encrypted_payload:
            notification.state, notification.next_attempt_at, notification.error_code = (
                "queued",
                timezone.now(),
                "",
            )
            notification.attempts = 0
            notification.save()
            from .tasks import request_delivery

            transaction.on_commit(request_delivery)
            services.audit("notification_retry", notification.identity, request.user)
            messages.success(request, bilingual("Notification queued for another delivery attempt."))
    return redirect("accounts:manage_requests")


@manager_required
def manage_feedback(request):
    form = forms.ReplyForm(request.POST or None)
    if request.method == "POST" and form.is_valid():
        with transaction.atomic():
            services.locked_policy()
            ticket = get_object_or_404(
                SupportTicket.objects.select_for_update(), pk=form.cleaned_data["ticket_id"]
            )
            changed = any(
                getattr(ticket, key) != form.cleaned_data[key]
                for key in ("category", "reply", "internal_note", "status")
            )
            if changed:
                previous_reply = ticket.reply
                for key in ("category", "reply", "internal_note", "status"):
                    setattr(ticket, key, form.cleaned_data[key])
                ticket.save()
                services.audit("feedback_updated", ticket.identity, request.user)
                if ticket.reply and ticket.reply != previous_reply:
                    services.notify(
                        ticket.identity,
                        "feedback_reply",
                        f"reply:{ticket.pk}:{ticket.updated_at.isoformat()}",
                        "Your feedback has a reply. Sign in and open /feedback/ to read it.",
                    )
        messages.success(request, bilingual("Reply and ticket status saved."))
        return redirect("accounts:manage_feedback")
    tickets = Paginator(
        SupportTicket.objects.select_related("identity").order_by("-created_at", "-pk"), 20
    ).get_page(request.GET.get("page"))
    bound_visible = False
    for ticket in tickets:
        if form.is_bound and request.POST.get("ticket_id") == str(ticket.pk):
            ticket.reply_form = form
            bound_visible = True
        else:
            ticket.reply_form = forms.ReplyForm(
                initial={
                    "ticket_id": ticket.pk,
                    "category": ticket.category,
                    "reply": ticket.reply,
                    "internal_note": ticket.internal_note,
                    "status": ticket.status,
                }
            )
        ticket.reply_form.auto_id = f"ticket-{ticket.pk}-%s"
    return render(
        request,
        "accounts/manage_feedback.html",
        {"tickets": tickets, "form": form, "bound_ticket_visible": bound_visible},
    )


@sensitive_post_parameters()
def owner_security(request):
    if not request.user.is_authenticated:
        return redirect("accounts:login")
    profile = services.require_manager(request.user)
    import time

    from .security import confirm_owner_code, owner_secret, owner_setup_qr

    form = forms.OTPForm(request.POST or None)
    secret = owner_secret(profile)
    if request.method == "POST":
        throttle = limited(request, "owner_mfa", str(profile.pk), 10)
        if throttle is not None:
            return throttle
        if form.is_valid() and confirm_owner_code(profile, form.cleaned_data["code"]):
            request.session["owner_mfa_version"] = profile.session_version
            request.session["owner_mfa_at"] = time.time()
            return redirect("accounts:manage_requests")
        form.add_error(None, "Invalid or already used authenticator code.")
    return render(
        request,
        "accounts/security.html",
        {
            "form": form,
            "setup_secret": secret if not profile.totp_confirmed else None,
            "setup_qr": owner_setup_qr(profile, secret) if not profile.totp_confirmed else None,
        },
    )
