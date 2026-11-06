from django.contrib.auth import logout
from django.http import HttpResponse
from django.shortcuts import redirect, render
from django.utils.deprecation import MiddlewareMixin

from .models import AccessProfile
from .services import usable


class ProductAccessMiddleware(MiddlewareMixin):
    def process_view(self, request, view_func, args, kwargs):
        profile = (
            AccessProfile.objects.select_related("user").filter(user=request.user).first()
            if request.user.is_authenticated
            else None
        )
        request.access_profile = profile
        match = request.resolver_match
        if profile:
            legacy_owner_session = (
                profile.is_owner
                and profile.user.is_staff
                and request.session.get("_auth_user_backend") == "django.contrib.auth.backends.ModelBackend"
            )
            if not usable(profile) or (
                not legacy_owner_session and request.session.get("access_version") != profile.session_version
            ):
                logout(request)
                request.access_profile = None
                if match.namespace == "accounts":
                    return None
                return redirect("accounts:login")
        if match.namespace == "accounts":
            try:
                if int(request.META.get("CONTENT_LENGTH", 0) or 0) > 16384:
                    return HttpResponse("Request too large. Attachments are not accepted.", status=413)
            except ValueError:
                return HttpResponse("Invalid request.", status=400)
            return None
        if not profile:
            return None  # Existing production users and public routes retain their established behavior.
        if profile.is_owner or profile.is_approver:
            if match.namespace == "admin":
                from .security import mfa_current

                if not mfa_current(request):
                    return redirect("accounts:security")
            return None
        if request.path == "/":
            return redirect("simulator:getting_started")
        if match.namespace == "docs_viewer":
            return None
        allowed = {
            "scenario_list",
            "getting_started",
            "scenario_detail",
            "algorithm_transparency",
            "scenario_simulate",
            "attempt_export",
            "attempt_plot_embed",
        }
        if match.namespace in {"simulator", "simulator_alias", "learn"} and match.url_name in allowed:
            if "twin_assessment_id" in request.GET or "assessment_id" in request.GET:
                return HttpResponse("Patient context is not available for this account.", status=403)
            return None
        if match.url_name in {
            "api_help_item",
            "api_help_search",
            "api_glossary",
            "api_drugs",
        } and request.method in {"GET", "HEAD"}:
            return None
        return render(request, "accounts/restricted.html", status=403)

    def process_response(self, request, response):
        match = getattr(request, "resolver_match", None)
        if match and match.namespace == "accounts":
            response["Cache-Control"] = "no-store, private"
            response["Referrer-Policy"] = "same-origin"
        return response
