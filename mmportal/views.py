from __future__ import annotations

from django.http import HttpRequest, HttpResponse
from django.shortcuts import redirect
from django.views.generic import TemplateView

from clinic.views import dashboard


class DocsIndexView(TemplateView):
    template_name = "docs/index.html"


def home(request: HttpRequest) -> HttpResponse:
    """Public landing page.

    - Anonymous users: redirect to docs viewer (public-safe demo)
    - Authenticated users: show the clinic dashboard
    """

    if request.user.is_authenticated:
        return dashboard(request)

    return redirect("docs_viewer:index")
