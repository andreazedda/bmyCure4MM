"""mmportal URL configuration."""

from __future__ import annotations

from django.conf import settings
from django.conf.urls.static import static
from django.contrib import admin
from django.http import JsonResponse
from django.urls import include, path

from simulator.api import drugs, glossary
from simulator.api_design import design_report
from simulator.api_help import help_item, help_search
from simulator.api_ux import audit as ux_audit

urlpatterns = [
    path("admin/", admin.site.urls),
    path("", include(("clinic.urls", "clinic"), namespace="clinic")),
    path("chem/", include(("chemtools.urls", "chemtools"), namespace="chemtools")),
    path("sim/", include(("simulator.urls", "simulator"), namespace="simulator")),
    path("simulator/", include(("simulator.urls", "simulator"), namespace="simulator_alias")),
    path("docs/", include(("docs_viewer.urls", "docs_viewer"), namespace="docs_viewer")),
    path("api/", include("clinic.api")),
    path("api/glossary/", glossary, name="api_glossary"),
    path("api/drugs/", drugs, name="api_drugs"),
    path("api/help/search/", help_search, name="api_help_search"),
    path("api/help/<slug:slug>/", help_item, name="api_help_item"),
    path("api/ux/audit/", ux_audit, name="api_ux_audit"),
    path("api/design-report/", design_report, name="api_design_report"),
    path("api/chem/", include("chemtools.api")),
    path("healthz/", lambda request: JsonResponse({"ok": True})),
]

if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
