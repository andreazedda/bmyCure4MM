"""mmportal URL configuration."""

from __future__ import annotations

from django.conf import settings
from django.conf.urls.static import static
from django.contrib import admin
from django.http import JsonResponse
from django.urls import include, path

from simulator.api import drugs, glossary
from .views import DocsIndexView

urlpatterns = [
    path("admin/", admin.site.urls),
    path("", include(("clinic.urls", "clinic"), namespace="clinic")),
    path("chem/", include(("chemtools.urls", "chemtools"), namespace="chemtools")),
    path("sim/", include(("simulator.urls", "simulator"), namespace="simulator")),
    path("docs/", DocsIndexView.as_view(), name="docs"),
    path("api/", include("clinic.api")),
    path("api/glossary/", glossary, name="api_glossary"),
    path("api/drugs/", drugs, name="api_drugs"),
    path("api/chem/", include("chemtools.api")),
    path("healthz/", lambda request: JsonResponse({"ok": True})),
]

if settings.DEBUG:
    urlpatterns += static(settings.MEDIA_URL, document_root=settings.MEDIA_ROOT)
