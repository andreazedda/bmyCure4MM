from __future__ import annotations

from django.views.generic import TemplateView


class DocsIndexView(TemplateView):
    template_name = "docs/index.html"
