from __future__ import annotations

from django.db.models import Q
from django.http import Http404, HttpResponse, JsonResponse
from django.utils.crypto import salted_hmac
from django.utils.http import http_date
from django.utils.text import slugify
from django.views.decorators.http import require_GET

from .models_help import HelpArticle
from .presets import PRESETS

FIELD_SLUGS = {
    "preset",
    "baseline_tumor_cells",
    "baseline_healthy_cells",
    "creatinine_clearance",
    "neuropathy_grade",
    "anc",
    "platelets",
    "pregnancy",
    "lenalidomide_dose",
    "bortezomib_dose",
    "daratumumab_dose",
    "time_horizon",
    "cohort_size",
    "tumor_growth_rate",
    "healthy_growth_rate",
    "interaction_strength",
    "use_twin",
    "seed",
}


def infer_type(slug: str) -> str:
    if slug.startswith("kpi_"):
        return "kpi"
    if slug in FIELD_SLUGS:
        return "field"
    if slug in PRESETS:
        return "preset"
    return "article"


def help_item(request, slug: str):
    lang = request.GET.get("lang", "en")
    try:
        article = HelpArticle.objects.get(slug=slug)
    except HelpArticle.DoesNotExist as exc:  # pragma: no cover - 404 path
        raise Http404() from exc
    payload = article.as_lang(lang)
    payload["type"] = infer_type(article.slug)
    payload["slug"] = article.slug

    # ETag for caching
    etag = salted_hmac(
        "help", f"{slug}:{lang}:{article.updated.isoformat()}"
    ).hexdigest()
    if request.headers.get("If-None-Match") == etag:
        return HttpResponse(status=304)

    resp = JsonResponse(payload)
    resp["ETag"] = etag
    resp["Cache-Control"] = "max-age=600, public"
    resp["Last-Modified"] = http_date(article.updated.timestamp())
    return resp


@require_GET
def help_search(request):
    lang = request.GET.get("lang", "en")
    query = (request.GET.get("q") or "").strip()
    qs = HelpArticle.objects.all()
    items: list[dict[str, str]] = []

    if query:
        lowered = query.lower()
        # Search articles
        for article in qs:
            data = article.as_lang(lang)
            title = data["title"]
            if lowered in article.slug.lower() or lowered in title.lower():
                items.append(
                    {
                        "slug": article.slug,
                        "title": title,
                        "type": infer_type(article.slug),
                    }
                )

        # Search presets - prioritize exact matches
        for slug, preset in PRESETS.items():
            label = preset.get("label", slug)
            if lowered in slug.lower() or lowered in label.lower():
                items.insert(
                    0,
                    {
                        "slug": slugify(slug),
                        "title": label,
                        "type": "preset",
                    },
                )
                break
    else:
        # No query: return first 20 articles
        for article in qs.order_by("slug")[:20]:
            data = article.as_lang(lang)
            items.append(
                {
                    "slug": article.slug,
                    "title": data["title"],
                    "type": infer_type(article.slug),
                }
            )

    # Deduplication
    seen = set()
    dedup = []
    for it in items:
        if it["slug"] in seen:
            continue
        seen.add(it["slug"])
        dedup.append(it)

    # Cutoff to 20 results
    return JsonResponse({"results": dedup[:20]})
