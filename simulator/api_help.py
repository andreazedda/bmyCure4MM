from __future__ import annotations

from django.db.models import Q
from django.http import Http404, HttpResponse, JsonResponse
from django.utils.crypto import salted_hmac
from django.utils.http import http_date
from django.views.decorators.http import require_GET

from .forms import SIMULATION_FORM_HELP_TEXT_EN, SIMULATION_FORM_HELP_TEXT_IT
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


HELP_ALIASES = {
    # Tour / UI aliases -> canonical slugs
    "twin": "use_twin",
    "dose_ranges_lenalidomide": "lenalidomide_dose",
    "dose_ranges_bortezomib": "bortezomib_dose",
    "dose_ranges_daratumumab": "daratumumab_dose",
}


def _field_payload(slug: str, lang: str) -> dict[str, str] | None:
    help_map = SIMULATION_FORM_HELP_TEXT_IT if lang == "it" else SIMULATION_FORM_HELP_TEXT_EN
    text = help_map.get(slug)
    if not text:
        return None

    title = slug.replace("_", " ").strip().title()
    if slug == "use_twin":
        title = "Gemello Paziente (Patient Twin)" if lang == "it" else "Patient Twin"

    body_parts: list[str] = [f"<p>{text}</p>"]

    if slug == "use_twin":
        if lang == "it":
            body_parts.append(
                "<p><strong>Nota:</strong> il Twin non appare nella lista pazienti come oggetto separato. "
                "È un calcolo fatto a runtime usando uno <em>Assessment</em> (snapshot lab/clinico).</p>"
            )
            body_parts.append(
                "<h6>Come usarlo (in breve)</h6>"
                "<ol>"
                "<li>Clinica → Pazienti → apri un paziente → <strong>Add Assessment</strong>.</li>"
                "<li>Simulatore → sezione <strong>Advanced &amp; Gemello</strong> → abilita <strong>Use Twin</strong>.</li>"
                "<li>Seleziona l'Assessment e metti <strong>Auto</strong> per derivare la biologia dai lab.</li>"
                "</ol>"
            )
        else:
            body_parts.append(
                "<p><strong>Note:</strong> the Twin is not a separate object in the patient list. "
                "It is computed at runtime from an <em>Assessment</em> (lab/clinical snapshot).</p>"
            )
            body_parts.append(
                "<h6>How to use it</h6>"
                "<ol>"
                "<li>Clinic → Patients → open a patient → <strong>Add Assessment</strong>.</li>"
                "<li>Simulator → <strong>Advanced &amp; Twin</strong> → enable <strong>Use Twin</strong>.</li>"
                "<li>Select the Assessment and choose <strong>Auto</strong> to derive biology from labs.</li>"
                "</ol>"
            )

    return {"title": title, "body": "\n".join(body_parts)}


def _static_article_payload(slug: str, lang: str) -> dict[str, str] | None:
    if slug != "quickstart":
        return None

    if lang == "it":
        title = "Guida rapida"
        body = (
            "<p>Per fare un test completo (incluso Patient Twin) in pochi minuti:</p>"
            "<ol>"
            "<li>Crea un paziente: Clinica → Pazienti → Nuovo.</li>"
            "<li>Aggiungi almeno un Assessment al paziente (serve al Twin).</li>"
            "<li>Vai al simulatore (/sim/…) e apri il form.</li>"
            "<li>In <strong>Advanced &amp; Gemello</strong>, abilita <strong>Use Twin</strong>, seleziona l'Assessment e scegli <strong>Auto</strong>.</li>"
            "<li>Esegui la simulazione e controlla nei risultati <code>twin_params.json</code>.</li>"
            "</ol>"
        )
    else:
        title = "Quickstart"
        body = (
            "<p>End-to-end run (including Patient Twin) in a few minutes:</p>"
            "<ol>"
            "<li>Create a patient: Clinic → Patients → New.</li>"
            "<li>Add at least one Assessment to the patient (required for the Twin).</li>"
            "<li>Open the simulator (/sim/…) and the simulation form.</li>"
            "<li>Under <strong>Advanced &amp; Twin</strong>, enable <strong>Use Twin</strong>, pick the Assessment and choose <strong>Auto</strong>.</li>"
            "<li>Run the simulation and check the results for <code>twin_params.json</code>.</li>"
            "</ol>"
        )
    return {"title": title, "body": body}


def infer_type(slug: str) -> str:
    if slug.startswith("kpi_"):
        return "kpi"
    if slug in FIELD_SLUGS:
        return "field"
    if slug in PRESETS:
        return "preset"
    return "article"


def _preset_payload(slug: str, lang: str) -> dict[str, str]:
    preset = PRESETS[slug]
    title = str(preset.get("label") or slug)
    description = str(preset.get(f"description_{lang}") or preset.get("description_en") or "")
    story = preset.get(f"story_{lang}") or preset.get("story_en") or {}

    body_parts: list[str] = []
    if description:
        body_parts.append(f"<p>{description}</p>")

    if isinstance(story, dict) and story:
        story_title = story.get("title")
        if story_title:
            body_parts.append(f"<h5>{story_title}</h5>")

        for key, heading_en, heading_it in (
            ("background", "Background", "Contesto"),
            ("challenge", "Challenge", "Sfida"),
            ("why_these_drugs", "Why these drugs", "Perché questi farmaci"),
            ("expected_outcome", "Expected outcome", "Risultato atteso"),
        ):
            value = story.get(key)
            if not value:
                continue
            heading = heading_it if lang == "it" else heading_en
            body_parts.append(f"<h6>{heading}</h6><p>{value}</p>")

        learning_points = story.get("learning_points")
        if isinstance(learning_points, list) and learning_points:
            heading = "Learning points" if lang != "it" else "Punti chiave"
            items = "".join(f"<li>{point}</li>" for point in learning_points if point)
            if items:
                body_parts.append(f"<h6>{heading}</h6><ul>{items}</ul>")

    body = "\n".join(body_parts) if body_parts else "<p>(no details)</p>"
    return {"title": title, "body": body}


def help_item(request, slug: str):
    lang = request.GET.get("lang", "en")

    slug = HELP_ALIASES.get(slug, slug)

    static_payload = _static_article_payload(slug, lang)
    if static_payload:
        static_payload["type"] = infer_type(slug)
        static_payload["slug"] = slug
        resp = JsonResponse(static_payload)
        resp["Cache-Control"] = "max-age=600, public"
        return resp

    field_payload = _field_payload(slug, lang)
    if field_payload:
        field_payload["type"] = infer_type(slug)
        field_payload["slug"] = slug
        resp = JsonResponse(field_payload)
        resp["Cache-Control"] = "max-age=600, public"
        return resp

    if slug in PRESETS:
        payload = _preset_payload(slug, lang)
        payload["type"] = "preset"
        payload["slug"] = slug
        resp = JsonResponse(payload)
        resp["Cache-Control"] = "max-age=600, public"
        return resp

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
                        "slug": slug,
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
