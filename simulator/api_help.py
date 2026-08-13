from __future__ import annotations

from django.db.models import Q
from django.http import Http404, HttpResponse, JsonResponse
from django.utils.crypto import salted_hmac
from django.utils.http import http_date
from django.views.decorators.http import require_GET

from mmportal.governance import EpistemicLabel, governance_metadata

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


GOVERNED_STATIC_HELP = {
    "quickstart": {
        "en": (
            "Research quickstart",
            "<p>Choose a synthetic scenario, record a hypothetical configuration and seed, run the model, then inspect lineage, uncertainty, and limitations.</p>",
        ),
        "it": (
            "Avvio rapido di ricerca",
            "<p>Scegli uno scenario sintetico, registra configurazione ipotetica e seed, esegui il modello, poi esamina lineage, incertezza e limiti.</p>",
        ),
    },
    "doses": {
        "en": (
            "Hypothetical exposure inputs",
            "<p>Displayed bounds delimit the configured model domain. They are not clinical safety ranges and do not authorize a dose change.</p>",
        ),
        "it": (
            "Input ipotetici di esposizione",
            "<p>I limiti mostrati delimitano il dominio del modello. Non sono intervalli di sicurezza clinica e non autorizzano modifiche di dose.</p>",
        ),
    },
    "horizon_cohort": {
        "en": (
            "Simulation horizon and virtual cohort",
            "<p>Long horizons can increase numerical stiffness. Cohort size controls virtual sampling; report seed and uncertainty.</p>",
        ),
        "it": (
            "Orizzonte e coorte virtuale",
            "<p>Orizzonti lunghi possono aumentare la rigidità numerica. La dimensione della coorte controlla il campionamento virtuale; riporta seed e incertezza.</p>",
        ),
    },
    "optimization_lab": {
        "en": (
            "Exploratory parameter search",
            "<p>The search compares model objectives and heuristic constraints. Pareto membership is not evidence of clinical safety, efficacy, or superiority.</p>",
        ),
        "it": (
            "Ricerca esplorativa dei parametri",
            "<p>La ricerca confronta obiettivi del modello e vincoli euristici. L'appartenenza al fronte Pareto non prova sicurezza, efficacia o superiorità clinica.</p>",
        ),
    },
    "kpi_tumor_reduction": {
        "en": (
            "Simulated tumor-state change",
            "<p>Computed as 1 - endpoint/start for the model state. It is SIMULATED and is not a validated clinical response category.</p>",
        ),
        "it": (
            "Variazione simulata dello stato tumorale",
            "<p>Calcolata come 1 - fine/inizio per lo stato del modello. È SIMULATED e non è una categoria di risposta clinica validata.</p>",
        ),
    },
    "kpi_healthy_loss": {
        "en": (
            "Simulated healthy-cell-state change",
            "<p>Computed from a simplified model state. Thresholds are heuristic flags, not validated toxicity endpoints or dose instructions.</p>",
        ),
        "it": (
            "Variazione simulata dello stato delle cellule sane",
            "<p>Calcolata da uno stato semplificato del modello. Le soglie sono segnali euristici, non endpoint di tossicità validati o istruzioni di dose.</p>",
        ),
    },
    "kpi_auc": {
        "en": (
            "Simulated exposure summary",
            "<p>AUC is derived from the configured exposure model. Interpretation requires the schedule, units, assumptions, and model version.</p>",
        ),
        "it": (
            "Sintesi simulata dell'esposizione",
            "<p>AUC deriva dal modello di esposizione configurato. L'interpretazione richiede calendario, unità, assunzioni e versione del modello.</p>",
        ),
    },
    "kpi_time_to_recurrence": {
        "en": (
            "Model-relative recurrence time",
            "<p>The first configured threshold crossing after the model nadir. It is a simulated horizon result, not patient prognosis.</p>",
        ),
        "it": (
            "Tempo di recidiva relativo al modello",
            "<p>Primo superamento della soglia configurata dopo il nadir del modello. È un risultato simulato, non una prognosi del paziente.</p>",
        ),
    },
}


def _static_article_payload(slug: str, lang: str) -> dict[str, object] | None:
    article = GOVERNED_STATIC_HELP.get(slug)
    if article is None:
        return None
    title, body = article["it" if lang == "it" else "en"]
    return {"title": title, "body": body}


def _add_governance(payload: dict[str, object], *, output_kind: str) -> None:
    payload["governance"] = governance_metadata(
        epistemic_label=EpistemicLabel.HYPOTHETICAL,
        output_kind=output_kind,
    )


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
    if lang == "it":
        body = (
            "<p>Preset ipotetico di input e calendario del modello. Il nome non "
            "costituisce una selezione terapeutica e i valori non sono limiti di sicurezza clinica.</p>"
        )
    else:
        body = (
            "<p>Hypothetical model-input and schedule preset. Its name is not a "
            "treatment selection and its values are not clinical safety limits.</p>"
        )
    return {"title": title, "body": body}


def help_item(request, slug: str):
    lang = request.GET.get("lang", "en")

    slug = HELP_ALIASES.get(slug, slug)

    static_payload = _static_article_payload(slug, lang)
    if static_payload:
        static_payload["type"] = infer_type(slug)
        static_payload["slug"] = slug
        _add_governance(static_payload, output_kind="governed_help_article")
        resp = JsonResponse(static_payload)
        resp["Cache-Control"] = "max-age=600, public"
        return resp

    field_payload = _field_payload(slug, lang)
    if field_payload:
        field_payload["type"] = infer_type(slug)
        field_payload["slug"] = slug
        _add_governance(field_payload, output_kind="model_input_help")
        resp = JsonResponse(field_payload)
        resp["Cache-Control"] = "max-age=600, public"
        return resp

    if slug in PRESETS:
        payload = _preset_payload(slug, lang)
        payload["type"] = "preset"
        payload["slug"] = slug
        _add_governance(payload, output_kind="hypothetical_preset_help")
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
    _add_governance(payload, output_kind="legacy_help_article")

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
    return JsonResponse({
        "governance": governance_metadata(
            epistemic_label=EpistemicLabel.DERIVED,
            output_kind="help_search_index",
        ),
        "results": dedup[:20],
    })
