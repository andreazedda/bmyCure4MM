from __future__ import annotations

from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.db.models import Q
from django.shortcuts import get_object_or_404, redirect, render
from django.urls import reverse

from django.apps import apps

Regimen = apps.get_model("clinic", "Regimen")

from . import explain, forms, models
from .presets import PRESETS
from .models_help import HelpArticle
from .pharmaco import registry as pharmaco_registry
from .permissions import is_editor


def scenario_list(request):
    scenarios = models.Scenario.objects.filter(active=True).prefetch_related("recommended_regimens")
    context = {
        "scenarios": scenarios,
        "twin_assessment_id": (request.GET.get("twin_assessment_id") or "").strip(),
    }
    return render(request, "simulator/scenario_list.html", context)


@login_required
def getting_started(request):
    """Getting started page with tutorials and practice scenarios."""
    return render(request, "simulator/getting_started.html")


@login_required
def interactive_tutorial(request):
    """Interactive step-by-step tutorial for creating a patient and running first simulation."""
    from clinic.models import Patient, Assessment
    
    # Check if demo patients exist
    demo_patients = Patient.objects.filter(mrn__startswith='DEMO').prefetch_related('assessments')
    has_demo_data = demo_patients.exists()
    
    context = {
        'has_demo_data': has_demo_data,
        'demo_patients': demo_patients,
        'demo_patients_count': demo_patients.count(),
    }
    return render(request, "simulator/interactive_tutorial.html", context)


@login_required
def visibility_diagnostics(request):
    """Explain what assessments a user can/can't use (and why)."""
    from clinic.models import Assessment
    from .access import DEMO_MRN_PREFIX, accessible_assessments

    all_qs = Assessment.objects.select_related("patient")
    total_assessments = all_qs.count()

    is_privileged = request.user.is_staff or is_editor(request.user)
    if is_privileged:
        accessible_count = total_assessments
        inaccessible_count = 0
        inaccessible_sample = []
    else:
        access_q = Q(patient__owner=request.user) | Q(patient__mrn__startswith=DEMO_MRN_PREFIX)
        accessible_count = accessible_assessments(request.user, base_qs=all_qs).count()
        inaccessible_count = all_qs.exclude(access_q).count()
        inaccessible_sample = list(all_qs.exclude(access_q).order_by("-date")[:50])

    context = {
        "is_privileged": is_privileged,
        "total_assessments": total_assessments,
        "accessible_count": accessible_count,
        "inaccessible_count": inaccessible_count,
        "inaccessible_sample": inaccessible_sample,
        "demo_prefix": DEMO_MRN_PREFIX,
    }
    return render(request, "simulator/visibility_diagnostics.html", context)


@login_required
def scenario_detail(request, pk: int):
    scenario = get_object_or_404(
        models.Scenario.objects.prefetch_related("recommended_regimens", "attempts__selected_regimen", "attempts__user"),
        pk=pk,
        active=True,
    )
    attempts = scenario.attempts.all()
    form = forms.SimulationAttemptForm(request.POST or None)
    if request.method == "POST" and form.is_valid():
        attempt = form.save(commit=False)
        attempt.scenario = scenario
        attempt.user = request.user
        regimen = attempt.selected_regimen
        if regimen and scenario.recommended_regimens.filter(pk=regimen.pk).exists():
            attempt.is_guideline_aligned = True
        attempt.save()
        if attempt.is_guideline_aligned:
            messages.success(
                request,
                "Simulation submitted. Great choice—aligned with the guideline set for this scenario.",
            )
        else:
            messages.warning(
                request,
                "Simulation submitted. Review the guideline notes below to compare approaches.",
            )
        return redirect(reverse("simulator:scenario_detail", args=[scenario.pk]))

    recommended_regimens = scenario.recommended_regimens.order_by("name")
    regimen_names = ", ".join(regimen.name for regimen in recommended_regimens) if recommended_regimens else "No guideline regimen linked yet."
    editor = is_editor(request.user)
    available_regimens = (
        Regimen.objects.exclude(pk__in=recommended_regimens.values_list("pk", flat=True)).order_by("name")
        if editor
        else Regimen.objects.none()
    )
    latest_simulation = (
        scenario.attempts.exclude(results_summary={}).order_by("-submitted").first()
    )
    latest_summary = latest_simulation.results_summary if latest_simulation else None
    latest_results = latest_simulation.results if latest_simulation else {}
    game_mode = request.GET.get("game") == "1"
    from .game import compute_game_metrics
    game = compute_game_metrics(latest_summary) if game_mode else None
    latest_warnings: list[str] = []
    if latest_summary:
        healthy_loss = latest_summary.get("healthy_loss")
        if healthy_loss is not None:
            if healthy_loss > 0.3:
                latest_warnings.append("Healthy cell loss exceeded 30% in the most recent simulation—dose reduction recommended.")
            elif healthy_loss > 0.2:
                latest_warnings.append("Healthy cell loss above 20% warrants close monitoring for toxicity.")
        tumor_reduction = latest_summary.get("tumor_reduction")
        if tumor_reduction is not None and tumor_reduction < 0:
            latest_warnings.append("Latest simulation predicted tumor regrowth (negative reduction). Adjust parameters and re-run.")
    def _profile_with_ranges(drug: str):
        profile = pharmaco_registry.get_drug_profile(drug)
        if not profile:
            return None
        profile_copy = dict(profile)
        dose = profile_copy.get("dose_range", {})
        span = f"{dose.get('min')}–{dose.get('max')} {profile_copy.get('unit', '')}".strip()
        profile_copy["range_en"] = f"Allowed: {span}"
        profile_copy["range_it"] = f"Consentito: {span}"
        return profile_copy

    drug_profiles = {
        "lenalidomide": _profile_with_ranges("lenalidomide"),
        "bortezomib": _profile_with_ranges("bortezomib"),
        "daratumumab": _profile_with_ranges("daratumumab"),
    }
    guide_slugs = ["quickstart", "optimization_lab"]
    guides: dict[str, dict[str, dict[str, str]]] = {}
    articles = {article.slug: article for article in HelpArticle.objects.filter(slug__in=guide_slugs)}
    for slug in guide_slugs:
        article = articles.get(slug)
        if not article:
            continue
        guides[slug] = {
            "en": article.as_lang("en"),
            "it": article.as_lang("it"),
        }
    preset_descriptions = {
        key: {
            "en": preset.get("description_en", ""),
            "it": preset.get("description_it", ""),
            "story_en": preset.get("story_en", {}),
            "story_it": preset.get("story_it", {}),
        }
        for key, preset in PRESETS.items()
    }
    help_index = [
        {"slug": article.slug, "title_en": article.title_en, "title_it": article.title_it}
        for article in HelpArticle.objects.order_by("slug")
    ]

    twin_assessment_id = (request.GET.get("twin_assessment_id") or "").strip()
    sim_initial = {}
    if twin_assessment_id.isdigit():
        sim_initial["twin_assessment_id"] = int(twin_assessment_id)
    prefill_twin = bool(sim_initial)

    context = {
        "scenario": scenario,
        "attempts": attempts,
        "form": form,
        "recommended_regimens": recommended_regimens,
        "regimen_names": regimen_names,
        "is_editor": editor,
        "available_regimens": available_regimens,
        "simulation_parameter_form": forms.SimulationParameterForm(user=request.user, initial=sim_initial),
        "prefill_twin": prefill_twin,
        "latest_simulation": latest_simulation,
        "latest_simulation_summary": latest_summary,
        "latest_simulation_results": latest_results,
        "latest_simulation_warnings": latest_warnings,
        "game_mode": game_mode,
        "game": game,
        "drug_profiles": drug_profiles,
        "sim_form_help_it": forms.SIMULATION_FORM_HELP_TEXT_IT,
        "sim_form_help_en": forms.SIMULATION_FORM_HELP_TEXT_EN,
        "kpi": explain.KPI,
        "guide_articles": guides,
        "help_index": help_index,
        "preset_descriptions": preset_descriptions,
    }
    return render(request, "simulator/scenario_detail.html", context)
