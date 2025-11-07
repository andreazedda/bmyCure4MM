from __future__ import annotations

from typing import Iterable

from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.contrib.auth.mixins import LoginRequiredMixin, UserPassesTestMixin
from django.db import transaction
from django.db.models import Count, Q
from django.http import (
    HttpRequest,
    HttpResponse,
    HttpResponseBadRequest,
    JsonResponse,
)
from django.shortcuts import get_object_or_404, redirect, render
from django.template.loader import render_to_string
from django.urls import reverse, reverse_lazy
from django.views import View
from django.views.generic import TemplateView

from django.core.exceptions import PermissionDenied

from clinic.models import Regimen

from . import explain, forms, models, optim
from .pharmaco import registry as pharmaco_registry
from .permissions import is_editor


class EditorRequiredMixin(UserPassesTestMixin):
    """Mixin ensuring the current user can manage simulator data."""

    raise_exception = True

    def test_func(self) -> bool:
        return is_editor(self.request.user)


class RegimenListView(LoginRequiredMixin, TemplateView):
    template_name = "simulator/regimen_list.html"

    def get_queryset(self):
        queryset = Regimen.objects.all()
        query = self.request.GET.get("q", "").strip()
        if query:
            queryset = queryset.filter(
                Q(name__icontains=query) | Q(components__icontains=query) | Q(notes__icontains=query)
            )
        line = self.request.GET.get("line", "").strip()
        if line:
            queryset = queryset.filter(line__iexact=line)
        intent = self.request.GET.get("intent", "").strip()
        if intent:
            queryset = queryset.filter(intent__iexact=intent)
        return queryset.order_by("name")

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        queryset = self.get_queryset()
        context.update(
            {
                "regimens": queryset,
                "filter_query": self.request.GET.get("q", "").strip(),
                "filter_line": self.request.GET.get("line", "").strip(),
                "filter_intent": self.request.GET.get("intent", "").strip(),
                "line_choices": _distinct_values(Regimen.objects.exclude(line=""), "line"),
                "intent_choices": _distinct_values(Regimen.objects.exclude(intent=""), "intent"),
                "is_editor": is_editor(self.request.user),
            }
        )
        return context


class RegimenCreateView(LoginRequiredMixin, EditorRequiredMixin, TemplateView):
    template_name = "simulator/regimen_form.html"
    form_class = forms.RegimenForm

    def get(self, request: HttpRequest, *args, **kwargs) -> HttpResponse:
        form = self.form_class()
        return render(request, self.template_name, {"form": form, "is_create": True})

    def post(self, request: HttpRequest, *args, **kwargs) -> HttpResponse:
        form = self.form_class(request.POST)
        if form.is_valid():
            regimen = form.save()
            messages.success(request, f"Regimen “{regimen.name}” created.")
            return redirect(reverse("simulator:regimen_list"))
        return render(request, self.template_name, {"form": form, "is_create": True})


class RegimenUpdateView(LoginRequiredMixin, EditorRequiredMixin, TemplateView):
    template_name = "simulator/regimen_form.html"
    form_class = forms.RegimenForm

    def dispatch(self, request: HttpRequest, *args, **kwargs):
        self.regimen = get_object_or_404(Regimen, pk=kwargs["pk"])
        return super().dispatch(request, *args, **kwargs)

    def get(self, request: HttpRequest, *args, **kwargs) -> HttpResponse:
        form = self.form_class(instance=self.regimen)
        return render(
            request,
            self.template_name,
            {"form": form, "regimen": self.regimen, "is_create": False},
        )

    def post(self, request: HttpRequest, *args, **kwargs) -> HttpResponse:
        form = self.form_class(request.POST, instance=self.regimen)
        if form.is_valid():
            regimen = form.save()
            messages.success(request, f"Regimen “{regimen.name}” updated.")
            return redirect(reverse("simulator:regimen_list"))
        return render(
            request,
            self.template_name,
            {"form": form, "regimen": self.regimen, "is_create": False},
        )


class RegimenDeleteView(LoginRequiredMixin, EditorRequiredMixin, View):
    template_name = "simulator/regimen_confirm_delete.html"
    success_url = reverse_lazy("simulator:regimen_list")

    def dispatch(self, request: HttpRequest, *args, **kwargs):
        self.regimen = get_object_or_404(Regimen, pk=kwargs["pk"])
        return super().dispatch(request, *args, **kwargs)

    def get(self, request: HttpRequest, *args, **kwargs) -> HttpResponse:
        return render(request, self.template_name, {"regimen": self.regimen})

    def post(self, request: HttpRequest, *args, **kwargs) -> HttpResponse:
        name = self.regimen.name
        self.regimen.delete()
        messages.success(request, f"Regimen “{name}” deleted.")
        if request.headers.get("HX-Request"):
            response = HttpResponse()
            response["HX-Redirect"] = str(self.success_url)
            return response
        return redirect(self.success_url)


class ScenarioManageListView(LoginRequiredMixin, EditorRequiredMixin, TemplateView):
    template_name = "simulator/scenario_manage_list.html"

    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        scenarios = (
            models.Scenario.objects.all()
            .prefetch_related("recommended_regimens")
            .annotate(attempt_count=Count("attempts"))
            .order_by("title")
        )
        context.update(
            {
                "scenarios": scenarios,
                "total_attempts": models.SimulationAttempt.objects.count(),
            }
        )
        return context


class ScenarioCreateView(LoginRequiredMixin, EditorRequiredMixin, TemplateView):
    template_name = "simulator/scenario_form.html"
    form_class = forms.ScenarioForm

    def get(self, request: HttpRequest, *args, **kwargs) -> HttpResponse:
        form = self.form_class()
        return render(request, self.template_name, {"form": form, "is_create": True})

    def post(self, request: HttpRequest, *args, **kwargs) -> HttpResponse:
        form = self.form_class(request.POST)
        if form.is_valid():
            scenario = form.save()
            messages.success(request, f"Scenario “{scenario.title}” created.")
            return redirect(reverse("simulator:scenario_manage_list"))
        return render(request, self.template_name, {"form": form, "is_create": True})


class ScenarioUpdateView(LoginRequiredMixin, EditorRequiredMixin, TemplateView):
    template_name = "simulator/scenario_form.html"
    form_class = forms.ScenarioForm

    def dispatch(self, request: HttpRequest, *args, **kwargs):
        self.scenario = get_object_or_404(models.Scenario, pk=kwargs["pk"])
        return super().dispatch(request, *args, **kwargs)

    def get(self, request: HttpRequest, *args, **kwargs) -> HttpResponse:
        form = self.form_class(instance=self.scenario)
        return render(
            request,
            self.template_name,
            {"form": form, "scenario": self.scenario, "is_create": False},
        )

    def post(self, request: HttpRequest, *args, **kwargs) -> HttpResponse:
        form = self.form_class(request.POST, instance=self.scenario)
        if form.is_valid():
            scenario = form.save()
            messages.success(request, f"Scenario “{scenario.title}” updated.")
            return redirect(reverse("simulator:scenario_manage_list"))
        return render(
            request,
            self.template_name,
            {"form": form, "scenario": self.scenario, "is_create": False},
        )


@login_required
def duplicate_scenario(request: HttpRequest, pk: int) -> HttpResponse:
    if not is_editor(request.user):
        raise PermissionDenied
    if request.method != "POST":
        return HttpResponseBadRequest("Invalid request method.")
    scenario = get_object_or_404(models.Scenario, pk=pk)
    with transaction.atomic():
        duplicate = models.Scenario.objects.create(
            title=f"{scenario.title} (copy)",
            clinical_stage=scenario.clinical_stage,
            summary=scenario.summary,
            risk_stratification=scenario.risk_stratification,
            lab_snapshot=scenario.lab_snapshot,
            guideline_notes=scenario.guideline_notes,
            expected_response=scenario.expected_response,
            active=scenario.active,
        )
        duplicate.recommended_regimens.set(scenario.recommended_regimens.all())
    messages.success(request, f"Scenario duplicated as “{duplicate.title}”.")
    return redirect(reverse("simulator:scenario_edit", args=[duplicate.pk]))


@login_required
def add_regimen_to_scenario(request: HttpRequest, pk: int) -> HttpResponse:
    if not is_editor(request.user):
        raise PermissionDenied
    if request.method != "POST":
        return HttpResponseBadRequest("Invalid request method.")
    scenario = get_object_or_404(models.Scenario, pk=pk)
    regimen_ids = request.POST.getlist("regimen_id")
    if not regimen_ids:
        return JsonResponse({"error": "Select a regimen to add."}, status=400)
    regimens = Regimen.objects.filter(pk__in=regimen_ids)
    if not regimens:
        return JsonResponse({"error": "Select a valid regimen."}, status=400)
    scenario.recommended_regimens.add(*regimens)
    return _render_guideline_update(request, scenario)


@login_required
def remove_regimen_from_scenario(request: HttpRequest, pk: int) -> HttpResponse:
    if not is_editor(request.user):
        raise PermissionDenied
    if request.method != "POST":
        return HttpResponseBadRequest("Invalid request method.")
    scenario = get_object_or_404(models.Scenario, pk=pk)
    regimen_id = request.POST.get("regimen_id")
    try:
        regimen = scenario.recommended_regimens.get(pk=int(regimen_id))
    except (Regimen.DoesNotExist, TypeError, ValueError):
        return JsonResponse({"error": "Regimen not linked to this scenario."}, status=400)
    scenario.recommended_regimens.remove(regimen)
    return _render_guideline_update(request, scenario)


def _render_guideline_update(request: HttpRequest, scenario: models.Scenario) -> HttpResponse:
    regimens = scenario.recommended_regimens.order_by("name")
    available = Regimen.objects.exclude(pk__in=regimens.values_list("pk", flat=True)).order_by("name")
    html = render_to_string(
        "simulator/_guideline_regimens_update.html",
        {
            "scenario": scenario,
            "regimens": regimens,
            "available_regimens": available,
        },
        request=request,
    )
    return HttpResponse(html)


def _distinct_values(queryset: Iterable[Regimen], field: str) -> list[str]:
    values = (
        queryset.order_by(field)
        .values_list(field, flat=True)
        .distinct()
    )
    return [value for value in values if value]


@login_required
def simulate_scenario(request: HttpRequest, pk: int) -> HttpResponse:
    if request.method != "POST":
        return HttpResponseBadRequest("Invalid request method.")
    scenario = get_object_or_404(models.Scenario, pk=pk, active=True)
    form = forms.SimulationParameterForm(request.POST)

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
    if not form.is_valid():
        html = render_to_string(
            "simulator/_simulation_form.html",
            {
                "form": form,
                "scenario": scenario,
                "warnings": getattr(form, "warnings", []),
                "drug_profiles": drug_profiles,
                "helptext_it": forms.SIMULATION_FORM_HELP_TEXT_IT,
                "helptext_en": forms.SIMULATION_FORM_HELP_TEXT_EN,
                "is_editor": is_editor(request.user),
            },
            request=request,
        )
        response = HttpResponse(html, status=400)
        response["HX-Retarget"] = "#simulation-form-container"
        response["HX-Reswap"] = "outerHTML"
        return response

    parameter_warnings = list(getattr(form, "warnings", []))
    attempt = models.SimulationAttempt.objects.create(
        scenario=scenario,
        user=request.user,
        parameters=form.cleaned_data,
        cohort_size=form.cleaned_data.get("cohort_size", 1),
        seed=form.cleaned_data.get("seed"),
    )
    try:
        attempt.run_model()
    except Exception as exc:  # pragma: no cover - defensive guard for solver failures
        attempt.delete()
        return JsonResponse({"error": f"Simulation solver failed: {exc}"}, status=400)

    summary_warnings: list[str] = []
    healthy_loss = attempt.results_summary.get("healthy_loss")
    tumor_reduction = attempt.results_summary.get("tumor_reduction")
    if healthy_loss is not None:
        if healthy_loss > 0.3:
            summary_warnings.append(
                "Healthy cell loss exceeded 30%, indicating a potentially toxic regimen—consider reducing doses."
            )
        elif healthy_loss > 0.2:
            summary_warnings.append("Healthy cell loss is above 20%; monitor for immunosuppression.")
    if tumor_reduction is not None and tumor_reduction < 0:
        summary_warnings.append(
            "Tumor reduction is negative (tumor growth). Adjust doses or shorten the horizon to explore alternatives."
        )

    context = {
        "scenario": scenario,
        "attempt": attempt,
        "summary": attempt.results_summary,
        "results": attempt.results,
        "warnings": parameter_warnings + summary_warnings,
        "kpi": explain.KPI,
    }
    html = render_to_string("simulator/_simulation_results.html", context, request=request)
    return HttpResponse(html)


# ────────────────────────────────────────────────────────────────────────────────
# Experiments / Optimization Lab Views
# ────────────────────────────────────────────────────────────────────────────────


class ExperimentsView(LoginRequiredMixin, EditorRequiredMixin, TemplateView):
    """
    Research-grade optimization lab.
    
    Provides interface for multi-objective optimization experiments:
    - Select scenario
    - Configure trials count
    - Run Optuna MOTPE sampler
    - View Pareto frontier results
    """
    template_name = "simulator/experiments.html"
    
    def get_context_data(self, **kwargs):
        context = super().get_context_data(**kwargs)
        context["scenarios"] = models.Scenario.objects.all().order_by("title")
        return context


@login_required
def run_experiment(request: HttpRequest, pk: int) -> HttpResponse:
    """
    Run multi-objective optimization experiment.
    
    POST params:
        - n_trials: int (default 100) - number of Optuna trials
    
    Returns:
        HTMX partial: _pareto_panel.html with Pareto frontier table
    """
    if not is_editor(request.user):
        raise PermissionDenied("Editor access required")
    
    if request.method != "POST":
        return HttpResponseBadRequest("POST required")
    
    scenario = get_object_or_404(models.Scenario, pk=pk)
    
    # Extract trial count from form
    try:
        n_trials = int(request.POST.get("n_trials", 100))
        if n_trials < 10:
            n_trials = 10
        elif n_trials > 500:
            n_trials = 500
    except (ValueError, TypeError):
        n_trials = 100
    
    # Run optimization
    result = optim.run_study(
        scenario=scenario,
        user_id=request.user.id,
        n_trials=n_trials,
        seed=42
    )
    
    # Prepare context for Pareto panel
    context = {
        "scenario": scenario,
        "result": result,
        "pareto_solutions": result.get("pareto", []),
    }
    
    html = render_to_string("simulator/_pareto_panel.html", context, request=request)
    return HttpResponse(html)
