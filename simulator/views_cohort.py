"""
Views for virtual cohort simulations (in-silico clinical trials).
"""
from __future__ import annotations

import csv
import io
from typing import Any, Dict

from django.contrib.auth.decorators import login_required
from django.contrib.auth.mixins import LoginRequiredMixin
from django.core.cache import cache
from django.core.exceptions import PermissionDenied
from django.http import HttpRequest, HttpResponse, HttpResponseBadRequest
from django.shortcuts import get_object_or_404, render
from django.template.loader import render_to_string
from django.views.generic import TemplateView

from . import cohort, forms, models
from .permissions import is_editor


class CohortView(LoginRequiredMixin, TemplateView):
    """
    Virtual cohort runner interface.
    
    Allows users to configure and run in-silico clinical trials
    with synthetic patient populations.
    """
    template_name = "simulator/cohort.html"
    
    def get_context_data(self, **kwargs) -> Dict[str, Any]:
        context = super().get_context_data(**kwargs)
        context["scenarios"] = models.Scenario.objects.filter(active=True).order_by("title")
        return context


@login_required
def cohort_run(request: HttpRequest) -> HttpResponse:
    """
    Run cohort simulation (HTMX endpoint).
    
    POST params:
        - scenario_id: int
        - n_patients: int (10-1000)
        - seed: int
        - lenalidomide_dose, bortezomib_dose, daratumumab_dose: float
        - time_horizon: int
        - interaction_strength: float
        - len_on_days, bor_weekly, dara_interval: schedule params
        
    Returns:
        HTMX partial: _cohort_results.html with aggregate stats and plots
    """
    if request.method != "POST":
        return HttpResponseBadRequest("POST required")
    
    # Validate inputs
    try:
        scenario_id = int(request.POST.get("scenario_id", 0))
        n_patients = int(request.POST.get("n_patients", 100))
        seed = int(request.POST.get("seed", 42))
        
        # Validate ranges
        if not (10 <= n_patients <= 1000):
            return HttpResponseBadRequest("n_patients must be between 10 and 1000")
        
        scenario = get_object_or_404(models.Scenario, pk=scenario_id, active=True)
        
        # Extract regimen parameters
        regimen_params = {
            "lenalidomide_dose": float(request.POST.get("lenalidomide_dose", 25.0)),
            "bortezomib_dose": float(request.POST.get("bortezomib_dose", 1.3)),
            "daratumumab_dose": float(request.POST.get("daratumumab_dose", 16.0)),
            "time_horizon": int(request.POST.get("time_horizon", 168)),
            "interaction_strength": float(request.POST.get("interaction_strength", 0.1)),
            "len_on_days": int(request.POST.get("len_on_days", 21)),
            "bor_weekly": int(request.POST.get("bor_weekly", 0)),
            "dara_interval": int(request.POST.get("dara_interval", 7)),
        }
        
        # Validate regimen param ranges (same as SimulationParameterForm)
        if not (0 <= regimen_params["lenalidomide_dose"] <= 40):
            return HttpResponseBadRequest("lenalidomide_dose must be 0-40 mg")
        if not (0 <= regimen_params["bortezomib_dose"] <= 1.6):
            return HttpResponseBadRequest("bortezomib_dose must be 0-1.6 mg/m²")
        if not (0 <= regimen_params["daratumumab_dose"] <= 16):
            return HttpResponseBadRequest("daratumumab_dose must be 0-16 mg/kg")
        if not (56 <= regimen_params["time_horizon"] <= 365):
            return HttpResponseBadRequest("time_horizon must be 56-365 days")
        if not (0 <= regimen_params["interaction_strength"] <= 0.2):
            return HttpResponseBadRequest("interaction_strength must be 0-0.2")
        
    except (ValueError, TypeError) as e:
        return HttpResponseBadRequest(f"Invalid input: {e}")
    
    # Run cohort simulation
    result = cohort.run_cohort(
        scenario=scenario,
        n=n_patients,
        regimen_params=regimen_params,
        user_id=request.user.id,
        seed=seed,
    )
    
    # Cache result for CSV export (30 min TTL)
    cache.set(f"cohort_result_{result['cohort_id']}", result, 1800)
    
    # Render results partial
    context = {
        "scenario": scenario,
        "result": result,
        "n_patients": n_patients,
        "regimen_params": regimen_params,
    }
    
    html = render_to_string("simulator/_cohort_results.html", context, request=request)
    return HttpResponse(html)


@login_required
def cohort_export(request: HttpRequest, cohort_id: str) -> HttpResponse:
    """
    Export cohort results as CSV.
    
    Retrieves cached cohort result and generates CSV with:
    - One row per patient (patient_id, efficacy, toxicity, recurrence, auc)
    - Final summary row with aggregates
    """
    result = cache.get(f"cohort_result_{cohort_id}")
    
    if not result:
        return HttpResponseBadRequest("Cohort result not found or expired. Please re-run the cohort.")
    
    # Generate CSV
    output = io.StringIO()
    writer = csv.writer(output)
    
    # Header
    writer.writerow([
        "patient_id",
        "tumor_reduction",
        "healthy_loss",
        "time_to_recurrence",
        "auc_lenalidomide",
        "auc_bortezomib",
        "auc_daratumumab",
        "auc_total",
        "effective",
    ])
    
    # Patient rows
    for summary in result["summaries"]:
        auc = summary.get("auc", {})
        auc_total = sum([
            auc.get("lenalidomide", 0.0),
            auc.get("bortezomib", 0.0),
            auc.get("daratumumab", 0.0),
        ])
        
        writer.writerow([
            summary.get("patient_id", ""),
            f"{summary.get('tumor_reduction', 0.0):.6f}",
            f"{summary.get('healthy_loss', 0.0):.6f}",
            summary.get("time_to_recurrence", ""),
            f"{auc.get('lenalidomide', 0.0):.6f}",
            f"{auc.get('bortezomib', 0.0):.6f}",
            f"{auc.get('daratumumab', 0.0):.6f}",
            f"{auc_total:.6f}",
            summary.get("effective", False),
        ])
    
    # Summary row
    agg = result["aggregates"]
    writer.writerow([])
    writer.writerow(["SUMMARY STATISTICS"])
    writer.writerow(["efficacy_mean", f"{agg['efficacy_mean']:.6f}"])
    writer.writerow(["efficacy_p05", f"{agg['efficacy_p05']:.6f}"])
    writer.writerow(["efficacy_p95", f"{agg['efficacy_p95']:.6f}"])
    writer.writerow(["toxicity_mean", f"{agg['toxicity_mean']:.6f}"])
    writer.writerow(["toxicity_p05", f"{agg['toxicity_p05']:.6f}"])
    writer.writerow(["toxicity_p95", f"{agg['toxicity_p95']:.6f}"])
    writer.writerow(["recurrence_rate", f"{agg['recurrence_rate']:.6f}"])
    writer.writerow(["auc_mean", f"{agg['auc_mean']:.6f}"])
    writer.writerow(["auc_p95", f"{agg['auc_p95']:.6f}"])
    writer.writerow(["seed", result["seed"]])
    writer.writerow(["n_patients", result["n"]])
    
    # Build response
    response = HttpResponse(output.getvalue(), content_type="text/csv")
    response["Content-Disposition"] = f'attachment; filename="cohort_{cohort_id}.csv"'
    
    return response
