from __future__ import annotations

from django.contrib import messages
from django.contrib.auth.decorators import login_required
from django.shortcuts import get_object_or_404, redirect, render
from django.urls import reverse

from . import forms, models


def scenario_list(request):
    scenarios = models.Scenario.objects.filter(active=True).prefetch_related("recommended_regimens")
    context = {"scenarios": scenarios}
    return render(request, "simulator/scenario_list.html", context)


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

    recommended_regimens = scenario.recommended_regimens.all()
    regimen_names = ", ".join(regimen.name for regimen in recommended_regimens) if recommended_regimens else "No guideline regimen linked yet."
    context = {
        "scenario": scenario,
        "attempts": attempts,
        "form": form,
        "recommended_regimens": recommended_regimens,
        "regimen_names": regimen_names,
    }
    return render(request, "simulator/scenario_detail.html", context)
