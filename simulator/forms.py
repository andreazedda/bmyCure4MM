from __future__ import annotations

from django import forms

from clinic.models import Assessment, Regimen

from . import models


class SimulationAttemptForm(forms.ModelForm):
    """Form used for clinician responses within a scenario."""

    selected_regimen = forms.ModelChoiceField(
        label="Chosen regimen",
        queryset=Regimen.objects.all(),
        required=False,
        help_text="Select the primary regimen you would start for this case.",
        widget=forms.Select(attrs={"class": "form-select"}),
    )
    predicted_response = forms.ChoiceField(
        label="Expected response",
        choices=[("", "Select response...")] + list(Assessment.RESPONSE_CHOICES),
        required=False,
        widget=forms.Select(attrs={"class": "form-select"}),
    )
    confidence = forms.IntegerField(
        min_value=0,
        max_value=100,
        initial=70,
        help_text="How confident are you in this plan? (%)",
        widget=forms.NumberInput(attrs={"class": "form-control", "min": "0", "max": "100"}),
    )
    notes = forms.CharField(
        widget=forms.Textarea(attrs={"rows": 4, "class": "form-control"}),
        required=False,
        label="Clinical reasoning",
    )

    class Meta:
        model = models.SimulationAttempt
        fields = ["selected_regimen", "predicted_response", "confidence", "notes"]
