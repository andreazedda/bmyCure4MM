from __future__ import annotations

from django import forms

from . import models


class AssessmentForm(forms.ModelForm):
    """Quick entry form for laboratory assessment."""

    date = forms.DateField(widget=forms.DateInput(attrs={"type": "date"}))

    class Meta:
        model = models.Assessment
        fields = [
            "date",
            "m_protein_g_dl",
            "kappa_mg_l",
            "lambda_mg_l",
            "flc_ratio",
            "hemoglobin_g_dl",
            "calcium_mg_dl",
            "creatinine_mg_dl",
            "beta2m_mg_l",
            "albumin_g_dl",
            "ldH_u_l",
            "r_iss",
            "response",
            "notes",
        ]

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for name, field in self.fields.items():
            widget = field.widget
            if isinstance(widget, forms.Select):
                widget.attrs.setdefault("class", "form-select")
            elif isinstance(widget, forms.CheckboxInput):
                widget.attrs.setdefault("class", "form-check-input")
            else:
                widget.attrs.setdefault("class", "form-control")

    def clean(self):
        cleaned = super().clean()
        non_negative_fields = [
            "m_protein_g_dl",
            "kappa_mg_l",
            "lambda_mg_l",
            "flc_ratio",
            "hemoglobin_g_dl",
            "calcium_mg_dl",
            "creatinine_mg_dl",
            "beta2m_mg_l",
            "albumin_g_dl",
            "ldH_u_l",
        ]
        for field in non_negative_fields:
            value = cleaned.get(field)
            if value is not None and value < 0:
                self.add_error(field, "Must be ≥ 0")
        return cleaned


class PatientTherapyForm(forms.ModelForm):
    """Manage therapy course entries."""

    start_date = forms.DateField(widget=forms.DateInput(attrs={"type": "date"}))
    end_date = forms.DateField(
        widget=forms.DateInput(attrs={"type": "date"}),
        required=False,
    )

    class Meta:
        model = models.PatientTherapy
        fields = [
            "regimen",
            "start_date",
            "end_date",
            "outcome",
            "adverse_events",
        ]
