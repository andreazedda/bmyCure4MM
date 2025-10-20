from __future__ import annotations


from django import forms

from mmportal.forms_mixins import BootstrapValidationMixin

from . import models


class PatientForm(BootstrapValidationMixin, forms.ModelForm):
    """Form for creating and editing patients."""

    birth_date = forms.DateField(
        widget=forms.DateInput(attrs={"type": "date"}),
        help_text="Patient's date of birth"
    )
    diagnosis_date = forms.DateField(
        widget=forms.DateInput(attrs={"type": "date"}),
        help_text="Date of multiple myeloma diagnosis"
    )

    class Meta:
        model = models.Patient
        fields = [
            "mrn",
            "first_name", 
            "last_name",
            "birth_date",
            "sex",
            "diagnosis_date",
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
            elif isinstance(widget, forms.Textarea):
                widget.attrs.setdefault("class", "form-control")
                widget.attrs.setdefault("rows", "3")
            else:
                widget.attrs.setdefault("class", "form-control")

    def clean_mrn(self):
        mrn = self.cleaned_data.get("mrn")
        if mrn:
            mrn = mrn.strip().upper()
            # Check for uniqueness if this is a new patient or MRN changed
            if self.instance.pk is None or self.instance.mrn != mrn:
                if models.Patient.objects.filter(mrn=mrn).exists():
                    raise forms.ValidationError("A patient with this MRN already exists.")
        return mrn


class AssessmentForm(BootstrapValidationMixin, forms.ModelForm):
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


class PatientTherapyForm(BootstrapValidationMixin, forms.ModelForm):
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
