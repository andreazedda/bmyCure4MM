from __future__ import annotations

from django.db import models
from django.utils import timezone


class Patient(models.Model):
    """Core patient demographics for MM management."""

    SEX_CHOICES = [
        ("M", "Male"),
        ("F", "Female"),
        ("O", "Other"),
    ]

    mrn = models.CharField("MRN", max_length=32, unique=True)
    first_name = models.CharField(max_length=64)
    last_name = models.CharField(max_length=64)
    birth_date = models.DateField()
    sex = models.CharField(max_length=1, choices=SEX_CHOICES)
    diagnosis_date = models.DateField()
    notes = models.TextField(blank=True)

    class Meta:
        ordering = ["last_name", "first_name"]

    def __str__(self) -> str:
        return f"{self.last_name}, {self.first_name} ({self.mrn})"

    @property
    def age(self) -> int:
        today = timezone.now().date()
        return today.year - self.birth_date.year - (
            (today.month, today.day) < (self.birth_date.month, self.birth_date.day)
        )


class CytogeneticAbnormality(models.Model):
    """Catálogo delle anomalie citogenetiche rilevanti."""

    code = models.CharField(max_length=32, unique=True)
    description = models.CharField(max_length=255)

    class Meta:
        ordering = ["code"]

    def __str__(self) -> str:
        return f"{self.code} - {self.description}"


class PatientCytogenetics(models.Model):
    """Review storica delle anomalie per paziente."""

    patient = models.ForeignKey(Patient, on_delete=models.CASCADE, related_name="cytogenetics")
    abnormality = models.ForeignKey(
        CytogeneticAbnormality,
        on_delete=models.CASCADE,
        related_name="patient_records",
    )
    detected_on = models.DateField()
    method = models.CharField(max_length=128, blank=True)

    class Meta:
        ordering = ["-detected_on"]
        unique_together = ("patient", "abnormality", "detected_on")

    def __str__(self) -> str:
        return f"{self.patient} - {self.abnormality.code}"


class Assessment(models.Model):
    """Key laboratory values per IMWG assessment."""

    RESPONSE_CHOICES = [
        ("sCR", "Stringent CR"),
        ("CR", "Complete Response"),
        ("VGPR", "Very Good Partial Response"),
        ("PR", "Partial Response"),
        ("SD", "Stable Disease"),
        ("PD", "Progressive Disease"),
    ]

    R_ISS_CHOICES = [
        ("I", "R-ISS I"),
        ("II", "R-ISS II"),
        ("III", "R-ISS III"),
    ]

    patient = models.ForeignKey(Patient, on_delete=models.CASCADE, related_name="assessments")
    date = models.DateField()
    m_protein_g_dl = models.DecimalField("M-Protein (g/dL)", max_digits=5, decimal_places=2, null=True, blank=True)
    kappa_mg_l = models.DecimalField("Kappa (mg/L)", max_digits=6, decimal_places=2, null=True, blank=True)
    lambda_mg_l = models.DecimalField("Lambda (mg/L)", max_digits=6, decimal_places=2, null=True, blank=True)
    flc_ratio = models.DecimalField("FLC Ratio", max_digits=6, decimal_places=2, null=True, blank=True)
    hemoglobin_g_dl = models.DecimalField("Hemoglobin (g/dL)", max_digits=4, decimal_places=1, null=True, blank=True)
    calcium_mg_dl = models.DecimalField("Calcium (mg/dL)", max_digits=4, decimal_places=1, null=True, blank=True)
    creatinine_mg_dl = models.DecimalField("Creatinine (mg/dL)", max_digits=4, decimal_places=2, null=True, blank=True)
    beta2m_mg_l = models.DecimalField("Beta-2 Microglobulin (mg/L)", max_digits=5, decimal_places=2, null=True, blank=True)
    albumin_g_dl = models.DecimalField("Albumin (g/dL)", max_digits=4, decimal_places=1, null=True, blank=True)
    ldH_u_l = models.DecimalField("LDH (U/L)", max_digits=6, decimal_places=1, null=True, blank=True)
    r_iss = models.CharField("R-ISS", max_length=4, choices=R_ISS_CHOICES, blank=True)
    response = models.CharField(max_length=4, choices=RESPONSE_CHOICES, blank=True)
    notes = models.TextField(blank=True)

    class Meta:
        ordering = ["-date"]
        verbose_name = "Assessment"
        verbose_name_plural = "Assessments"

    def __str__(self) -> str:
        return f"{self.patient} assessment on {self.date}"


class Regimen(models.Model):
    """Therapy regimen definitions."""

    name = models.CharField(max_length=128)
    line = models.CharField(max_length=64, help_text="Line of therapy (e.g. frontline, salvage).")
    components = models.TextField(help_text="Comma-separated agents or regimen details.")
    intent = models.CharField(max_length=128, blank=True)
    notes = models.TextField(blank=True)

    class Meta:
        ordering = ["name"]

    def __str__(self) -> str:
        return self.name


class PatientTherapy(models.Model):
    """Therapy timeline for each patient."""

    patient = models.ForeignKey(Patient, on_delete=models.CASCADE, related_name="therapies")
    regimen = models.ForeignKey(Regimen, on_delete=models.CASCADE, related_name="patient_courses")
    start_date = models.DateField()
    end_date = models.DateField(null=True, blank=True)
    outcome = models.CharField(max_length=128, blank=True)
    adverse_events = models.TextField(blank=True)

    class Meta:
        ordering = ["-start_date"]

    def __str__(self) -> str:
        return f"{self.patient} - {self.regimen.name}"
