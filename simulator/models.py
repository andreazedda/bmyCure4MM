from __future__ import annotations

from django.conf import settings
from django.db import models

from clinic.models import Assessment, Regimen


class Scenario(models.Model):
    """Clinical scenarios used to practice MM treatment decision making."""

    STAGE_CHOICES = [
        ("newly_diagnosed", "Newly Diagnosed"),
        ("relapsed_refractory", "Relapsed/Refractory"),
        ("maintenance", "Maintenance"),
        ("supportive", "Supportive Care"),
    ]

    title = models.CharField(max_length=200)
    clinical_stage = models.CharField(
        max_length=32,
        choices=STAGE_CHOICES,
        default="newly_diagnosed",
    )
    summary = models.TextField(
        help_text="Concise description of the case including age, stage and key features."
    )
    risk_stratification = models.CharField(
        max_length=128,
        blank=True,
        help_text="e.g., High-risk cytogenetics, R-ISS II.",
    )
    lab_snapshot = models.JSONField(
        blank=True,
        default=dict,
        help_text="Dictionary of key laboratory values shown to the learner.",
    )
    recommended_regimens = models.ManyToManyField(
        Regimen,
        related_name="training_scenarios",
        blank=True,
    )
    guideline_notes = models.TextField(
        blank=True,
        help_text="Rationale, supporting evidence or pearls displayed after completion.",
    )
    expected_response = models.CharField(
        max_length=4,
        choices=Assessment.RESPONSE_CHOICES,
        blank=True,
        help_text="Target IMWG response. Optional; used for feedback.",
    )
    created = models.DateTimeField(auto_now_add=True)
    updated = models.DateTimeField(auto_now=True)
    active = models.BooleanField(default=True)

    class Meta:
        ordering = ["title"]

    def __str__(self) -> str:
        return self.title

    def lab_items(self) -> list[tuple[str, str]]:
        snapshot = self.lab_snapshot or {}
        return list(snapshot.items())


class SimulationAttempt(models.Model):
    """Captures a clinician's decisions for a given scenario."""

    scenario = models.ForeignKey(
        Scenario,
        on_delete=models.CASCADE,
        related_name="attempts",
    )
    user = models.ForeignKey(
        settings.AUTH_USER_MODEL,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
    )
    selected_regimen = models.ForeignKey(
        Regimen,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
    )
    predicted_response = models.CharField(
        max_length=4,
        choices=Assessment.RESPONSE_CHOICES,
        blank=True,
    )
    confidence = models.PositiveSmallIntegerField(
        default=70,
        help_text="Confidence in chosen plan (0-100%).",
    )
    notes = models.TextField(blank=True)
    submitted = models.DateTimeField(auto_now_add=True)
    is_guideline_aligned = models.BooleanField(default=False)

    class Meta:
        ordering = ["-submitted"]

    def __str__(self) -> str:
        user_display = self.user.get_username() if self.user else "anonymous"
        return f"{self.scenario} attempt by {user_display} @ {self.submitted:%Y-%m-%d}"

