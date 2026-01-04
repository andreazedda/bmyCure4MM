from __future__ import annotations

from pathlib import Path

from django.conf import settings
from django.core.files.storage import default_storage
from django.db import models


class ChemJob(models.Model):
    """Tracks executions of cheminformatics utilities."""

    PARAM = "PARAM"
    BIND = "BIND"
    SIM = "SIM"

    KIND_CHOICES = [
        (PARAM, "Drug Parameters"),
        (BIND, "Binding Visualizer"),
        (SIM, "Similarity Search"),
    ]

    kind = models.CharField(max_length=10, choices=KIND_CHOICES)
    created = models.DateTimeField(auto_now_add=True)
    input_a = models.CharField(max_length=255, blank=True)
    input_b = models.CharField(max_length=255, blank=True)
    out_html = models.FileField(upload_to="chem/html/", blank=True, null=True)
    out_csv = models.FileField(upload_to="chem/csv/", blank=True, null=True)
    log = models.TextField(blank=True)
    user = models.ForeignKey(
        settings.AUTH_USER_MODEL,
        on_delete=models.SET_NULL,
        null=True,
        blank=True,
        related_name="chem_jobs",
    )
    
    # Progress tracking fields
    progress_percent = models.IntegerField(default=0)
    progress_message = models.CharField(max_length=255, blank=True)
    
    # API preferences (stored as JSON)
    api_preferences = models.JSONField(default=dict, blank=True)

    class Meta:
        ordering = ["-created"]

    def __str__(self) -> str:
        return f"{self.get_kind_display()} job @ {self.created:%Y-%m-%d %H:%M}"

    @property
    def media_directory(self) -> Path:
        return Path(settings.MEDIA_ROOT) / "chem" / str(self.pk)

    def thumbnail_url(self) -> str | None:
        relative = f"chem/{self.pk}/thumb.png"
        if default_storage.exists(relative):
            return settings.MEDIA_URL + relative
        return None

    def status_label(self) -> tuple[str, str]:
        log_text = (self.log or "").strip()
        if log_text.startswith("ERROR"):
            return "Failed", "danger"
        if self.out_html or self.out_csv:
            return "Completed", "success"
        return "Queued", "secondary"
    
    def update_progress(self, percent: int, message: str = "") -> None:
        """Update job progress."""
        self.progress_percent = min(100, max(0, percent))
        self.progress_message = message
        self.save(update_fields=['progress_percent', 'progress_message'])
