"""Celery configuration for mmportal."""

from __future__ import annotations

import os

from celery import Celery

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "mmportal.settings")

app = Celery("mmportal")
app.config_from_object("django.conf:settings", namespace="CELERY")
app.autodiscover_tasks()


@app.task(bind=True)
def debug_task(self):  # pragma: no cover - helper task
    print(f"Request: {self.request!r}")

