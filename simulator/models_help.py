from __future__ import annotations

from django.db import models


class HelpArticle(models.Model):
    slug = models.SlugField(unique=True)
    title_en = models.CharField(max_length=200)
    body_en = models.TextField()
    title_it = models.CharField(max_length=200)
    body_it = models.TextField()
    updated = models.DateTimeField(auto_now=True)

    class Meta:
        ordering = ["slug"]

    def __str__(self) -> str:  # pragma: no cover - admin readability
        return self.slug

    def as_lang(self, lang: str) -> dict[str, str]:
        if lang == "it":
            return {"title": self.title_it, "body": self.body_it}
        return {"title": self.title_en, "body": self.body_en}
