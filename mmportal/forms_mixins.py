from __future__ import annotations


class BootstrapValidationMixin:
    """Mixin that tags invalid fields with Bootstrap-compatible error classes."""

    error_css_class = "is-invalid"

    def full_clean(self) -> None:  # type: ignore[override]
        super().full_clean()  # type: ignore[misc]
        if not getattr(self, "is_bound", False):
            return
        for name in getattr(self, "errors", {}):
            if name in getattr(self, "fields", {}):
                widget = self.fields[name].widget  # type: ignore[index]
                existing = widget.attrs.get("class", "")
                if self.error_css_class not in existing.split():
                    widget.attrs["class"] = (f"{existing} {self.error_css_class}").strip()
