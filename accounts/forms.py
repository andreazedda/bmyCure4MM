from django import forms
from django.core.exceptions import ValidationError
from django.utils.html import format_html

from .models import SupportTicket
from .services import normalize_email


class StrictForm(forms.Form):
    """Reject privilege/scope injection instead of silently accepting extra keys."""

    extra_fields = {"csrfmiddlewaretoken"}

    def __init__(self, *args, **kwargs):
        kwargs.setdefault("label_suffix", "")
        super().__init__(*args, **kwargs)
        for name, field in self.fields.items():
            if field.label is None:
                field.label = name.replace("_", " ").capitalize()
        labels_it = {
            "Email": "Email",
            "Nickname": "Nickname",
            "Reason": "Motivo della richiesta",
            "Professional role": "Ruolo professionale (facoltativo)",
            "Affiliation": "Affiliazione (facoltativa)",
            "Password": "Password",
            "Password confirm": "Conferma password",
            "Code": "Codice autenticatore",
            "I have read the privacy notice and research-only limits. I will not enter real patient data.": "Ho letto l’informativa e i limiti di ricerca. Non inserirò dati reali di pazienti.",
        }
        for field in self.fields.values():
            if field.label in labels_it and field.label != labels_it[field.label]:
                field.widget.attrs.update(
                    {
                        "aria-label": field.label,
                        "data-label-en": field.label,
                        "data-label-it": labels_it[field.label],
                    }
                )
                field.label = format_html(
                    '<span class="t-en">{}</span><span class="t-it">{}</span>',
                    field.label,
                    labels_it[field.label],
                )
            field.widget.attrs.setdefault(
                "class",
                "form-check-input" if isinstance(field.widget, forms.CheckboxInput) else "form-control",
            )

    def clean(self):
        cleaned = super().clean()
        if self.is_bound and set(self.data) - set(self.fields) - self.extra_fields:
            raise ValidationError("Unexpected form fields.")
        return cleaned


class RequestForm(StrictForm):
    email = forms.EmailField(max_length=254)
    nickname = forms.CharField(max_length=60)
    reason = forms.CharField(max_length=1000, widget=forms.Textarea(attrs={"rows": 3}))
    professional_role = forms.ChoiceField(
        required=False,
        choices=[
            ("", "Optional"),
            ("researcher", "Researcher"),
            ("clinician", "Clinician"),
            ("engineer", "Engineer"),
            ("student", "Student"),
            ("other", "Other"),
        ],
    )
    affiliation = forms.CharField(max_length=120, required=False)
    acknowledgement = forms.BooleanField(
        label="I have read the privacy notice and research-only limits. I will not enter real patient data."
    )

    def clean_email(self):
        return normalize_email(self.cleaned_data["email"])


class LoginForm(StrictForm):
    email = forms.EmailField(max_length=254)
    password = forms.CharField(widget=forms.PasswordInput(attrs={"autocomplete": "current-password"}))


class EmailForm(StrictForm):
    email = forms.EmailField(max_length=254)


class CapabilityForm(StrictForm):
    token = forms.CharField(max_length=100, widget=forms.HiddenInput)


class ActivationForm(CapabilityForm):
    password = forms.CharField(
        widget=forms.PasswordInput(attrs={"autocomplete": "new-password"}), min_length=12
    )
    password_confirm = forms.CharField(widget=forms.PasswordInput(attrs={"autocomplete": "new-password"}))

    def clean(self):
        cleaned = super().clean()
        if cleaned.get("password") != cleaned.get("password_confirm"):
            raise ValidationError("The passwords do not match.")
        return cleaned


class DecisionForm(StrictForm):
    identity_id = forms.UUIDField(widget=forms.HiddenInput)
    action = forms.ChoiceField(
        choices=[(v, v) for v in ("approve", "reject", "resend", "suspend", "revoke", "renew")]
    )
    days = forms.IntegerField(min_value=1, max_value=14, initial=14, required=False)


class SupportTicketForm(StrictForm):
    submission_key = forms.UUIDField(widget=forms.HiddenInput)
    category = forms.ChoiceField(choices=SupportTicket.CATEGORIES)
    message = forms.CharField(
        max_length=2000,
        widget=forms.Textarea(attrs={"rows": 5}),
        help_text="Do not include real patient data, passwords, tokens or identifying medical details.",
    )


class ReplyForm(StrictForm):
    ticket_id = forms.UUIDField(widget=forms.HiddenInput)
    category = forms.ChoiceField(choices=SupportTicket.CATEGORIES)
    reply = forms.CharField(max_length=2000, required=False, widget=forms.Textarea(attrs={"rows": 3}))
    internal_note = forms.CharField(max_length=1000, required=False, widget=forms.Textarea(attrs={"rows": 2}))
    status = forms.ChoiceField(choices=SupportTicket._meta.get_field("status").choices)


class SwitchForm(StrictForm):
    requests_open = forms.BooleanField(required=False)
    activations_open = forms.BooleanField(required=False)


class RetryForm(StrictForm):
    notification_id = forms.UUIDField()


class OTPForm(StrictForm):
    code = forms.RegexField(
        r"^[0-9]{6}$",
        max_length=6,
        widget=forms.TextInput(attrs={"autocomplete": "one-time-code", "inputmode": "numeric"}),
    )
