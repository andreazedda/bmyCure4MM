"""Transactional email presentation; no tracking, external assets or ticket contents."""

import re
from email.utils import formataddr, parseaddr
from urllib.parse import urlsplit

from django.conf import settings
from django.template.loader import render_to_string

from .services import public_base

# Subject, explanation and action. Locale is runtime configuration, independent of identity.
COPY = {
    "verify": (
        (
            "Conferma la tua email",
            "Conferma il tuo indirizzo per inviare la richiesta ad Andrea. Dopo la verifica riceverai l'esito della revisione: l'accesso richiede la sua approvazione.",
            "Conferma email",
        ),
        (
            "Confirm your email",
            "Confirm your address to submit your request to Andrea. Email verification is followed by review; access requires his approval.",
            "Confirm email",
        ),
    ),
    "owner": (
        (
            "Imposta la password del proprietario",
            "Il tuo account proprietario è pronto. Scegli la password personale della piattaforma; al primo accesso potrai configurare l'autenticatore e gestire richieste e ticket.",
            "Imposta la mia password",
        ),
        (
            "Set your owner password",
            "Your owner account is ready. Choose your personal platform password, then sign in to set up your authenticator and manage requests and tickets.",
            "Set my password",
        ),
    ),
    "activate": (
        (
            "La tua richiesta è stata approvata",
            "Andrea ha approvato il tuo accesso. Imposta la password e accedi alla guida iniziale per esplorare la piattaforma e inviare suggerimenti.",
            "Imposta la mia password",
        ),
        (
            "Your access request is approved",
            "Andrea has approved your access. Set your password and follow the getting-started guide to explore the platform and send suggestions.",
            "Set my password",
        ),
    ),
    "reset": (
        (
            "Reimposta la tua password",
            "Hai richiesto una nuova password per bmyCure4MM. Usa il pulsante per sceglierla. Lo stato e la scadenza del tuo accesso restano quelli assegnati dal proprietario.",
            "Scegli una nuova password",
        ),
        (
            "Reset your password",
            "Use the button to choose a new bmyCure4MM password. Your approval status and access expiry remain governed by the owner.",
            "Choose a new password",
        ),
    ),
    "verified_request": (
        (
            "Una richiesta verificata attende la tua revisione",
            "Una persona ha confermato la propria email. Puoi leggere il motivo della richiesta e decidere se approvare l'accesso dalla console.",
            "Esamina la richiesta",
        ),
        (
            "A verified request needs your review",
            "An applicant has confirmed their email. Read the request and approve or reject access in the console.",
            "Review requests",
        ),
    ),
    "rejected": (
        (
            "Aggiornamento sulla richiesta di accesso",
            "La richiesta non è stata approvata. Per chiarimenti sull'accesso puoi rispondere a questa email.",
            "",
        ),
        (
            "An update on your access request",
            "Your request was not approved. You can reply to this email if you need clarification.",
            "",
        ),
    ),
    "feedback": (
        (
            "Hai un nuovo ticket da leggere",
            "Un utente ha inviato un suggerimento o una richiesta di supporto. Leggi il messaggio e rispondi dalla console: il contenuto resta privato.",
            "Leggi i ticket",
        ),
        (
            "A new ticket needs your attention",
            "A user has sent a suggestion or support request. Read and reply in the console; ticket contents remain private.",
            "Read tickets",
        ),
    ),
    "feedback_ack": (
        (
            "Il tuo ticket è stato ricevuto",
            "Il messaggio è salvato nella piattaforma. Puoi seguirne lo stato e leggere la risposta di Andrea nella pagina dei tuoi ticket.",
            "Apri i miei ticket",
        ),
        (
            "We received your ticket",
            "Your message is saved. Follow its status and read Andrea's response on your tickets page.",
            "Open my tickets",
        ),
    ),
    "feedback_reply": (
        (
            "Hai una risposta al tuo ticket",
            "Andrea ha pubblicato una risposta al tuo messaggio. Accedi alla piattaforma per leggerla nella pagina dei tuoi ticket.",
            "Leggi la risposta",
        ),
        (
            "Your ticket has a reply",
            "Andrea has replied to your message. Sign in to read the response on your tickets page.",
            "Read the reply",
        ),
    ),
    "expiry": (
        (
            "Il tuo accesso scade entro 24 ore",
            "Se desideri continuare a usare la piattaforma, contatta Andrea attraverso un ticket prima della scadenza.",
            "Apri i miei ticket",
        ),
        (
            "Your access expires within 24 hours",
            "If you would like to continue using the platform, contact Andrea through a ticket before your access expires.",
            "Open my tickets",
        ),
    ),
    "renewed": (
        (
            "Il tuo accesso è stato rinnovato",
            "Andrea ha rinnovato il tuo accesso. Accedi di nuovo alla piattaforma per continuare.",
            "Accedi alla piattaforma",
        ),
        (
            "Your access has been renewed",
            "Andrea has renewed your access. Sign in again to continue using the platform.",
            "Sign in",
        ),
    ),
    "suspended": (
        (
            "Il tuo accesso è sospeso",
            "L'accesso alla piattaforma è stato sospeso dal proprietario. Per chiarimenti puoi rispondere a questa email.",
            "",
        ),
        (
            "Your access is suspended",
            "The owner has suspended your access. You can reply to this email for clarification.",
            "",
        ),
    ),
    "revoked": (
        (
            "Il tuo accesso è stato revocato",
            "L'accesso alla piattaforma è stato revocato dal proprietario. Per chiarimenti puoi rispondere a questa email.",
            "",
        ),
        (
            "Your access has been revoked",
            "The owner has revoked your access. You can reply to this email for clarification.",
            "",
        ),
    ),
}


def presentation(item, payload, capability=None):
    italian = settings.ACCOUNT_EMAIL_LANGUAGE == "it"
    locale = "it" if italian else "en"
    fallback = (
        ("Aggiornamento sul tuo account", payload["body"], ""),
        ("Account update", payload["body"], ""),
    )
    title, explanation, action_label = COPY.get(item.kind, fallback)[0 if italian else 1]
    base = public_base()
    action_url = ""
    if capability:
        match = re.search(r"https?://[^\s]+#token=[A-Za-z0-9_-]+", payload["body"])
        if match and urlsplit(match.group()).netloc == urlsplit(base).netloc:
            action_url = match.group()
    elif action_label:
        path = {
            "verified_request": "/manage/access-requests/",
            "feedback": "/manage/feedback/",
            "renewed": "/accounts/login/",
        }.get(item.kind, "/feedback/")
        action_url = base + path
    context = {
        "locale": locale,
        "italian": italian,
        "title": title,
        "explanation": explanation,
        "nickname": item.identity.nickname,
        "action_label": action_label,
        "action_url": action_url,
        "expires_at": capability.expires_at if capability else None,
        "access_expires_at": item.identity.expires_at
        if item.kind in {"activate", "renewed", "expiry"}
        else None,
        "origin": base,
    }
    return (
        title + " · bmyCure4MM",
        render_to_string("accounts/email/notification.txt", context),
        render_to_string("accounts/email/notification.html", context),
        formataddr(("bmyCure4MM", parseaddr(settings.DEFAULT_FROM_EMAIL)[1])),
    )
