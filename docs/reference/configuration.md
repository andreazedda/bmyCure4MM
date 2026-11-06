---
title: Configuration and Runtime Settings
language: en
status: CURRENT_PARTIAL
owner: Andrea Zedda
audience: developers and operators
last_verified_at: 2026-11-06
---

# Configuration (env, settings, runtime)

The tables below are a quick reference, not the complete runtime configuration. The [current account/SMTP contract](../operations/AUTHENTICATION_AND_ACCESS.md#runtime-configuration-and-email), [actual architecture](../architecture/CURRENT_SYSTEM.md) and [disposable handbook setup](../development/DEVELOPMENT.md) govern onboarding and operations.

`DJANGO_ALLOWED_HOSTS` and `DJANGO_CSRF_TRUSTED_ORIGINS` are the preferred environment names; legacy `ALLOWED_HOSTS` / `CSRF_TRUSTED_ORIGINS` are accepted fallbacks. The actual cluster uses SQLite on its data PVC. `BMY_SOURCE_SHA` identifies the image source and `CONTAINER_IMAGE_DIGEST` its immutable artifact; neither is the same as the local checkout HEAD. Never print the full environment to diagnose a mismatch.

`CELERY_TASK_ALWAYS_EAGER` is not parsed as a production environment switch. Test settings set it explicitly in Python; see [manual troubleshooting](../development/TROUBLESHOOTING.md).

=== "IT"
    ## Variabili d’ambiente principali

    ## Se cerchi…

    | Voglio… | Vai a… |
    | --- | --- |
    | avvio stack | `Guides → Operations` |
    | troubleshooting | `Guides → Troubleshooting` |
    | endpoint principali | `Reference → Endpoints` |

    ## Terminologia

    - **env var**: variabile d’ambiente letta da `mmportal/settings.py`.
    - **feature flag**: toggle che abilita/disabilita feature.

    ### Django

    - `DJANGO_DEBUG`: `1` (default) / `0`
    - `DJANGO_SECRET_KEY`: obbligatoria in produzione (`DJANGO_DEBUG=0`)
    - `ALLOWED_HOSTS`: lista separata da virgole (solo se `DJANGO_DEBUG=0` o per deploy)
    - `CSRF_TRUSTED_ORIGINS`: lista separata da virgole (solo se serve)

    ### Feature flags

    - `PREDLAB_V2`: `1` abilita funzionalità sperimentali (vedi `mmportal/settings.py`)

    ## File di riferimento

    - Template env: `.env.example`
    - Settings: `mmportal/settings.py`
    - URL routing: `mmportal/urls.py`

    !!! warning "Produzione"
        Con `DJANGO_DEBUG=0`, l’app valida `DJANGO_SECRET_KEY` e rifiuta valori deboli o “di default”.

=== "EN"
    ## Main environment variables

    ## Quick find

    | I want… | Go to… |
    | --- | --- |
    | run the stack | `Guides → Operations` |
    | troubleshooting | `Guides → Troubleshooting` |
    | main endpoints | `Reference → Endpoints` |

    ## Terminology

    - **env var**: environment variable read by `mmportal/settings.py`.
    - **feature flag**: toggle that enables/disables a feature.

    ### Django

    - `DJANGO_DEBUG`: `1` (default) / `0`
    - `DJANGO_SECRET_KEY`: required in production (`DJANGO_DEBUG=0`)
    - `ALLOWED_HOSTS`: comma-separated list (needed when `DJANGO_DEBUG=0` or in deploys)
    - `CSRF_TRUSTED_ORIGINS`: comma-separated list (only if needed)

    ### Feature flags

    - `PREDLAB_V2`: `1` enables experimental features (see `mmportal/settings.py`)

    ## Reference files

    - env template: `.env.example`
    - settings: `mmportal/settings.py`
    - URL routing: `mmportal/urls.py`

    !!! warning "Production"
        With `DJANGO_DEBUG=0`, the app validates `DJANGO_SECRET_KEY` and rejects weak or “default-like” values.
