# Configuration (env, settings, runtime)

=== "IT"
    ## Variabili d’ambiente principali

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
