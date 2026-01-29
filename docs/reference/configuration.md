# Configuration (env, settings, runtime)

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

