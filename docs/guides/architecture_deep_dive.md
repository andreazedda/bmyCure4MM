# Architecture Deep Dive

=== "IT"
    Questa pagina raccoglie note “engineering” (pattern, sicurezza, performance, deploy) a complemento di `Guides → Architecture`.

    ## Se cerchi…

    | Voglio… | Vai a… |
    | --- | --- |
    | overview architettura | `Guides → Architecture` |
    | run/deploy operativo | `Guides → Operations` |
    | configurazione env | `Reference → Configuration` |
    | oggetti DB | `Guides → Database` + `Reference → Database Objects` |

    ## Terminologia

    - **N+1**: antipattern di query ripetute (riducibile con `select_related/prefetch_related`).
    - **Broker**: coda messaggi (es. Redis) usata da Celery.
    - **Worker**: processo che esegue job asincroni (Celery).

    ## Struttura repo (alta-level)

    ```text
    mmportal/        # Django project (settings/urls/wsgi)
    clinic/          # patient management
    simulator/       # scenarios + PK/PD simulation + optimization
    chemtools/       # chem-informatics jobs/tools
    docs_viewer/     # in-app docs viewer
    templates/       # base templates
    mmportal/static/ # static assets
    media/           # runtime artifacts (plots/csv/json)
    logs/            # log files
    ```

    ## Design patterns

    | Area | Pattern | Perché |
    | --- | --- | --- |
    | Models | “fat models” | logica vicino ai dati, KPI/utility riusabili |
    | Views | HTMX partials | UI reattiva senza SPA |
    | Tasks | Celery jobs | isolare computazioni lente (chemtools, batch) |

    ## Sicurezza (pratica)

    - Auth Django standard
    - Validazione input (forms/serializers)
    - Secrets via env (`DJANGO_SECRET_KEY`, ecc.)
    - In produzione: HTTPS, ALLOWED_HOSTS, CSRF_TRUSTED_ORIGINS

    ## Performance

    | Area | Tecnica | Note |
    | --- | --- | --- |
    | DB | indici + select_related/prefetch_related | riduce N+1 |
    | Async | Celery + broker | sposta lavoro fuori request |
    | Static | `collectstatic` | produzione più veloce |

    ## Logging & monitoring

    File principali (vedi `mmportal/settings.py`):

    - `logs/django.log`
    - `logs/activity.log`
    - `logs/embed_debug.log`
    - `logs/celery_tasks.log`

    ## Deploy checklist

    - `DJANGO_DEBUG=0`
    - `DJANGO_SECRET_KEY` forte e non di default
    - configurare DB (Postgres consigliato)
    - configurare static/media (WhiteNoise o reverse proxy)
    - backup strategy + error tracking (es. Sentry)

=== "EN"
    This page collects engineering notes (patterns, security, performance, deployment) as a companion to `Guides → Architecture`.

    ## Quick find

    | I want… | Go to… |
    | --- | --- |
    | architecture overview | `Guides → Architecture` |
    | operational run/deploy | `Guides → Operations` |
    | env configuration | `Reference → Configuration` |
    | DB objects | `Guides → Database` + `Reference → Database Objects` |

    ## Terminology

    - **N+1**: repeated query antipattern (reduced via `select_related/prefetch_related`).
    - **Broker**: message queue backend (e.g., Redis) for Celery.
    - **Worker**: process executing background jobs (Celery).

    ## Repo structure (high level)

    ```text
    mmportal/        # Django project (settings/urls/wsgi)
    clinic/          # patient management
    simulator/       # scenarios + PK/PD simulation + optimization
    chemtools/       # chem-informatics jobs/tools
    docs_viewer/     # in-app docs viewer
    templates/       # base templates
    mmportal/static/ # static assets
    media/           # runtime artifacts (plots/csv/json)
    logs/            # log files
    ```

    ## Design patterns

    | Area | Pattern | Why |
    | --- | --- | --- |
    | Models | fat models | keep logic close to data, reuse KPI utilities |
    | Views | HTMX partials | reactive UX without a full SPA |
    | Tasks | Celery jobs | isolate slow computations from requests |

    ## Security

    - Django auth
    - input validation (forms/serializers)
    - secrets via env (`DJANGO_SECRET_KEY`, etc.)
    - in production: HTTPS, ALLOWED_HOSTS, CSRF_TRUSTED_ORIGINS

    ## Performance

    | Area | Technique | Notes |
    | --- | --- | --- |
    | DB | indexes + select_related/prefetch_related | reduce N+1 |
    | Async | Celery + broker | move heavy work off request path |
    | Static | `collectstatic` | faster production serving |

    ## Logging

    Key files (see `mmportal/settings.py`):

    - `logs/django.log`
    - `logs/activity.log`
    - `logs/embed_debug.log`
    - `logs/celery_tasks.log`

    ## Deployment checklist

    - `DJANGO_DEBUG=0`
    - strong `DJANGO_SECRET_KEY`
    - production DB (usually Postgres)
    - static/media strategy (WhiteNoise or reverse proxy)
    - backups + error tracking (e.g., Sentry)
