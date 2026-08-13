# Operations (run & deploy)

=== "IT"
    Questa pagina documenta i modi principali per avviare bmyCure4MM e gli elementi operativi (servizi, static/media, logging).

    ## Se cerchi…

    | Voglio… | Vai a… |
    | --- | --- |
    | variabili d’ambiente | `Reference → Configuration` |
    | endpoint e routing | `Reference → Endpoints` |
    | troubleshooting | `Guides → Troubleshooting` |

    ## Terminologia

    - **Run mode**: modalità di esecuzione (dev, docker, prod).
    - **Static**: asset serviti (CSS/JS) generati con `collectstatic`.
    - **Media**: output runtime (artifact) salvati in `media/`.
    - **Worker**: processo asincrono (Celery) per job lunghi.

    ## Run modes (tabella)

    | Modalità | Quando usarla | Comando |
    | --- | --- | --- |
    | Dev locale | sviluppo/debug | `python3 manage.py runserver` |
    | Docker Compose | stack completo (redis/celery) | `docker compose up --build` |
    | Produzione | deploy | `collectstatic` + WSGI + worker |

    ## Avvio in sviluppo

    1. Crea un virtualenv e installa dipendenze:

    ```bash
    python3.11 -m pip install uv==0.12.3
    uv sync --frozen --extra chemistry
    ```

    2. Variabili d’ambiente minime (dev):

    ```bash
    export DJANGO_DEBUG=1
    export DJANGO_SECRET_KEY='dev-key-not-for-production'
    ```

    3. Migrazioni + runserver:

    ```bash
    python3 manage.py migrate
    python3 manage.py runserver
    ```

    ## Avvio con Docker Compose

    Vedi `docker-compose.yml`. Tipicamente include:

    - Web (Django)
    - Redis (broker Celery)
    - Worker Celery

    ```bash
    docker compose up --build
    ```

    ## Static & media

    | Tipo | Path | Note |
    | --- | --- | --- |
    | static | `STATIC_ROOT=staticfiles/` | `collectstatic` in produzione |
    | media | `MEDIA_ROOT=media/` | artifact e upload |

    In produzione:

    ```bash
    python3 manage.py collectstatic
    ```

    ## Celery

    Script presenti:
    - `start_celery.sh`
    - `manage_services.sh`

    !!! tip "Dev"
        Per debugging rapido puoi impostare `CELERY_TASK_ALWAYS_EAGER=1` (se previsto) per eseguire task inline.

    ## Logging

    Directory: `logs/`

    | File | Contenuto |
    | --- | --- |
    | `logs/django.log` | log Django |
    | `logs/activity.log` | activity logging middleware |
    | `logs/embed_debug.log` | debug embed |
    | `logs/celery_tasks.log` | task celery |

    Config in `mmportal/settings.py` (`LOGGING`).

=== "EN"
    This page describes how to run bmyCure4MM and operational concerns (services, static/media, logging).

    ## Quick find

    | I want… | Go to… |
    | --- | --- |
    | environment variables | `Reference → Configuration` |
    | routing/endpoints | `Reference → Endpoints` |
    | troubleshooting | `Guides → Troubleshooting` |

    ## Terminology

    - **Run mode**: execution mode (dev, docker, prod).
    - **Static**: served assets (CSS/JS) produced by `collectstatic`.
    - **Media**: runtime outputs (artifacts) stored in `media/`.
    - **Worker**: background process (Celery) for long jobs.

    ## Run modes

    | Mode | When | Command |
    | --- | --- | --- |
    | Local dev | development/debug | `python3 manage.py runserver` |
    | Docker Compose | full stack (redis/celery) | `docker compose up --build` |
    | Production | deployment | `collectstatic` + WSGI + workers |

    ## Local development

    ```bash
    python3.11 -m pip install uv==0.12.3
    uv sync --frozen --extra chemistry

    export DJANGO_DEBUG=1
    export DJANGO_SECRET_KEY='dev-key-not-for-production'

    python3 manage.py migrate
    python3 manage.py runserver
    ```

    ## Docker Compose

    ```bash
    docker compose up --build
    ```

    ## Static & media

    | Type | Path | Notes |
    | --- | --- | --- |
    | static | `STATIC_ROOT=staticfiles/` | run `collectstatic` in production |
    | media | `MEDIA_ROOT=media/` | artifacts and uploads |

    ## Logging

    Logs live in `logs/` and are configured in `mmportal/settings.py` (`LOGGING`).
