# Operations (run & deploy)

Questa pagina documenta i modi principali per avviare bmyCure4MM e gli elementi operativi (servizi, static/media, logging).

## Avvio in sviluppo

1. Crea un virtualenv e installa dipendenze:

```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
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

- Static: `STATIC_URL=/static/`, `STATIC_ROOT=staticfiles/`, `STATICFILES_DIRS=[mmportal/static]`
- Media: `MEDIA_URL=/media/`, `MEDIA_ROOT=media/`

In produzione:

```bash
python3 manage.py collectstatic
```

## Celery

Script presenti:
- `start_celery.sh`
- `manage_services.sh`

!!! tip "Dev"
    Per debugging rapido puoi impostare `CELERY_TASK_ALWAYS_EAGER=1` (se previsto dal setup) per far girare task in-process.

## Logging

Directory: `logs/`

File tipici:
- `logs/django.log`
- `logs/activity.log`
- `logs/embed_debug.log`
- `logs/celery_tasks.log`

La configurazione è in `mmportal/settings.py` (`LOGGING`).

