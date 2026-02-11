# Deploy su bmycure4mm.duckdns.org (Docker + Nginx + Let's Encrypt)

Questa repo è pubblica: **non committare mai** `.env`, certificati, chiavi o database.

## Prerequisiti (server)

- DNS DuckDNS già puntato all'IP del server: `bmycure4mm.duckdns.org`
- Porte aperte: **80** e **443** (e SSH)
- Docker + Docker Compose plugin installati

## 1) Clone e configurazione `.env`

Sul server:

```bash
git clone https://github.com/andreazedda/bmyCure4MM.git
cd bmyCure4MM
cp .env.example .env
```

Modifica `.env` (valori tipici per prod):

- `DJANGO_DEBUG=0`
- `DJANGO_SECRET_KEY=<genera una chiave forte>`
- `DOMAIN=bmycure4mm.duckdns.org`
- `ALLOWED_HOSTS=bmycure4mm.duckdns.org`
- `CSRF_TRUSTED_ORIGINS=https://bmycure4mm.duckdns.org`
- `DJANGO_SECURE_SSL_REDIRECT=1` (solo dopo aver attivato HTTPS)
- `DJANGO_SESSION_COOKIE_SECURE=1`
- `DJANGO_CSRF_COOKIE_SECURE=1`

Generazione SECRET_KEY (esegui localmente o sul server):

```bash
python -c 'from django.core.management.utils import get_random_secret_key; print(get_random_secret_key())'
```

## 2) Avvio stack HTTP (per ottenere il certificato)

Da `bmyCure4MM/deploy`:

```bash
cd deploy
docker compose --env-file ../.env -f docker-compose.prod.yml up -d --build
```

Verifica HTTP:

- `http://bmycure4mm.duckdns.org` deve rispondere (senza TLS per ora)

## 3) Ottenere certificato Let's Encrypt

Sempre in `deploy/`:

```bash
docker compose --env-file ../.env -f docker-compose.prod.yml run --rm certbot \
  certonly --webroot -w /var/www/certbot \
  -d "$DOMAIN" \
  --email "$CERTBOT_EMAIL" --agree-tos --no-eff-email
```

## 4) Abilitare HTTPS

Avvia Nginx con config HTTPS tramite override:

```bash
docker compose --env-file ../.env -f docker-compose.prod.yml -f docker-compose.ssl.yml up -d nginx
```

Poi abilita redirect HTTPS in Django (in `.env`):

- `DJANGO_SECURE_SSL_REDIRECT=1`
- `DJANGO_SESSION_COOKIE_SECURE=1`
- `DJANGO_CSRF_COOKIE_SECURE=1`

Applica:

```bash
docker compose --env-file ../.env -f docker-compose.prod.yml up -d web celery
```

## 5) Rinnovo certificati

Consigliato: cron sul server (es. 2 volte al giorno) per provare rinnovo e reload di Nginx.
Esempio (crontab):

```bash
0 3,15 * * * cd /path/to/bmyCure4MM/deploy && \
  docker compose --env-file ../.env -f docker-compose.prod.yml run --rm certbot renew && \
  docker compose --env-file ../.env -f docker-compose.prod.yml -f docker-compose.ssl.yml exec -T nginx nginx -s reload
```

## Note su dati e persistenza

- SQLite persistente: il compose usa `DJANGO_SQLITE_PATH=/data/db.sqlite3` su volume dedicato.
- Media e static:
  - `media_data` per upload/report
  - `static_data` per `collectstatic`

## Troubleshooting

- Se Nginx HTTPS non parte, verifica che esistano i file in `/etc/letsencrypt/live/$DOMAIN/` dentro il volume.
- Se Django rifiuta di avviarsi in prod: controlla `DJANGO_SECRET_KEY` (obbligatoria con `DJANGO_DEBUG=0`).
