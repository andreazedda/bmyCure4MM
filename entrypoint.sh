#!/bin/sh
set -eu

# Optional one-time setup tasks. Controlled via env flags so dev stays fast.

if [ "${DJANGO_MIGRATE:-0}" = "1" ]; then
  echo "[entrypoint] Running migrations..."
  python manage.py migrate --noinput
fi

if [ "${DJANGO_COLLECTSTATIC:-0}" = "1" ]; then
  echo "[entrypoint] Collecting static files..."
  python manage.py collectstatic --noinput
fi

exec "$@"
