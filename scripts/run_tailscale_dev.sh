#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

if [[ ! -x "./venv/bin/python" ]]; then
    echo "Error: ./venv/bin/python not found or not executable. Create the virtualenv first." >&2
    exit 1
fi

if [[ -f ".env" ]]; then
    set -a
    # shellcheck disable=SC1091
    source ".env"
    set +a
fi

TAILSCALE_IP=""
if command -v tailscale >/dev/null 2>&1; then
    TAILSCALE_IP="$(tailscale ip -4 2>/dev/null | head -n 1 || true)"
fi

echo "Django dev server binding: 0.0.0.0:8000"
echo "Local URL: http://127.0.0.1:8000/"
if [[ -n "${TAILSCALE_IP}" ]]; then
    echo "Tailscale IP: ${TAILSCALE_IP}"
    echo "Tailscale URL: http://${TAILSCALE_IP}:8000/"
else
    echo "Tailscale IP: unavailable"
    echo "Tailscale URL: configure DJANGO_ALLOWED_HOSTS and DJANGO_CSRF_TRUSTED_ORIGINS, then use http://<TAILSCALE_IP>:8000/"
fi
echo "Reminder: 0.0.0.0 is not the browser URL."

./venv/bin/python manage.py check
./venv/bin/python manage.py runserver 0.0.0.0:8000