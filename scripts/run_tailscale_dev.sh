#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "${REPO_ROOT}"

REQUESTED_PORT="${PORT:-}"

if [[ ! -x "./venv/bin/python" ]]; then
    echo "ERROR: ./venv/bin/python not found or not executable. Create the virtualenv first." >&2
    exit 1
fi

if [[ -f ".env" ]]; then
    set -a
    # shellcheck disable=SC1091
    source ".env"
    set +a
fi

PORT="${REQUESTED_PORT:-${PORT:-8001}}"
if ! [[ "${PORT}" =~ ^[0-9]+$ ]] || (( PORT < 1 || PORT > 65535 )); then
    echo "ERROR: PORT must be an integer between 1 and 65535; got '${PORT}'." >&2
    exit 1
fi

trim_value() {
    local value="$1"
    value="${value#"${value%%[![:space:]]*}"}"
    value="${value%"${value##*[![:space:]]}"}"
    printf '%s' "${value}"
}

append_csv_value() {
    local current="$1"
    local raw_value="$2"
    local value
    value="$(trim_value "${raw_value}")"
    if [[ -z "${value}" || "${value}" == "*" ]]; then
        printf '%s' "${current}"
        return 0
    fi
    case ",${current}," in
        *,"${value}",*) printf '%s' "${current}" ;;
        ,) printf '%s' "${value}" ;;
        *) printf '%s,%s' "${current}" "${value}" ;;
    esac
}

append_csv_values() {
    local current="$1"
    local raw_list="$2"
    local old_ifs="${IFS}"
    local item
    IFS=','
    for item in ${raw_list}; do
        current="$(append_csv_value "${current}" "${item}")"
    done
    IFS="${old_ifs}"
    printf '%s' "${current}"
}

TAILSCALE_IP=""
if command -v tailscale >/dev/null 2>&1; then
    TAILSCALE_IP="$(tailscale ip -4 2>/dev/null | head -n 1 || true)"
fi

DJANGO_ALLOWED_HOSTS=""
DJANGO_ALLOWED_HOSTS="$(append_csv_value "${DJANGO_ALLOWED_HOSTS}" "localhost")"
DJANGO_ALLOWED_HOSTS="$(append_csv_value "${DJANGO_ALLOWED_HOSTS}" "127.0.0.1")"
DJANGO_ALLOWED_HOSTS="$(append_csv_value "${DJANGO_ALLOWED_HOSTS}" "0.0.0.0")"
if [[ -n "${TAILSCALE_IP}" ]]; then
    DJANGO_ALLOWED_HOSTS="$(append_csv_value "${DJANGO_ALLOWED_HOSTS}" "${TAILSCALE_IP}")"
fi
DJANGO_ALLOWED_HOSTS="$(append_csv_values "${DJANGO_ALLOWED_HOSTS}" "${DJANGO_EXTRA_ALLOWED_HOSTS:-}")"

DJANGO_CSRF_TRUSTED_ORIGINS=""
DJANGO_CSRF_TRUSTED_ORIGINS="$(append_csv_value "${DJANGO_CSRF_TRUSTED_ORIGINS}" "http://localhost:${PORT}")"
DJANGO_CSRF_TRUSTED_ORIGINS="$(append_csv_value "${DJANGO_CSRF_TRUSTED_ORIGINS}" "http://127.0.0.1:${PORT}")"
if [[ -n "${TAILSCALE_IP}" ]]; then
    DJANGO_CSRF_TRUSTED_ORIGINS="$(append_csv_value "${DJANGO_CSRF_TRUSTED_ORIGINS}" "http://${TAILSCALE_IP}:${PORT}")"
fi
DJANGO_CSRF_TRUSTED_ORIGINS="$(append_csv_values "${DJANGO_CSRF_TRUSTED_ORIGINS}" "${DJANGO_EXTRA_CSRF_ORIGINS:-}")"

export DJANGO_ALLOWED_HOSTS
export DJANGO_CSRF_TRUSTED_ORIGINS

echo "Django dev server binding: 0.0.0.0:${PORT}"
echo "Local desktop URL: http://127.0.0.1:${PORT}/"
if [[ -n "${TAILSCALE_IP}" ]]; then
    echo "Tailscale IP: ${TAILSCALE_IP}"
    echo "Phone/Tailscale URL: http://${TAILSCALE_IP}:${PORT}/"
    echo "Research cockpit example: http://${TAILSCALE_IP}:${PORT}/research/patient/4/cockpit/"
    echo "Developer console example: http://${TAILSCALE_IP}:${PORT}/research/developer/"
else
    echo "Tailscale IP: unavailable"
    echo "Phone/Tailscale URL: http://<TAILSCALE_IP>:${PORT}/ or http://<TAILSCALE_MAGICDNS_NAME>:${PORT}/"
    echo "Research cockpit example: http://<TAILSCALE_IP>:${PORT}/research/patient/<patient_id>/cockpit/"
    echo "Developer console example: http://<TAILSCALE_IP>:${PORT}/research/developer/"
    echo "Set DJANGO_EXTRA_ALLOWED_HOSTS and DJANGO_EXTRA_CSRF_ORIGINS for MagicDNS or manually discovered Tailscale hosts."
fi
echo "DJANGO_ALLOWED_HOSTS=${DJANGO_ALLOWED_HOSTS}"
echo "DJANGO_CSRF_TRUSTED_ORIGINS=${DJANGO_CSRF_TRUSTED_ORIGINS}"
echo "Note: do not open http://0.0.0.0:${PORT}/ in the browser; 0.0.0.0 is only the bind address."

./venv/bin/python manage.py check
./venv/bin/python manage.py runserver "0.0.0.0:${PORT}"