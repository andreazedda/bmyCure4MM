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

find_tailscale_cli() {
    local candidate=""
    if [[ -n "${TAILSCALE_CLI:-}" && -x "${TAILSCALE_CLI}" ]]; then
        printf '%s' "${TAILSCALE_CLI}"
        return 0
    fi
    if command -v tailscale >/dev/null 2>&1; then
        command -v tailscale
        return 0
    fi
    for candidate in \
        "/Applications/Tailscale.app/Contents/MacOS/tailscale" \
        "/Applications/Tailscale.app/Contents/Resources/tailscale/bin/tailscale"
    do
        if [[ -x "${candidate}" ]]; then
            printf '%s' "${candidate}"
            return 0
        fi
    done
    return 1
}

print_port_blocker() {
    local port="$1"
    local listener_names=""
    local listener_name=""
    local pid=""
    local command_line=""
    local all_loopback_listeners=1

    listener_names="$(lsof -nP -iTCP:"${port}" -sTCP:LISTEN 2>/dev/null | awk 'NR > 1 {print $9}')"
    while IFS= read -r listener_name; do
        [[ -z "${listener_name}" ]] && continue
        case "${listener_name}" in
            127.0.0.1:${port}|[[]::1[]]:${port}) ;;
            *) all_loopback_listeners=0 ;;
        esac
    done <<< "${listener_names}"

    if (( all_loopback_listeners )); then
        echo "Port ${port} is already occupied by a loopback-only listener. Stop that process before Tailscale access can work." >&2
    else
        echo "Port ${port} is already occupied by another listener. Stop that process or choose a different port." >&2
    fi

    while IFS= read -r pid; do
        [[ -z "${pid}" ]] && continue
        command_line="$(ps -p "${pid}" -o command= 2>/dev/null | sed 's/^[[:space:]]*//')"
        printf 'PID: %s\n' "${pid}" >&2
        printf 'Command: %s\n' "${command_line:-unknown}" >&2
    done < <(lsof -tiTCP:"${port}" -sTCP:LISTEN 2>/dev/null)

    if [[ -n "${listener_names}" ]]; then
        printf 'Listener(s):\n%s\n' "${listener_names}" >&2
    fi
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
        ,,) printf '%s' "${value}" ;;
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

TAILSCALE_IP="${TAILSCALE_IP:-}"
TAILSCALE_CLI_PATH="$(find_tailscale_cli || true)"
if [[ -z "${TAILSCALE_IP}" && -n "${TAILSCALE_CLI_PATH}" ]]; then
    TAILSCALE_IP="$("${TAILSCALE_CLI_PATH}" ip -4 2>/dev/null | head -n 1 || true)"
fi

if lsof -nP -iTCP:"${PORT}" -sTCP:LISTEN >/dev/null 2>&1; then
    print_port_blocker "${PORT}"
    exit 1
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

echo "Bind address: 0.0.0.0:${PORT}"
echo "Local URL: http://127.0.0.1:${PORT}/"
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
if [[ -n "${TAILSCALE_IP}" ]]; then
    echo "0.0.0.0 is bind address; phone URL is http://${TAILSCALE_IP}:${PORT}/"
else
    echo "0.0.0.0 is bind address; phone URL is http://<TAILSCALE_IP>:${PORT}/"
fi
if [[ "${NORELOAD:-0}" == "1" ]]; then
    echo "Autoreload: disabled via NORELOAD=1"
fi
echo "Note: do not open http://0.0.0.0:${PORT}/ in the browser; 0.0.0.0 is only the bind address."

./venv/bin/python manage.py check
RUNSERVER_ARGS=(manage.py runserver "0.0.0.0:${PORT}")
if [[ "${NORELOAD:-0}" == "1" ]]; then
    RUNSERVER_ARGS+=(--noreload)
fi

./venv/bin/python "${RUNSERVER_ARGS[@]}"