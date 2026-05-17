#!/usr/bin/env bash

set -euo pipefail

PORT="${PORT:-8001}"
if ! [[ "${PORT}" =~ ^[0-9]+$ ]] || (( PORT < 1 || PORT > 65535 )); then
    echo "ERROR: PORT must be an integer between 1 and 65535; got '${PORT}'." >&2
    exit 1
fi

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

listener_output() {
    lsof -nP -iTCP:"${PORT}" -sTCP:LISTEN 2>/dev/null || true
}

listener_names() {
    printf '%s\n' "$1" | awk 'NR > 1 {print $9}'
}

listener_state() {
    local names="$1"
    local name=""
    local all_loopback=1

    if [[ -z "${names}" ]]; then
        printf 'not listening'
        return 0
    fi

    while IFS= read -r name; do
        [[ -z "${name}" ]] && continue
        case "${name}" in
            127.0.0.1:${PORT}|[[]::1[]]:${PORT}) ;;
            *) all_loopback=0 ;;
        esac
    done <<< "${names}"

    if (( all_loopback )); then
        printf 'loopback-only'
    else
        printf 'network-reachable'
    fi
}

curl_status() {
    local url="$1"
    local header=""
    if header="$(curl -I -m 5 -sS "${url}" | head -n 1)"; then
        printf '%s' "${header}"
    else
        printf 'curl failed'
    fi
}

TAILSCALE_CLI_PATH="$(find_tailscale_cli || true)"
TAILSCALE_IP="${TAILSCALE_IP:-}"
if [[ -z "${TAILSCALE_IP}" && -n "${TAILSCALE_CLI_PATH}" ]]; then
    TAILSCALE_IP="$("${TAILSCALE_CLI_PATH}" ip -4 2>/dev/null | head -n 1 || true)"
fi

LISTENER_OUTPUT="$(listener_output)"
LISTENER_NAMES="$(listener_names "${LISTENER_OUTPUT}")"

PHONE_URL="http://<TAILSCALE_IP>:${PORT}/"
TAILSCALE_CURL_STATUS="not available"
if [[ -n "${TAILSCALE_IP}" ]]; then
    PHONE_URL="http://${TAILSCALE_IP}:${PORT}/"
    TAILSCALE_CURL_STATUS="$(curl_status "${PHONE_URL}")"
fi

printf 'Port: %s\n' "${PORT}"
if [[ -n "${TAILSCALE_IP}" ]]; then
    printf 'Tailscale IP: %s\n' "${TAILSCALE_IP}"
else
    printf 'Tailscale IP: unavailable\n'
fi
printf 'Current listener:\n'
if [[ -n "${LISTENER_NAMES}" ]]; then
    printf '%s\n' "${LISTENER_OUTPUT}"
else
    printf 'none\n'
fi
printf 'Listener classification: %s\n' "$(listener_state "${LISTENER_NAMES}")"
printf 'Local curl status: %s\n' "$(curl_status "http://127.0.0.1:${PORT}/")"
printf 'Tailscale-IP curl status: %s\n' "${TAILSCALE_CURL_STATUS}"
printf 'ALLOWED_HOSTS suggestion: DJANGO_ALLOWED_HOSTS="localhost,127.0.0.1,0.0.0.0,%s"\n' "${TAILSCALE_IP:-<TAILSCALE_IP>}"
printf 'CSRF suggestion: DJANGO_CSRF_TRUSTED_ORIGINS="http://localhost:%s,http://127.0.0.1:%s,http://%s:%s"\n' "${PORT}" "${PORT}" "${TAILSCALE_IP:-<TAILSCALE_IP>}" "${PORT}"
printf 'Phone URL: %s\n' "${PHONE_URL}"
printf 'Cockpit URL: %sresearch/patient/4/cockpit/\n' "${PHONE_URL}"