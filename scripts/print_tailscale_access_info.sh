#!/usr/bin/env bash
set -euo pipefail

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

port_status() {
    local port="$1"
    if lsof -nP -iTCP:"${port}" -sTCP:LISTEN >/dev/null 2>&1; then
        echo "in use"
    else
        echo "free"
    fi
}

TAILSCALE_IP="${TAILSCALE_IP:-}"
TAILSCALE_CLI_PATH="$(find_tailscale_cli || true)"
if [[ -z "${TAILSCALE_IP}" && -n "${TAILSCALE_CLI_PATH}" ]]; then
    TAILSCALE_IP="$("${TAILSCALE_CLI_PATH}" ip -4 2>/dev/null | head -n 1 || true)"
fi

HOSTNAME_VALUE="$(hostname 2>/dev/null || printf 'unknown')"
PORT_8000_STATUS="$(port_status 8000)"
PORT_8001_STATUS="$(port_status 8001)"

printf 'Hostname: %s\n' "${HOSTNAME_VALUE}"
if [[ -n "${TAILSCALE_IP}" ]]; then
    printf 'Tailscale IP: %s\n' "${TAILSCALE_IP}"
else
    printf 'Tailscale IP: unavailable; run tailscale status and tailscale ip -4 on the host.\n'
fi
printf 'Port 8000: %s\n' "${PORT_8000_STATUS}"
printf 'Port 8001: %s\n' "${PORT_8001_STATUS}"
printf 'Suggested run command: scripts/run_tailscale_dev.sh\n'
printf 'Port 8000 override: PORT=8000 scripts/run_tailscale_dev.sh\n'
printf 'Detailed diagnostic: PORT=8001 scripts/check_tailscale_django_access.sh\n'
if [[ -n "${TAILSCALE_IP}" ]]; then
    printf 'Suggested phone URL: http://%s:8001/\n' "${TAILSCALE_IP}"
else
    printf 'Suggested phone URL: http://<TAILSCALE_IP>:8001/\n'
fi
printf 'Note: 0.0.0.0 is only the server bind address, not the phone URL.\n'
