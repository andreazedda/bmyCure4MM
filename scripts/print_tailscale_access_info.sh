#!/usr/bin/env bash
set -euo pipefail

port_status() {
    local port="$1"
    if lsof -nP -iTCP:"${port}" -sTCP:LISTEN >/dev/null 2>&1; then
        echo "in use"
    else
        echo "free"
    fi
}

TAILSCALE_IP=""
if command -v tailscale >/dev/null 2>&1; then
    TAILSCALE_IP="$(tailscale ip -4 2>/dev/null | head -n 1 || true)"
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
if [[ -n "${TAILSCALE_IP}" ]]; then
    printf 'Suggested phone URL: http://%s:8001/\n' "${TAILSCALE_IP}"
else
    printf 'Suggested phone URL: http://<TAILSCALE_IP>:8001/\n'
fi
printf 'Note: 0.0.0.0 is only the server bind address, not the phone URL.\n'
