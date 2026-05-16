# Tailscale Localhost Access

## Concept

`0.0.0.0` is a server bind address. It tells Django to listen on all local network interfaces.

It is not a browser URL. From a phone or another device on the same Tailscale tailnet, use the host machine's Tailscale IP or MagicDNS hostname:

```text
http://<TAILSCALE_IP>:8001/
http://<TAILSCALE_MAGICDNS_NAME>:8001/
```

Never open `http://0.0.0.0:8001/` in the browser.

## Preconditions

- Both devices are connected to the same Tailscale tailnet.
- You are logged into a Django account that is allowed to view the target page.
- The local Django development server is running.
- The firewall allows inbound TCP on the selected port through the Tailscale interface.
- `DJANGO_ALLOWED_HOSTS` includes the Tailscale IP or MagicDNS hostname.
- `DJANGO_CSRF_TRUSTED_ORIGINS` includes the full origin with scheme and port.

This mode is for local development and research validation only. It does not change production, Kubernetes, mathematical model, clinical model, counterfactual, calibration, or patient-data behavior.

## Commands

Get the Tailscale IP:

```bash
tailscale status
tailscale ip -4
```

Start the default port `8001`:

```bash
scripts/run_tailscale_dev.sh
```

Start port `8000` when it is free:

```bash
PORT=8000 scripts/run_tailscale_dev.sh
```

Manual fallback:

```bash
TAILSCALE_IP="$(tailscale ip -4 | head -n1)"
export DJANGO_ALLOWED_HOSTS="localhost,127.0.0.1,0.0.0.0,${TAILSCALE_IP}"
export DJANGO_CSRF_TRUSTED_ORIGINS="http://localhost:8001,http://127.0.0.1:8001,http://${TAILSCALE_IP}:8001"
./venv/bin/python manage.py check
./venv/bin/python manage.py runserver 0.0.0.0:8001
```

For MagicDNS, add the hostname as a local extra when using the helper:

```bash
DJANGO_EXTRA_ALLOWED_HOSTS="<TAILSCALE_MAGICDNS_NAME>" \
DJANGO_EXTRA_CSRF_ORIGINS="http://<TAILSCALE_MAGICDNS_NAME>:8001" \
scripts/run_tailscale_dev.sh
```

Use placeholders in committed docs and examples. Do not commit real tailnet names or real private host details.

## URLs

Desktop:

```text
http://127.0.0.1:8001/
```

Phone:

```text
http://<TAILSCALE_IP>:8001/
```

Research cockpit:

```text
http://<TAILSCALE_IP>:8001/research/patient/4/cockpit/
```

Developer console:

```text
http://<TAILSCALE_IP>:8001/research/developer/
```

The same patterns work with a MagicDNS hostname:

```text
http://<TAILSCALE_MAGICDNS_NAME>:8001/
```

## Status Helper

For a read-only local summary:

```bash
scripts/print_tailscale_access_info.sh
```

It prints the host name, detected Tailscale IP if available, whether ports `8000` and `8001` are already listening, and the suggested command/URL.

## Troubleshooting

If the phone cannot connect:

- Verify `tailscale status` on both devices.
- Verify `tailscale ip -4` on the host machine.
- Verify the server is bound to `0.0.0.0`, not only `127.0.0.1`.
- Verify the phone URL uses the Tailscale IP or MagicDNS hostname, not `0.0.0.0`.
- Verify `DJANGO_ALLOWED_HOSTS` includes the Tailscale IP or MagicDNS hostname.
- Verify `DJANGO_CSRF_TRUSTED_ORIGINS` includes the full origin with scheme and port.
- Verify the macOS firewall allows incoming Python/Django connections.
- Verify no other server occupies the selected port.
- Try `PORT=8001 scripts/run_tailscale_dev.sh` if `8000` is busy.

On Linux with `ufw`, the equivalent tailnet-local firewall rule is:

```bash
sudo ufw allow in on tailscale0 to any port 8001 proto tcp
```

## Security Boundary

- Tailnet-local development only.
- No public port forwarding.
- No public tunnel.
- No production use of Django `runserver`.
- No `ALLOWED_HOSTS=*` default.
- No PHI, raw clinical documents, media artifacts, SQLite databases, screenshots, or generated reports in Git.
- Do not share the access URL outside the tailnet.
