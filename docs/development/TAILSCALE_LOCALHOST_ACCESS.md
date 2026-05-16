# Tailscale Localhost Access

## 1. Purpose

This guide adds a tailnet-local developer access mode for the Django development server so a local instance can be reached from another Tailscale-connected device.

This is for local research and development only. It does not change the clinical, twin, or what-if logic.

## 2. Security boundary

- This mode is for tailnet-local access only.
- Do not expose the Django development server with public port forwarding.
- Do not use `runserver 0.0.0.0:8000` as a production deployment method.
- Do not share access URLs outside the tailnet.
- Keep PHI out of Git, public artifacts, and copied logs.

## 3. Required preconditions

- The host machine is signed in to Tailscale.
- The client device is connected to the same tailnet.
- The repository virtual environment exists at `./venv/bin/python`.
- Django host and CSRF allowlists include the Tailscale IP or MagicDNS hostname.

Verify Tailscale and Django basics:

```bash
tailscale status
tailscale ip -4

./venv/bin/python manage.py check
./venv/bin/python manage.py runserver 0.0.0.0:8000
```

## 4. How to get the Tailscale IP

Run:

```bash
tailscale status
tailscale ip -4
```

Use the IPv4 address from `tailscale ip -4` as `<TAILSCALE_IP>`.

If MagicDNS is enabled for your tailnet, you can also use the device hostname such as `<TAILSCALE_MAGICDNS_HOST>`.

## 5. Required .env variables

The repository supports the preferred environment variables below and still accepts the legacy aliases `ALLOWED_HOSTS` and `CSRF_TRUSTED_ORIGINS`.

Example local shell exports:

```bash
export DJANGO_ALLOWED_HOSTS="localhost,127.0.0.1,0.0.0.0,<TAILSCALE_IP>,<TAILSCALE_MAGICDNS_HOST>"
export DJANGO_CSRF_TRUSTED_ORIGINS="http://localhost:8000,http://127.0.0.1:8000,http://<TAILSCALE_IP>:8000,http://<TAILSCALE_MAGICDNS_HOST>:8000"
```

Example local `.env` entries:

```dotenv
DJANGO_ALLOWED_HOSTS=localhost,127.0.0.1,0.0.0.0,100.x.y.z,my-machine.tailnet-name.ts.net
DJANGO_CSRF_TRUSTED_ORIGINS=http://localhost:8000,http://127.0.0.1:8000,http://100.x.y.z:8000,http://my-machine.tailnet-name.ts.net:8000
```

Notes:

- `DJANGO_ALLOWED_HOSTS` and `DJANGO_CSRF_TRUSTED_ORIGINS` are comma-separated.
- Whitespace is stripped and empty entries are ignored.
- `localhost` and `127.0.0.1` stay allowed by default.
- The repository does not auto-load `.env` during direct `manage.py` commands. Either export the variables in your shell or use `./scripts/run_tailscale_dev.sh`, which sources `.env` if present.

## 6. Run command

Direct command:

```bash
./venv/bin/python manage.py runserver 0.0.0.0:8000
```

Helper script:

```bash
./scripts/run_tailscale_dev.sh
```

The helper script:

- checks for `./venv/bin/python`
- sources `.env` if present
- prints the detected Tailscale IPv4 address when the `tailscale` CLI is installed
- runs `./venv/bin/python manage.py check`
- starts Django on `0.0.0.0:8000`

## 7. Access URL from another device

Do not use `http://0.0.0.0:8000/` in the browser.

Use one of these instead:

```text
http://<TAILSCALE_IP>:8000/
http://<TAILSCALE_MAGICDNS_HOST>:8000/
```

Example access pattern:

```text
http://<TAILSCALE_IP>:8000/
```

## 8. Troubleshooting

If the browser cannot connect:

- verify both devices are connected to the same tailnet
- verify `tailscale status`
- verify the host command is binding to `0.0.0.0:8000`
- verify `DJANGO_ALLOWED_HOSTS` contains the Tailscale IP or MagicDNS name
- verify `DJANGO_CSRF_TRUSTED_ORIGINS` contains the full origin including scheme and port
- verify the OS firewall allows inbound TCP/8000 on the Tailscale interface

Firewall notes:

- Linux with ufw:

```bash
sudo ufw allow in on tailscale0 to any port 8000 proto tcp
```

- macOS: allow incoming connections for Python or Django if prompted by the firewall
- Windows: allow Python through Windows Defender Firewall for private networks
