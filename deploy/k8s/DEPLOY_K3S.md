# Deploy su k3s (cluster) via Cloudflare Tunnel

Questo cluster espone bmyCure4MM tramite Cloudflare Tunnel sul dominio pubblico `bmycure4mm.clusterlab.uk`.
Il tunnel deve inoltrare le richieste HTTP al service `web` del namespace `bmycure4mm` sulla porta `8001`.

## Prerequisiti

- Hai accesso a `sudo kubectl` (in questo nodo sembra ok)
- L'immagine container è pubblicata su GHCR: `ghcr.io/andreazedda/bmycure4mm:latest`
- Un token GHCR (PAT GitHub) con permesso `read:packages` per permettere al cluster di fare pull (se l'immagine non è pubblica)
- Esiste un tunnel Cloudflare configurato per `bmycure4mm.clusterlab.uk`
- Il tunnel inoltra l'header `Host: bmycure4mm.clusterlab.uk`
- Cloudflare inoltra `X-Forwarded-Proto: https` se vuoi attivare redirect SSL e cookie sicuri in Django

## 1) Pubblicare l'immagine su GHCR

Nel repo è presente workflow GitHub Actions: `.github/workflows/ghcr-image.yml`.

Serve fare push su `master` per triggerare la build.

## 2) Creare secret in cluster (senza committare nulla)

Esegui sul nodo (NON stampa la chiave):

```bash
sudo kubectl create namespace bmycure4mm --dry-run=client -o yaml | sudo kubectl apply -f -
DJANGO_SECRET_KEY="$(python3 -c 'from django.core.management.utils import get_random_secret_key; print(get_random_secret_key())')"
sudo kubectl -n bmycure4mm delete secret bmycure4mm-secrets --ignore-not-found
sudo kubectl -n bmycure4mm create secret generic bmycure4mm-secrets --from-literal=DJANGO_SECRET_KEY="$DJANGO_SECRET_KEY"
unset DJANGO_SECRET_KEY
```

Se l'immagine su GHCR non è pubblica, crea anche un `imagePullSecret` (non viene committato):

```bash
GHCR_USER="<github-username>"
GHCR_PAT="<github-pat-with-read:packages>"

sudo kubectl -n bmycure4mm delete secret ghcr-creds --ignore-not-found
sudo kubectl -n bmycure4mm create secret docker-registry ghcr-creds \
	--docker-server=ghcr.io \
	--docker-username="$GHCR_USER" \
	--docker-password="$GHCR_PAT"

unset GHCR_PAT
```

## 3) Applicare i manifest

```bash
sudo kubectl apply -f deploy/k8s/bmycure4mm.yaml
```

Il manifest imposta già:

- `ALLOWED_HOSTS=bmycure4mm.clusterlab.uk`
- `CSRF_TRUSTED_ORIGINS=https://bmycure4mm.clusterlab.uk`
- probe HTTP con header `Host: bmycure4mm.clusterlab.uk`

## 4) Configurare Cloudflare Tunnel

Nel pannello Cloudflare Zero Trust, oppure nella config locale del tunnel, crea una route pubblica:

- Hostname: `bmycure4mm.clusterlab.uk`
- Service: `http://web.bmycure4mm.svc.cluster.local:8001`

Se il tunnel gira fuori dal cluster e non risolve il DNS Kubernetes interno, punta invece a un endpoint raggiungibile dal processo `cloudflared` (ad esempio IP del nodo + porta esposta internamente).

## 5) Verifiche dopo il deploy

Controlli cluster:

```bash
sudo kubectl -n bmycure4mm get pods -o wide
sudo kubectl -n bmycure4mm rollout status deploy/web
sudo kubectl -n bmycure4mm logs deploy/web --tail=200
```

Controlli applicativi:

```bash
curl -I https://bmycure4mm.clusterlab.uk/
curl https://bmycure4mm.clusterlab.uk/healthz/
```

Verifica anche un login o un form POST reale per escludere errori CSRF dopo il cambio host.

## 6) Hardening HTTPS in Django

Quando il tunnel conferma sempre `X-Forwarded-Proto: https`, puoi attivare in `deploy/k8s/bmycure4mm.yaml`:

- `DJANGO_SECURE_SSL_REDIRECT: "1"`
- `DJANGO_SESSION_COOKIE_SECURE: "1"`
- `DJANGO_CSRF_COOKIE_SECURE: "1"`

Dopo l'attivazione, ricontrolla che non ci siano redirect loop.

## Debug

```bash
sudo kubectl -n bmycure4mm get pods -o wide
sudo kubectl -n bmycure4mm logs deploy/web --tail=200
```
