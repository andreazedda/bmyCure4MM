# Deploy su k3s (cluster) via Caddy (duckdns-demo)

Questo cluster non ha IngressClass/cert-manager; l'esposizione HTTP attuale passa da `duckdns-demo/caddy` (NodePort 31080) con routing per host.

## Prerequisiti

- Hai accesso a `sudo kubectl` (in questo nodo sembra ok)
- L'immagine container è pubblicata su GHCR: `ghcr.io/andreazedda/bmycure4mm:latest`

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

## 3) Applicare i manifest

```bash
sudo kubectl apply -f deploy/k8s/bmycure4mm.yaml
```

## 4) Switch routing su bmycure4mm.duckdns.org

Aggiorna `duckdns-demo/caddy-config` per far puntare `bmycure4mm.duckdns.org` al service `web` del namespace `bmycure4mm`.

## Debug

```bash
sudo kubectl -n bmycure4mm get pods -o wide
sudo kubectl -n bmycure4mm logs deploy/web --tail=200
sudo kubectl -n duckdns-demo logs deploy/caddy --tail=200
```
