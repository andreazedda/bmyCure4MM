# Archivio: vecchio deploy DuckDNS / Nginx / Let's Encrypt

Questa guida descrive il vecchio percorso di deploy basato su `bmycure4mm.duckdns.org`, Nginx e certificati Let's Encrypt sull'origine.
Non e' piu' il percorso di produzione attivo.

## Stato attuale

La produzione attuale usa:

- hostname pubblico: `bmycure4mm.clusterlab.uk`
- esposizione: Cloudflare Tunnel
- deploy applicativo: k3s

Per il runbook aggiornato usa `deploy/k8s/DEPLOY_K3S.md`.

## Perche' questa guida e' stata ritirata

Il vecchio setup presumeva:

- DNS DuckDNS puntato direttamente al server
- Nginx esposto su `80/443`
- emissione e rinnovo certificati con certbot sul server origin

Con Cloudflare Tunnel questi passaggi non sono piu' parte del flusso standard di produzione, quindi mantenerli come istruzioni attive sarebbe fuorviante.

## Note storiche

Se in futuro servirà ripristinare un deploy Docker standalone, andrà aggiornata almeno questa configurazione:

- `DOMAIN=bmycure4mm.clusterlab.uk`
- `ALLOWED_HOSTS=bmycure4mm.clusterlab.uk`
- `CSRF_TRUSTED_ORIGINS=https://bmycure4mm.clusterlab.uk`

La gestione TLS/origine andrà ridefinita in funzione del tunnel o di un eventuale reverse proxy nuovo.

## Sicurezza

La repo resta pubblica: non committare mai `.env`, chiavi, certificati, token Cloudflare o database.
