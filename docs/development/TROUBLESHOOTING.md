# Troubleshooting Guide - MM Portal

## Celery Worker Crashes (SIGSEGV)

### Problema
Job rimangono in stato "Queued" e nei log di Celery appare:
```
ERROR/MainProcess] Process 'ForkPoolWorker-X' pid:XXXXX exited with 'signal 11 (SIGSEGV)'
WorkerLostError: Worker exited prematurely: signal 11 (SIGSEGV)
```

### Causa
RDKit e altre librerie scientifiche non sono fork-safe su macOS. Quando Celery usa il pool di processi prefork (default), i worker crashano con segmentation fault.

### Soluzione
Usare `--pool=solo` invece del default prefork:

```bash
celery -A mmportal worker --loglevel=info --pool=solo
```

**Nota:** `--pool=solo` esegue i task in modo sincrono in un singolo thread, quindi non ci sarà parallelismo. Per applicazioni in produzione con alto carico, considera:
1. Deploy su Linux (dove fork funziona meglio)
2. Usare `--pool=threads` (meno performante ma più stabile)
3. Containerizzare con Docker

### Configurazione Automatica
Tutti gli script di avvio sono già configurati per usare `--pool=solo`:
- `python3 manage.py runserver` ✅
- `./start_dev.sh` ✅
- `./manage_services.sh start` ✅
- `docker-compose up` ✅

## Job Bloccati in "Queued"

### Diagnosi

**1. Verifica che Redis sia in esecuzione:**
```bash
redis-cli ping
# Output atteso: PONG
```

**2. Verifica che Celery worker sia attivo:**
```bash
ps aux | grep "celery.*worker"
```

**3. Controlla i log di Celery:**
```bash
tail -f logs/celery.log
```

### Soluzioni

**Se Redis non è in esecuzione:**
```bash
redis-server
```

**Se Celery non è in esecuzione:**
```bash
celery -A mmportal worker --loglevel=info --pool=solo
```

**Se i log mostrano errori:**
- Controlla che tutte le dipendenze siano installate: `pip install -r requirements.txt`
- Verifica che il virtual environment sia attivato
- Riavvia i servizi

## Barra di Progresso Non Si Aggiorna

### Causa Possibile 1: Job Non Processato
Se il job non viene processato da Celery, la barra rimane ferma.

**Soluzione:** Vedi sezione "Job Bloccati in Queued"

### Causa Possibile 2: Polling JavaScript Disabilitato
Il browser potrebbe bloccare il polling HTMX.

**Soluzione:** 
- Ricarica la pagina
- Controlla la console del browser per errori JavaScript
- Verifica che l'endpoint `/chem/job/<id>/status.json` risponda correttamente

## Modalità Sincrona (Senza Redis/Celery)

Se non vuoi usare Redis/Celery durante lo sviluppo:

**Opzione 1: Variabile d'ambiente**
```bash
export CELERY_TASK_ALWAYS_EAGER=True
python3 manage.py runserver 8001
```

**Opzione 2: Settings.py**
Decommenta in `mmportal/settings.py`:
```python
CELERY_TASK_ALWAYS_EAGER = True
CELERY_TASK_EAGER_PROPAGATES = True
```

**Nota:** In modalità sincrona i job bloccano la richiesta HTTP fino al completamento (potrebbero volerci minuti).

## Errori Comuni

### ImportError: No module named 'celery'
```bash
pip install celery redis
```

### ConnectionError: Error 61 connecting to localhost:6379
Redis non è in esecuzione:
```bash
redis-server
```

### AttributeError: 'ChemJob' object has no attribute 'update_progress'
Migrations non applicate:
```bash
python3 manage.py migrate chemtools
```

### SIGSEGV durante test
Usa pytest con `--forked` disabilitato o Django test runner:
```bash
python3 manage.py test chemtools.tests
```

## Performance Tips

### Sviluppo Locale
- Usa `--pool=solo` (già configurato)
- Un solo worker è sufficiente
- Redis in modalità default va bene

### Produzione
- Usa Linux (Ubuntu/Debian preferito)
- Considera `--pool=prefork` con 2-4 worker
- Redis persistente con configurazione ottimizzata
- Supervisord o systemd per gestione processi
- Gunicorn/uWSGI invece di runserver

## Log Files

- **Celery:** `logs/celery.log`
- **Redis:** `/var/log/redis/` o output terminale
- **Django:** Console standard output
- **Supervisord:** `logs/supervisord.log`

## Supporto

Se il problema persiste:
1. Controlla tutti i log
2. Verifica versioni dipendenze: `pip list | grep -E "celery|redis|rdkit"`
3. Testa con modalità sincrona per escludere problemi Celery
4. Apri issue su GitHub con log completi
