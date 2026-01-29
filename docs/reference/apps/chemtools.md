# App `chemtools`

## Scopo

Strumenti di chem-informatics integrati:

- parametri farmaco (utility)
- visualizzazione binding (HTML)
- similarity search (es. ligandi)

## Modello principale

`chemtools/models.py`:

- `ChemJob`: log, progress, output file, preferenze API (JSON)

## Esecuzione asincrona

Molti tool vengono eseguiti come job (spesso via Celery): vedi `chemtools/tasks.py`.

