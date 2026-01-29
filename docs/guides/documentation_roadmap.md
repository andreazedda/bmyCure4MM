# Documentation Roadmap (seria e professionale)

Questa roadmap definisce le prossime iterazioni per portare la documentazione a livello “software product”.

## Milestone 1 — Stabilizzazione e coerenza (1–2 giorni)

- Rendere **bilingue (IT/EN)** tutte le pagine “core”:
  - Architecture, Database, Operations, Development, Troubleshooting, Reference (config/endpoints)
- Standardizzare stile: `Guides → Documentation Standard`
- Aggiungere “Se cerchi…” come tabella all’inizio di ogni pagina core
- Eliminare ogni warning `mkdocs --strict`

## Milestone 2 — Copertura funzionale (2–4 giorni)

- Documentare in modo approfondito:
  - Difficulty scoring (`simulator/difficulty_scoring.py`)
  - Virtual patients (`simulator/virtual_patients.py`)
  - Prognosis (`simulator/prognosis.py`)
  - Regimen suggester (`simulator/regimen_suggester.py`)
- Aggiungere esempi numerici e “worked examples” (input → output) per ogni modulo
- Aggiungere figure SVG dedicate (sensitivity, score curves, prognosis curves)

## Milestone 3 — Reference completa (3–6 giorni)

- Endpoints: tabella completa per app (URL → view → template → payload)
- DB: mappa “query cookbook” (ORM + SQL) per use-case frequenti
- API: payload schema (JSON) per simulator/chemtools/docs_viewer
- “Runbooks” operativi (local/dev/docker, celery/redis, troubleshooting)

## Milestone 4 — Esperienza utente docs (continuo)

- Indice globale e cross-linking aggressivo
- “Tour guidati” per persona (clinico/dev/research) con tempo stimato
- Più mappe Mermaid (componenti, sequence, dataflow)
- “Gallery” estesa con esempi e casi reali (demo patients)

## KPI di qualità (misurabili)

- 0 warning in `mkdocs build --strict`
- 100% pagine core in IT/EN
- 100% nuove feature con:
  - overview + reference + example + diagram
- Search: parole chiave minime per ogni modulo (es. `tumor_reduction`, `healthy_loss`, `cohort`, `twin`)

