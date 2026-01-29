# Documentation Roadmap (seria e professionale)

=== "IT"
    Questa roadmap definisce le prossime iterazioni per portare la documentazione a livello “software product”.

    ## Milestone 1 — Stabilizzazione e coerenza (1–2 giorni)

    - Rendere **bilingue (IT/EN)** tutte le pagine “core”:
      - Architecture, Database, Operations, Development, Troubleshooting, Reference (config/endpoints)
    - Standardizzare stile: `Guides → Documentation Standard`
    - Aggiungere “Se cerchi…” come tabella all’inizio di ogni pagina core
    - Eliminare ogni warning `mkdocs --strict`

    ## Milestone 2 — Copertura funzionale (2–4 giorni)

    - Documentare in modo approfondito:
      - Difficulty scoring (`simulator/difficulty_scoring.py`) ✅ (prima iterazione)
      - Virtual patients (`simulator/virtual_patients.py`) ✅ (prima iterazione)
      - Prognosis (`simulator/prognosis.py`) ✅ (prima iterazione)
      - Regimen suggester (`simulator/regimen_suggester.py`) ✅ (prima iterazione)
    - Aggiungere esempi numerici e “worked examples” (input → output) per ogni modulo ✅ (base)
    - Aggiungere figure SVG dedicate (sensitivity, score curves, prognosis curves) ✅ (base)

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

=== "EN"
    This roadmap defines the next iterations to bring the documentation to a “software product” level.

    ## Milestone 1 — Stabilization and consistency (1–2 days)

    - Make all core pages **bilingual (IT/EN)**:
      - Architecture, Database, Operations, Development, Troubleshooting, Reference (config/endpoints)
    - Standardize style: `Guides → Documentation Standard`
    - Add a “Quick find / Se cerchi…” table at the start of each core page
    - Remove all warnings by keeping `mkdocs --strict` green

    ## Milestone 2 — Functional coverage (2–4 days)

    - Document deeply:
      - Difficulty scoring (`simulator/difficulty_scoring.py`) ✅ (first iteration)
      - Virtual patients (`simulator/virtual_patients.py`) ✅ (first iteration)
      - Prognosis (`simulator/prognosis.py`) ✅ (first iteration)
      - Regimen suggester (`simulator/regimen_suggester.py`) ✅ (first iteration)
    - Add numeric worked examples (inputs → outputs) for each module ✅ (baseline)
    - Add dedicated SVG figures (sensitivity, score curves, prognosis curves) ✅ (baseline)

    ## Milestone 3 — Full reference (3–6 days)

    - Endpoints: complete per-app table (URL → view → template → payload)
    - DB: query cookbook (ORM + SQL) for common use-cases
    - API: JSON payload schema for simulator/chemtools/docs_viewer
    - Operational runbooks (local/dev/docker, celery/redis, troubleshooting)

    ## Milestone 4 — Docs UX (ongoing)

    - Global index and aggressive cross-linking
    - Persona-based guided tours (clinician/dev/research) with time estimates
    - More Mermaid maps (components, sequences, dataflows)
    - Extended gallery with realistic demo cases

    ## Quality KPIs (measurable)

    - 0 warnings in `mkdocs build --strict`
    - 100% of core pages available in IT/EN
    - 100% of new features documented with:
      - overview + reference + example + diagram
    - Search coverage: minimum keywords for each module (e.g., `tumor_reduction`, `healthy_loss`, `cohort`, `twin`)
