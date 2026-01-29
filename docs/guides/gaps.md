# Documentation Coverage (code vs docs)

=== "IT"
    Questa pagina esplicita **cosa c’è nel codice** ma non è ancora spiegato in modo esaustivo nella documentazione, così sai subito dove mancano pezzi e cosa documentare per primi.

    ## Gap principali (alto impatto)

    ### 1) Coverage “quasi completa”: difficulty / prognosis / virtual patients / regimens

    Ora esistono guide dedicate:

    - `Guides → Difficulty Scoring` (da `simulator/difficulty_scoring.py`)
    - `Guides → Prognosis` (da `simulator/prognosis.py`)
    - `Guides → Virtual Patients` (da `simulator/virtual_patients.py`)
    - `Guides → Regimen Suggester` (da `simulator/regimen_suggester.py`)

    **Restano da migliorare** (se vuoi “perfetto”):
    - collegamenti diretti UI → funzione (view/template) per ogni output
    - esempi end-to-end con dati reali/dummy (JSON) salvati come artifact
    - validazione/benchmark delle assunzioni (documento separato “evidence”)

    ### 2) Decision support “UI mapping”

    Nel codice esistono molte regole “auditable” distribuite tra:

    - `clinic/` (input clinici + vincoli)
    - `simulator/` (validazione parametri + twin + scoring)

    **Manca** nella doc:
    - una tabella “UI field → Model field → file → regola → messaggio badge”
    - esempi screenshot-first (per PM/clinici) per ogni regola principale

    ## Gap medi (utile per dev)

    - API interne (HTMX endpoints) e payload: `simulator/views*.py`, `clinic/views.py`
    - gestione artifact e media: `SimulationAttempt.run_model` (già parzialmente coperto) + templates clinic
    - `docs_viewer`: whitelist e sicurezza (coperto da “Deep Dives”, ma manca quickstart operativo)

    ## Cosa abbiamo già coperto bene

    - schema DB + mappa Model↔Tabella: `Guides → Database`, `Reference → Database Objects`
    - modelli ODE/PK/PD + KPI + grafici toy: `Guides → Mathematical Models`
    - optimization theory e Pareto: `Guides → Optimization Theory`
    - difficulty/prognosis/virtual patients/regimens (base + teoria): vedi guide dedicate

=== "EN"
    This page lists what **exists in code** but is not yet fully explained in the docs, so you can prioritize documentation work.

    ## Main gaps (high impact)

    ### 1) “Mostly covered”: difficulty / prognosis / virtual patients / regimens

    Dedicated guides now exist:

    - `Guides → Difficulty Scoring` (from `simulator/difficulty_scoring.py`)
    - `Guides → Prognosis` (from `simulator/prognosis.py`)
    - `Guides → Virtual Patients` (from `simulator/virtual_patients.py`)
    - `Guides → Regimen Suggester` (from `simulator/regimen_suggester.py`)

    **Still worth improving** (if you want it “perfect”):
    - direct UI → function (view/template) links for each output
    - end-to-end worked examples with saved artifacts (JSON)
    - evidence/assumption validation notes (separate “evidence” doc)

    ### 2) Decision support “UI mapping”

    Many auditable rules are spread across:

    - `clinic/` (clinical inputs and constraints)
    - `simulator/` (parameter validation + twin + scoring)

    **Missing in docs:**
    - a single table: “UI field → Model field → file → rule → badge message”
    - screenshot-first examples (PM/clinician friendly) for each major rule

    ## Medium gaps (dev focused)

    - internal/HTMX endpoints and payloads: `simulator/views*.py`, `clinic/views.py`
    - artifact/media lifecycle: `SimulationAttempt.run_model` + clinic templates
    - `docs_viewer`: operational quickstart (security is covered in deep dives)

    ## Already well covered

    - DB schema + Model↔Table mapping: `Guides → Database`, `Reference → Database Objects`
    - ODE/PK/PD model + KPIs + toy figures: `Guides → Mathematical Models`
    - optimization + Pareto: `Guides → Optimization Theory`
    - difficulty/prognosis/virtual patients/regimens (baseline theory): see dedicated guides
