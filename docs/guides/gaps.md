# Documentation Coverage (code vs docs)

=== "IT"
    Questa pagina esplicita **cosa c’è nel codice** ma non è ancora spiegato in modo esaustivo nella documentazione, così sai subito dove mancano pezzi e cosa documentare per primi.

    ## Gap principali (alto impatto)

    ### 1) Difficulty scoring & expected outcomes
    Nel codice esiste una pipeline di scoring/predizioni:

    - `simulator/difficulty_scoring.py` (scoring componenti, stime response/toxicity/survival)
    - `simulator/scenario_extensions.py` (mixin “mathematical” e calcolo difficulty)

    **Manca** nella doc:
    - formule/razionale delle componenti di score
    - mapping “campo scenario” → “score”
    - esempi con numeri (input → score → badge)

    ### 2) Virtual patients generator
    Nel codice:
    - `simulator/virtual_patients.py` (archetipi, generator, parametri)

    **Manca** nella doc:
    - come vengono generati gli archetipi (distribuzioni, range)
    - come si collega a Scenario/Twin e ai preset

    ### 3) Prognosis module
    Nel codice:
    - `simulator/prognosis.py` (PFS/OS stimati, modificatori)

    **Manca** nella doc:
    - definizione delle metriche e assunzioni
    - esempi (R-ISS I vs III, effetto citogenetica/età)

    ### 4) Regimen suggester
    Nel codice:
    - `simulator/regimen_suggester.py` (db regimi + motore suggerimenti)

    **Manca** nella doc:
    - regole di suggerimento (feature→regimi)
    - come si integra con UI Clinic/Simulator

    ## Gap medi (utile per dev)

    - API interne (HTMX endpoints) e payload: `simulator/views*.py`, `clinic/views.py`
    - gestione artifact e media: `SimulationAttempt.run_model` (già parzialmente coperto) + templates clinic
    - `docs_viewer`: whitelist e sicurezza (coperto da “Deep Dives”, ma manca quickstart operativo)

    ## Cosa abbiamo già coperto bene

    - schema DB + mappa Model↔Tabella: `Guides → Database`, `Reference → Database Objects`
    - modelli ODE/PK/PD + KPI + grafici toy: `Guides → Mathematical Models`
    - optimization theory e Pareto: `Guides → Optimization Theory`

=== "EN"
    This page lists what **exists in code** but is not yet fully explained in the docs, so you can prioritize documentation work.

    ## Main gaps (high impact)

    ### 1) Difficulty scoring & expected outcomes
    In code:
    - `simulator/difficulty_scoring.py` (component scores + response/toxicity/survival estimation)
    - `simulator/scenario_extensions.py` (mathematical mixin + difficulty computation)

    **Missing in docs:**
    - the math/rationale behind each component
    - field mapping “Scenario → score”
    - worked examples (inputs → score → badges)

    ### 2) Virtual patients generator
    In code:
    - `simulator/virtual_patients.py` (archetypes + generator)

    **Missing in docs:**
    - how archetypes are generated (ranges/distributions)
    - how it connects to Scenario/Twin and presets

    ### 3) Prognosis module
    In code:
    - `simulator/prognosis.py` (estimated PFS/OS + modifiers)

    **Missing in docs:**
    - metric definitions and assumptions
    - examples (R-ISS I vs III, cytogenetics/age effect)

    ### 4) Regimen suggester
    In code:
    - `simulator/regimen_suggester.py` (regimen DB + suggestion engine)

    **Missing in docs:**
    - suggestion rules (features → regimens)
    - how it integrates into Clinic/Simulator UI

    ## Medium gaps (dev focused)

    - internal/HTMX endpoints and payloads: `simulator/views*.py`, `clinic/views.py`
    - artifact/media lifecycle: `SimulationAttempt.run_model` + clinic templates
    - `docs_viewer`: operational quickstart (security is covered in deep dives)

    ## Already well covered

    - DB schema + Model↔Table mapping: `Guides → Database`, `Reference → Database Objects`
    - ODE/PK/PD model + KPIs + toy figures: `Guides → Mathematical Models`
    - optimization + Pareto: `Guides → Optimization Theory`

