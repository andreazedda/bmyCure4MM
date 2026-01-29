# Glossary (terminology)

=== "IT"
    Questa pagina raccoglie la terminologia più usata nella piattaforma, con definizioni **leggibili da umani** e riferimenti a dove si vede in UI/codice.

    ## Clinica (MM)

    - **MM (Multiple Myeloma)**: neoplasia delle plasmacellule; nel simulatore è modellata come carico tumorale \(T(t)\).
    - **CRAB**: *Calcium, Renal, Anemia, Bone* (manifestazioni tipiche).
    - **R-ISS**: staging prognostico; nella piattaforma compare in `Assessment` e negli scenari.
    - **FLC ratio**: rapporto catene leggere; usato come indicatore di attività di malattia (Twin/score).
    - **LDH**: marker di aggressività/turnover; entra in risk/score.
    - **β2-microglobulin**: marker prognostico; entra in R-ISS e rischio.

    ## Modelli matematici

    - **PK (Pharmacokinetics)**: come \(C(t)\) cambia per input ed eliminazione.
    - **PD (Pharmacodynamics)**: come \(C(t)\) si traduce in effetto \(E(C)\).
    - **Emax**: modello saturante dell’effetto; vedi `Guides → Mathematical Models`.
    - **EC50**: concentrazione che produce 50% dell’effetto massimo.
    - **K (carrying capacity)**: “limite superiore” di crescita in modelli logistici.
    - **Half-life**: tempo di dimezzamento; determina \(k_{elim}\).
    - **AUC**: area sotto la curva di concentrazione; proxy di esposizione complessiva.
    - **KPI**: metriche riassuntive (es. `tumor_reduction`, `healthy_loss`).

    ## Difficulty / prognosis (stima)

    - **Difficulty score (0–100)**: score composito che stima “quanto è difficile” trattare uno scenario; vedi `Guides → Difficulty Scoring`.
    - **Frailty (frailty score)**: proxy di tollerabilità del paziente (età, ECOG, comorbidità, funzione renale, albumina).
    - **PFS / OS**: *Progression-Free Survival* / *Overall Survival*; vedi `Guides → Prognosis`.

    ## Ottimizzazione

    - **Multi-objective optimization**: ottimizzare più obiettivi contemporaneamente (efficacia/sicurezza/esposizione).
    - **Pareto front**: insieme di soluzioni non dominate (nessuna è “meglio in tutto” delle altre).
    - **Trial**: una proposta di parametri (dosi/horizon) valutata dal modello.
    - **Constraint**: vincolo “hard” (es. `healthy_loss <= 0.25`).

    ## Software / dati

    - **Django model**: classi Python ORM (es. `clinic.Patient`) che mappano tabelle DB.
    - **Migration**: versione dello schema DB (Django migrations).
    - **Artifact**: file generato (CSV/HTML/JSON) salvato sotto `MEDIA_ROOT`.
    - **JSONField**: campo strutturato; in SQLite è serializzato in `TEXT` con check JSON.

    ## Dove trovare “la cosa giusta”

    - Oggetti DB: `Reference → Database Objects`
    - Parametri simulatore: `Reference → Simulator Parameters`
    - Modelli e KPI: `Guides → Mathematical Models`
    - Ottimizzazione: `Guides → Optimization Theory`

=== "EN"
    This page collects the most common terms used in the platform, with **human-readable** definitions and pointers to where they appear in UI/code.

    ## Clinical (MM)

    - **MM (Multiple Myeloma)**: plasma cell malignancy; in the simulator it is modeled as tumor burden \(T(t)\).
    - **CRAB**: *Calcium, Renal, Anemia, Bone* (typical clinical manifestations).
    - **R-ISS**: prognostic staging; appears in `Assessment` and scenario logic.
    - **FLC ratio**: free light chain ratio; used as a disease activity marker (Twin/score).
    - **LDH**: aggressiveness/turnover marker; affects risk/score.
    - **β2-microglobulin**: prognostic marker; used in R-ISS and risk.

    ## Mathematical models

    - **PK (Pharmacokinetics)**: how concentration \(C(t)\) evolves with input and elimination.
    - **PD (Pharmacodynamics)**: how \(C(t)\) maps to an effect \(E(C)\).
    - **Emax**: saturating effect model; see `Guides → Mathematical Models`.
    - **EC50**: concentration that yields 50% of the maximum effect.
    - **K (carrying capacity)**: upper growth limit in logistic-like models.
    - **Half-life**: time to halve; determines \(k_{elim}\).
    - **AUC**: area under the concentration curve; proxy for total exposure.
    - **KPI**: summary metrics (e.g., `tumor_reduction`, `healthy_loss`).

    ## Difficulty / prognosis (estimation)

    - **Difficulty score (0–100)**: composite score estimating how hard a scenario is to treat; see `Guides → Difficulty Scoring`.
    - **Frailty (frailty score)**: tolerance proxy (age, ECOG, comorbidities, renal function, albumin).
    - **PFS / OS**: *Progression-Free Survival* / *Overall Survival*; see `Guides → Prognosis`.

    ## Optimization

    - **Multi-objective optimization**: optimize multiple objectives at once (efficacy/safety/exposure).
    - **Pareto front**: set of non-dominated solutions (no single solution is “best in everything”).
    - **Trial**: a candidate parameter vector (doses/horizon) evaluated by the model.
    - **Constraint**: a hard validity condition (e.g., `healthy_loss <= 0.25`).

    ## Software / data

    - **Django model**: ORM Python classes (e.g., `clinic.Patient`) mapping to DB tables.
    - **Migration**: schema versioning (Django migrations).
    - **Artifact**: generated file (CSV/HTML/JSON) stored under `MEDIA_ROOT`.
    - **JSONField**: structured field; in SQLite it is stored as `TEXT` with JSON checks.

    ## Where to find the right thing

    - DB objects: `Reference → Database Objects`
    - Simulator parameters: `Reference → Simulator Parameters`
    - Models and KPIs: `Guides → Mathematical Models`
    - Optimization: `Guides → Optimization Theory`
