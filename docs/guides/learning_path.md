# Learning Path (find info in 30 seconds)

=== "IT"
    Questa pagina è il “**GPS**” della documentazione: cosa leggere, in che ordine, e dove trovare rapidamente le informazioni.

    ## Se cerchi… vai qui

    | Voglio… | Vai a… |
    | --- | --- |
    | capire cosa fa il sistema | `Guides → Architecture` |
    | vedere **tabelle DB** e relazioni | `Guides → Database` |
    | trovare **un oggetto DB** (Model ↔ tabella ↔ file) | `Reference → Database Objects` |
    | capire il **modello ODE / PK/PD** e KPI | `Guides → Mathematical Models` |
    | capire l’**Optimization Lab** (Pareto, obiettivi) | `Guides → Optimization Theory` |
    | capire difficulty / prognosis / virtual patients | `Guides → Difficulty Scoring`, `Guides → Prognosis`, `Guides → Virtual Patients` |
    | vedere regole e soglie “auditabili” | `Guides → Decision Algorithm` |
    | avviare local/dev/docker | `Guides → Operations` + `Guides → Development` |
    | consultare documenti lunghi/tecnici | `Deep Dives` |

    ## Percorso “Clinico” (30–45 min)

    1. `IT → Simulatore` (workflow + KPI)
    2. `Guides → Mathematical Models` (solo sezioni “KPI” e “Patient Twin”)
    3. `Guides → Decision Algorithm` (soglie + regole principali)
    4. `IT → Optimization Lab` (come leggere la Pareto table)

    Se vuoi un percorso “hands-on” con dati demo: `Deep Dives → Learning (IT) → Guida apprendimento`.

    ## Percorso “Developer” (45–90 min)

    1. `Guides → Architecture`
    2. `Reference → Configuration`
    3. `Reference → Database Objects`
    4. `Guides → Database` (DDL e relazioni)
    5. `Guides → Mathematical Models` (pipeline/solver/outputs)

    ## Percorso “Ricerca/Modeling” (60–120 min)

    1. `Guides → Mathematical Models` (tutto)
    2. `Guides → Optimization Theory`
    3. `Deep Dives → Mathematical Models` (documento lungo già presente nel repo)

    ## Come cercare bene (MkDocs)

    La search di MkDocs Material è la via più veloce:

    - cerca per `PatientTherapy`, `clinic_patienttherapy`, `SimulationAttempt`, `simulator_simulationattempt`
    - cerca per KPI: `tumor_reduction`, `healthy_loss`, `AUC`, `time_to_recurrence`

    !!! tip "Se vedi ‘Syntax error in text’"
        È un errore di parsing Mermaid in una pagina. Vai su `Home` e verifica che *diagramma* e *formula* si vedano; se si vedono lì ma non altrove, il problema è localizzato a un blocco Mermaid specifico.

=== "EN"
    This page is the documentation “**GPS**”: what to read, in what order, and where to find things quickly.

    ## Quick find

    | I want… | Go to… |
    | --- | --- |
    | understand what the system does | `Guides → Architecture` |
    | see **DB tables** and relations | `Guides → Database` |
    | find a **DB object** (Model ↔ table ↔ file) | `Reference → Database Objects` |
    | understand the **ODE / PK/PD** model and KPIs | `Guides → Mathematical Models` |
    | understand the **Optimization Lab** (Pareto, objectives) | `Guides → Optimization Theory` |
    | understand difficulty / prognosis / virtual patients | `Guides → Difficulty Scoring`, `Guides → Prognosis`, `Guides → Virtual Patients` |
    | see auditable rules and thresholds | `Guides → Decision Algorithm` |
    | run local/dev/docker | `Guides → Operations` + `Guides → Development` |
    | read long/technical documents | `Deep Dives` |

    ## Clinician path (30–45 min)

    1. `EN → Simulator` (workflow + KPIs)
    2. `Guides → Mathematical Models` (KPI + Patient Twin sections first)
    3. `Guides → Decision Algorithm` (main thresholds/rules)
    4. `EN → Optimization Lab` (how to read the Pareto table)

    ## Developer path (45–90 min)

    1. `Guides → Architecture`
    2. `Reference → Configuration`
    3. `Reference → Database Objects`
    4. `Guides → Database` (DDL + relations)
    5. `Guides → Mathematical Models` (pipeline/solver/outputs)

    ## Research/modeling path (60–120 min)

    1. `Guides → Mathematical Models` (full)
    2. `Guides → Optimization Theory`
    3. `Deep Dives → Mathematical Models` (long form doc already in the repo)

    ## How to search effectively (MkDocs)

    MkDocs Material search is the fastest way:

    - search by `PatientTherapy`, `clinic_patienttherapy`, `SimulationAttempt`, `simulator_simulationattempt`
    - search by KPIs: `tumor_reduction`, `healthy_loss`, `AUC`, `time_to_recurrence`
