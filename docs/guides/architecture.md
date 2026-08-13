# Architecture

> **Governance boundary:** this is an `E1_research_prototype`. Components and
> data flows described below do not constitute clinical decision support,
> validated patient-specific prediction, or identified causal inference.

=== "IT"
    ## Cos’è bmyCure4MM

    bmyCure4MM è un’applicazione web **Django** orientata a:

    - **Clinic**: gestione pazienti/assessment/terapie (CRUD + timeline)
    - **Simulator**: scenari clinici e simulazioni (PK/PD, patient twin, reportistica)
    - **ChemTools**: tool di chem-informatics (job asincroni + output in media)
    - **Docs Viewer**: viewer interno per documentazione Markdown (analytics + feedback)

    La root app Django è `mmportal/`.

    ## Se cerchi…

    | Voglio… | Vai a… |
    | --- | --- |
    | schema DB e tabelle | `Guides → Database` |
    | mapping Model↔Table | `Reference → Database Objects` |
    | parametri simulazione | `Reference → Simulator Parameters` |
    | teoria modelli (ODE/PK/PD) | `Guides → Mathematical Models` |

    ## Terminologia

    - **App**: modulo Django (es. `clinic`, `simulator`).
    - **Artifact**: file generato (CSV/HTML/JSON) salvato in `media/`.
    - **HTMX**: aggiornamenti parziali della UI senza SPA.
    - **Attempt**: tentativo simulazione persistito (`SimulationAttempt`).

    ## Componenti (vista di insieme)

    ```mermaid
    flowchart LR
      subgraph Client
        B[Browser]
      end

      subgraph DjangoApp["Django (mmportal)"]
        URLS[URL routing]
        VIEWS[Views / DRF ViewSets]
        TPL[Templates]
        ORM[Django ORM]
      end

      subgraph Apps
        CLINIC[clinic]
        SIM[simulator]
        CHEM[chemtools]
        DV[docs_viewer]
      end

      subgraph Infra
        DB[(SQLite / Postgres)]
        REDIS[(Redis)]
        CELERY[Celery worker]
        FS[(media/ + static/)]
      end

      B -->|HTTP| URLS --> VIEWS
      VIEWS --> TPL
      VIEWS -->|ORM| ORM --> DB

      VIEWS --> CLINIC
      VIEWS --> SIM
      VIEWS --> CHEM
      VIEWS --> DV

      CHEM -->|enqueue| REDIS --> CELERY
      CELERY -->|read/write| DB
      CELERY -->|artifacts| FS
      SIM -->|plots/artifacts| FS
      TPL -->|static| FS
    ```

    !!! note "Frontend"
        La UI è principalmente server-rendered (Django templates) con interazioni incremental-update tramite **HTMX**.

    ## Apps (responsabilità)

    | App | Scopo | Model principali | Folder template |
    | --- | --- | --- | --- |
    | `clinic` | pazienti/assessment/terapie | `Patient`, `Assessment`, `PatientTherapy` | `clinic/templates/clinic/` |
    | `simulator` | scenari + simulazioni + twin | `Scenario`, `SimulationAttempt` | `simulator/templates/simulator/` |
    | `chemtools` | job chem-informatics | `ChemJob` | `chemtools/templates/chemtools/` |
    | `docs_viewer` | viewer doc interno | `DocumentView`, `DocumentFeedback` | `docs_viewer/templates/docs_viewer/` |

    ## Flussi principali

    ### Simulazione (Scenario → Attempt → Results)

    ```mermaid
    sequenceDiagram
      autonumber
      actor U as Utente
      participant W as Django view
      participant DB as DB
      participant M as Modello matematico
      participant FS as media/

      U->>W: Submit form (scenario + parametri)
      W->>DB: Crea/aggiorna SimulationAttempt
      W->>M: run_model(parameters)
      M-->>W: trajectory + summary + artifacts
      W->>DB: Salva results/results_summary/artifacts (JSON)
      W->>FS: (opzionale) salva HTML/plot in media/
      W-->>U: Render pagina risultati (HTMX o full page)
    ```

    ### ChemTools (job asincrono)

    ```mermaid
    sequenceDiagram
      autonumber
      actor U as Utente
      participant W as Django view
      participant DB as DB
      participant R as Redis (broker)
      participant C as Celery worker
      participant FS as media/

      U->>W: Avvia tool (es. similarity)
      W->>DB: Crea ChemJob (Queued)
      W->>R: enqueue task(job_id)
      C->>R: fetch task
      C->>DB: update_progress(...)
      C->>FS: scrive output (html/csv)
      C->>DB: aggiorna ChemJob (Completed/Failed)
      W-->>U: UI poll/refresh stato (HTMX)
    ```

    ## Dati e persistenza

    - Default database: `db.sqlite3` (sviluppo)
    - In produzione: consigliato Postgres

    ## Punti operativi

    - Static: `mmportal/static/` → `collectstatic` in produzione
    - Media: `media/` (artifact e upload)
    - Logs: `logs/` (Django + activity + Celery)

=== "EN"
    ## What is bmyCure4MM

    bmyCure4MM is a **Django** web application focused on:

    - **Clinic**: patient records, assessments, therapies
    - **Simulator**: scenarios + PK/PD simulations + patient twin
    - **ChemTools**: cheminformatics tools as background jobs
    - **Docs Viewer**: in-app Markdown documentation viewer

    Root Django project: `mmportal/`.

    ## Quick find

    | I want… | Go to… |
    | --- | --- |
    | DB schema and tables | `Guides → Database` |
    | Model↔Table mapping | `Reference → Database Objects` |
    | simulation parameters | `Reference → Simulator Parameters` |
    | ODE/PK/PD theory | `Guides → Mathematical Models` |

    ## Terminology

    - **App**: Django module (e.g., `clinic`, `simulator`).
    - **Artifact**: generated file (CSV/HTML/JSON) stored under `media/`.
    - **HTMX**: partial UI updates without a full SPA.
    - **Attempt**: persisted simulation attempt (`SimulationAttempt`).

    ## High-level components

    ```mermaid
    flowchart LR
      subgraph Client
        B[Browser]
      end

      subgraph DjangoApp["Django (mmportal)"]
        URLS[URL routing]
        VIEWS[Views / DRF ViewSets]
        TPL[Templates]
        ORM[Django ORM]
      end

      subgraph Apps
        CLINIC[clinic]
        SIM[simulator]
        CHEM[chemtools]
        DV[docs_viewer]
      end

      subgraph Infra
        DB[(SQLite / Postgres)]
        REDIS[(Redis)]
        CELERY[Celery worker]
        FS[(media/ + static/)]
      end

      B -->|HTTP| URLS --> VIEWS
      VIEWS --> TPL
      VIEWS -->|ORM| ORM --> DB

      VIEWS --> CLINIC
      VIEWS --> SIM
      VIEWS --> CHEM
      VIEWS --> DV

      CHEM -->|enqueue| REDIS --> CELERY
      CELERY -->|read/write| DB
      CELERY -->|artifacts| FS
      SIM -->|plots/artifacts| FS
      TPL -->|static| FS
    ```

    !!! note "Frontend"
        Mostly server-rendered templates with incremental updates via **HTMX**.

    ## Apps overview

    | App | Purpose | Main models | Templates |
    | --- | --- | --- | --- |
    | `clinic` | patients/assessments/therapies | `Patient`, `Assessment`, `PatientTherapy` | `clinic/templates/clinic/` |
    | `simulator` | scenarios + simulations + twin | `Scenario`, `SimulationAttempt` | `simulator/templates/simulator/` |
    | `chemtools` | background chem jobs | `ChemJob` | `chemtools/templates/chemtools/` |
    | `docs_viewer` | in-app docs viewer | `DocumentView`, `DocumentFeedback` | `docs_viewer/templates/docs_viewer/` |

    ## Core flows

    ```mermaid
    sequenceDiagram
      autonumber
      actor U as User
      participant W as Django view
      participant DB as DB
      participant M as Mathematical model
      participant FS as media/

      U->>W: Submit form (scenario + params)
      W->>DB: Create/update SimulationAttempt
      W->>M: run_model(parameters)
      M-->>W: trajectory + summary + artifacts
      W->>DB: Persist results/results_summary/artifacts (JSON)
      W->>FS: (optional) save HTML/plots in media/
      W-->>U: Render results (HTMX or full page)
    ```

    ## Persistence and ops

    - Default DB: SQLite (`db.sqlite3`)
    - Production: typically Postgres + `collectstatic`, background workers, logs.

    !!! tip "More details"
        For deeper engineering notes (patterns, security, performance, deployment), see `Guides → Architecture Deep Dive`.
