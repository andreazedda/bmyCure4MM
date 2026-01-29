# Database Objects (Model ↔ Table ↔ File)

=== "IT"
    Questa pagina è pensata per “orientarsi” velocemente: dato un nome (es. `PatientTherapy`) vuoi sapere:

    - dove sta nel codice
    - come si chiama la tabella nel DB
    - quali FK/relazioni usa
    - come interrogarlo (ORM e SQL)

    !!! info "Fonte"
        - Definizioni dei model: file `*/models*.py`
        - Tabelle effettive: `db.sqlite3` (vedi anche `Guides → Database`)

    ## Se cerchi…

    | Voglio… | Vai a… |
    | --- | --- |
    | ER completo + DDL SQLite | `Guides → Database` |
    | parametri e output simulazione | `Reference → Simulator Parameters` |
    | endpoint | `Reference → Endpoints` |

    ## Mappa rapida (app → oggetti)

    ### `clinic`

    | Django Model | File | Tabella | Relazioni chiave |
    | --- | --- | --- | --- |
    | `Patient` | `clinic/models.py` | `clinic_patient` | `owner -> auth_user` |
    | `Assessment` | `clinic/models.py` | `clinic_assessment` | `patient -> clinic_patient` |
    | `Regimen` | `clinic/models.py` | `clinic_regimen` | usato da Clinic + Simulator |
    | `PatientTherapy` | `clinic/models.py` | `clinic_patienttherapy` | `patient -> clinic_patient`, `regimen -> clinic_regimen` |
    | `CytogeneticAbnormality` | `clinic/models.py` | `clinic_cytogeneticabnormality` |  |
    | `PatientCytogenetics` | `clinic/models.py` | `clinic_patientcytogenetics` | `patient -> clinic_patient`, `abnormality -> clinic_cytogeneticabnormality` |

    ### `simulator`

    | Django Model | File | Tabella | Relazioni chiave |
    | --- | --- | --- | --- |
    | `Scenario` | `simulator/models.py` | `simulator_scenario` | M2M con `clinic_regimen` |
    | `SimulationAttempt` | `simulator/models.py` | `simulator_simulationattempt` | `scenario -> simulator_scenario`, `selected_regimen -> clinic_regimen`, `user -> auth_user` |
    | `HelpArticle` | `simulator/models_help.py` | `simulator_helparticle` |  |
    | `Scenario.recommended_regimens` (M2M) | `simulator/models.py` | `simulator_scenario_recommended_regimens` | join `scenario_id`, `regimen_id` |

    ### `chemtools`

    | Django Model | File | Tabella | Relazioni chiave |
    | --- | --- | --- | --- |
    | `ChemJob` | `chemtools/models.py` | `chemtools_chemjob` | `user -> auth_user` |

    ### `docs_viewer`

    | Django Model | File | Tabella | Relazioni chiave |
    | --- | --- | --- | --- |
    | `DocumentView` | `docs_viewer/models.py` | `docs_viewer_documentview` | `user -> auth_user` |
    | `DocumentFeedback` | `docs_viewer/models.py` | `docs_viewer_documentfeedback` | unique (`path`, `user`) |

    ## Relazioni DB (mini-ER)

    ```mermaid
    erDiagram
      AUTH_USER ||--o{ CLINIC_PATIENT : owner
      CLINIC_PATIENT ||--o{ CLINIC_ASSESSMENT : has
      CLINIC_PATIENT ||--o{ CLINIC_PATIENTTHERAPY : has
      CLINIC_REGIMEN ||--o{ CLINIC_PATIENTTHERAPY : regimen

      CLINIC_PATIENT ||--o{ CLINIC_PATIENTCYTOGENETICS : has
      CLINIC_CYTOGENETICABNORMALITY ||--o{ CLINIC_PATIENTCYTOGENETICS : abnormality

      SIMULATOR_SCENARIO ||--o{ SIMULATOR_SIMULATIONATTEMPT : attempts
      AUTH_USER ||--o{ SIMULATOR_SIMULATIONATTEMPT : user
      CLINIC_REGIMEN ||--o{ SIMULATOR_SIMULATIONATTEMPT : selected_regimen

      SIMULATOR_SCENARIO ||--o{ SIMULATOR_SCENARIO_RECOMMENDED_REGIMENS : m2m
      CLINIC_REGIMEN ||--o{ SIMULATOR_SCENARIO_RECOMMENDED_REGIMENS : m2m
    ```

    ## Esempi pratici (ORM)

    ### Caricare un paziente con assessment e terapie

    ```python
    from clinic.models import Patient

    patient = Patient.objects.get(mrn="MRN001")
    assessments = patient.assessments.all()
    therapies = patient.therapies.select_related("regimen").all()
    ```

    ### Caricare un attempt e i risultati (JSON)

    ```python
    from simulator.models import SimulationAttempt

    attempt = SimulationAttempt.objects.select_related("scenario").get(pk=1)
    params = attempt.parameters
    results = attempt.results
    summary = attempt.results_summary
    ```

    ## Esempi pratici (SQL / SQLite)

    Mostrare le tabelle “app”:

    ```sql
    SELECT name
    FROM sqlite_master
    WHERE type='table' AND name NOT LIKE 'sqlite_%'
    ORDER BY name;
    ```

    Vedere struttura di una tabella:

    ```sql
    .schema clinic_patienttherapy
    ```

    ## “Non trovo un oggetto”: regole di naming

    Per default Django usa:

    - tabella: `<app>_<model_lower>`
    - FK: `<field>_id`
    - M2M: `<app>_<model>_<field>` (o simile, a seconda dei nomi)

    Se hai dubbi, vai su `Guides → Database` dove trovi anche il DDL completo della SQLite inclusa nel repo.

=== "EN"
    This page is for quick orientation: given a name (e.g., `PatientTherapy`), you want to know:

    - where it lives in the codebase
    - which DB table it maps to
    - key foreign keys / relations
    - how to query it (ORM and SQL)

    !!! info "Source"
        - Model definitions: `*/models*.py`
        - Effective tables: `db.sqlite3` (see also `Guides → Database`)

    ## Quick find

    | I want… | Go to… |
    | --- | --- |
    | full ER + SQLite DDL | `Guides → Database` |
    | simulator inputs/outputs | `Reference → Simulator Parameters` |
    | endpoints | `Reference → Endpoints` |

    ## Quick map (app → objects)

    ### `clinic`

    | Django Model | File | Table | Key relations |
    | --- | --- | --- | --- |
    | `Patient` | `clinic/models.py` | `clinic_patient` | `owner -> auth_user` |
    | `Assessment` | `clinic/models.py` | `clinic_assessment` | `patient -> clinic_patient` |
    | `Regimen` | `clinic/models.py` | `clinic_regimen` | used by Clinic + Simulator |
    | `PatientTherapy` | `clinic/models.py` | `clinic_patienttherapy` | `patient -> clinic_patient`, `regimen -> clinic_regimen` |
    | `CytogeneticAbnormality` | `clinic/models.py` | `clinic_cytogeneticabnormality` |  |
    | `PatientCytogenetics` | `clinic/models.py` | `clinic_patientcytogenetics` | `patient -> clinic_patient`, `abnormality -> clinic_cytogeneticabnormality` |

    ### `simulator`

    | Django Model | File | Table | Key relations |
    | --- | --- | --- | --- |
    | `Scenario` | `simulator/models.py` | `simulator_scenario` | M2M with `clinic_regimen` |
    | `SimulationAttempt` | `simulator/models.py` | `simulator_simulationattempt` | `scenario -> simulator_scenario`, `selected_regimen -> clinic_regimen`, `user -> auth_user` |
    | `HelpArticle` | `simulator/models_help.py` | `simulator_helparticle` |  |
    | `Scenario.recommended_regimens` (M2M) | `simulator/models.py` | `simulator_scenario_recommended_regimens` | join `scenario_id`, `regimen_id` |

    ### `chemtools`

    | Django Model | File | Table | Key relations |
    | --- | --- | --- | --- |
    | `ChemJob` | `chemtools/models.py` | `chemtools_chemjob` | `user -> auth_user` |

    ### `docs_viewer`

    | Django Model | File | Table | Key relations |
    | --- | --- | --- | --- |
    | `DocumentView` | `docs_viewer/models.py` | `docs_viewer_documentview` | `user -> auth_user` |
    | `DocumentFeedback` | `docs_viewer/models.py` | `docs_viewer_documentfeedback` | unique (`path`, `user`) |

    ## DB relations (mini-ER)

    ```mermaid
    erDiagram
      AUTH_USER ||--o{ CLINIC_PATIENT : owner
      CLINIC_PATIENT ||--o{ CLINIC_ASSESSMENT : has
      CLINIC_PATIENT ||--o{ CLINIC_PATIENTTHERAPY : has
      CLINIC_REGIMEN ||--o{ CLINIC_PATIENTTHERAPY : regimen

      CLINIC_PATIENT ||--o{ CLINIC_PATIENTCYTOGENETICS : has
      CLINIC_CYTOGENETICABNORMALITY ||--o{ CLINIC_PATIENTCYTOGENETICS : abnormality

      SIMULATOR_SCENARIO ||--o{ SIMULATOR_SIMULATIONATTEMPT : attempts
      AUTH_USER ||--o{ SIMULATOR_SIMULATIONATTEMPT : user
      CLINIC_REGIMEN ||--o{ SIMULATOR_SIMULATIONATTEMPT : selected_regimen

      SIMULATOR_SCENARIO ||--o{ SIMULATOR_SCENARIO_RECOMMENDED_REGIMENS : m2m
      CLINIC_REGIMEN ||--o{ SIMULATOR_SCENARIO_RECOMMENDED_REGIMENS : m2m
    ```

    ## Practical examples (ORM)

    ### Load a patient with assessments and therapies

    ```python
    from clinic.models import Patient

    patient = Patient.objects.get(mrn="MRN001")
    assessments = patient.assessments.all()
    therapies = patient.therapies.select_related("regimen").all()
    ```

    ### Load an attempt and its results (JSON)

    ```python
    from simulator.models import SimulationAttempt

    attempt = SimulationAttempt.objects.select_related("scenario").get(pk=1)
    params = attempt.parameters
    results = attempt.results
    summary = attempt.results_summary
    ```

    ## Practical examples (SQL / SQLite)

    List app tables:

    ```sql
    SELECT name
    FROM sqlite_master
    WHERE type='table' AND name NOT LIKE 'sqlite_%'
    ORDER BY name;
    ```

    Inspect table structure:

    ```sql
    .schema clinic_patienttherapy
    ```

    ## “I can’t find an object”: naming rules

    Django defaults:

    - table: `<app>_<model_lower>`
    - FK: `<field>_id`
    - M2M: `<app>_<model>_<field>` (or similar, depending on naming)

    If you’re unsure, go to `Guides → Database` for the full SQLite DDL committed in the repo.
