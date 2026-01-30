# Decision Algorithm (transparency and audit)

=== "IT"
    Questa guida descrive l’algoritmo decisionale dichiarativo (rule-based) esposto in `simulator/decision_algorithm.py`.

    ## Se cerchi…

    | Voglio… | Vai a… |
    | --- | --- |
    | file sorgente algoritmo | `simulator/decision_algorithm.py` |
    | definizione KPI | `Guides → Mathematical Models` |
    | endpoint/UI che mostra badge | `Reference → Endpoints` |

    ## Terminologia

    - **Threshold**: classificazione di un KPI in una classe (good/moderate/poor).
    - **Rule**: condizione logica che produce un’azione e un rationale.
    - **Badge**: etichetta UI derivata da threshold/rule.

    ## Perché esiste

    L’obiettivo è rendere:

    - **trasparenti** soglie e regole
    - **auditabili** le raccomandazioni
    - **versionabili** i cambi di policy

    ## Struttura dati

    Il file espone un dizionario `DECISION_ALGORITHM` con:

    - `version` e `last_updated`
    - `thresholds`: classificazioni (efficacy/toxicity/durability)
    - `decision_rules`: regole attivabili da trigger condition
    - `risk_stratification`: schema R-ISS (descrittivo)
    - `high_risk_cytogenetics`: impatto e management per marker

    ```mermaid
    flowchart TD
      S[KPI summary] --> T[Threshold classification]
      S --> R[Evaluate decision rules]
      T --> O[Badges/labels]
      R --> A[Actions + rationale + evidence]
      O --> UI[UI rendering]
      A --> UI
    ```

    ## Soglie principali (KPI → badge)

    ### Efficacy

    Basato su `tumor_reduction`:

    - good: `>= 0.50`
    - moderate: `0.30–0.49`
    - poor: `< 0.30`
    - failure: `< 0` (crescita tumorale)

    ### Toxicity

    Basato su `healthy_loss`:

    - acceptable: `< 0.20`
    - borderline: `0.20–0.29`
    - high: `>= 0.30`

    ### Durability

    Basato su `time_to_recurrence` (giorni):

    - good: `>= 180`
    - moderate: `90–179`
    - poor: `< 90`

    ## Regole (decision_rules)

    Ogni regola contiene:

    - `trigger_condition` (espressione logica su KPI)
    - `priority`
    - `action` (EN/IT)
    - `rationale` (EN/IT)
    - `evidence` (lista riferimenti testuali)
    - (opzionale) `suggested_alternatives`

    Esempio di flusso:

    ```mermaid
    sequenceDiagram
      autonumber
      participant K as KPI summary
      participant E as Evaluator
      participant UI as UI
      K->>E: tumor_reduction, healthy_loss, time_to_recurrence
      E->>E: apply thresholds
      E->>E: evaluate trigger conditions
      E-->>UI: badges + suggested actions
    ```

    ## Versioning e date (importante)

    `ALGORITHM_UPDATED` è una data ISO (es. `2026-01-17`) e va aggiornata quando:

    - cambiano soglie
    - si aggiungono regole
    - si cambia wording clinico rilevante

    !!! tip "Buona pratica"
        Tratta l’algoritmo come “policy as code”: review, changelog, test automatici sulle condizioni critiche.

=== "EN"
    This guide describes the declarative (rule-based) decision algorithm exposed in `simulator/decision_algorithm.py`.

    ## Quick find

    | I want… | Go to… |
    | --- | --- |
    | algorithm source file | `simulator/decision_algorithm.py` |
    | KPI definitions | `Guides → Mathematical Models` |
    | endpoint/UI that shows badges | `Reference → Endpoints` |

    ## Terminology

    - **Threshold**: KPI classification into a class (good/moderate/poor).
    - **Rule**: logical condition producing an action and rationale.
    - **Badge**: UI label derived from thresholds/rules.

    ## Why it exists

    The goal is to make:

    - thresholds and rules **transparent**
    - recommendations **auditable**
    - policy changes **versionable**

    ## Data structure

    The file exposes a `DECISION_ALGORITHM` dictionary with:

    - `version` and `last_updated`
    - `thresholds`: classifications (efficacy/toxicity/durability)
    - `decision_rules`: rules driven by trigger conditions
    - `risk_stratification`: R-ISS schema (descriptive)
    - `high_risk_cytogenetics`: impact and management notes for markers

    ```mermaid
    flowchart TD
      S[KPI summary] --> T[Threshold classification]
      S --> R[Evaluate decision rules]
      T --> O[Badges/labels]
      R --> A[Actions + rationale + evidence]
      O --> UI[UI rendering]
      A --> UI
    ```

    ## Main thresholds (KPI → badge)

    ### Efficacy

    Based on `tumor_reduction`:

    - good: `>= 0.50`
    - moderate: `0.30–0.49`
    - poor: `< 0.30`
    - failure: `< 0` (tumor growth)

    ### Toxicity

    Based on `healthy_loss`:

    - acceptable: `< 0.20`
    - borderline: `0.20–0.29`
    - high: `>= 0.30`

    ### Durability

    Based on `time_to_recurrence` (days):

    - good: `>= 180`
    - moderate: `90–179`
    - poor: `< 90`

    ## Rules (decision_rules)

    Each rule includes:

    - `trigger_condition` (logical expression on KPIs)
    - `priority`
    - `action` (EN/IT)
    - `rationale` (EN/IT)
    - `evidence` (text references list)
    - (optional) `suggested_alternatives`

    Flow example:

    ```mermaid
    sequenceDiagram
      autonumber
      participant K as KPI summary
      participant E as Evaluator
      participant UI as UI
      K->>E: tumor_reduction, healthy_loss, time_to_recurrence
      E->>E: apply thresholds
      E->>E: evaluate trigger conditions
      E-->>UI: badges + suggested actions
    ```

    ## Versioning and dates

    `ALGORITHM_UPDATED` is an ISO date (e.g., `2026-01-17`) and should be updated when:

    - thresholds change
    - rules are added/removed
    - clinically relevant wording changes

    !!! tip "Best practice"
        Treat the algorithm as “policy as code”: review, changelog, and automated tests for critical conditions.
