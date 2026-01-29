# Optimization Theory (Optimization Lab)

=== "IT"
    Questa sezione documenta **come** l’Optimization Lab costruisce e valuta le soluzioni (multi-obiettivo / Pareto).

    ## Se cerchi…

    | Voglio… | Vai a… |
    | --- | --- |
    | come funziona Optuna nel repo | `simulator/optim.py` |
    | mapping obiettivi/vincoli | `simulator/objectives.py` |
    | spazio di ricerca | `simulator/search_space.py` |
    | esempi visivi | `Guides → Simulations Gallery` |

    ## Obiettivi e vincoli

    Nel codice (`simulator/objectives.py`) gli obiettivi sono:

    | Nome | Definizione | Direzione |
    | --- | --- | --- |
    | Efficacy | `tumor_reduction` | massimizza |
    | Safety | `-healthy_loss` | massimizza (equivale a minimizzare loss) |
    | Exposure | `-sum(AUC)` | massimizza (equivale a minimizzare esposizione) |

    Vincolo “hard”:

    - `healthy_loss <= 0.25` (soluzione “feasible”)

    ```mermaid
    flowchart TD
      S[Trial params] --> R[Run model]
      R --> K[KPI summary]
      K --> C{Constraints OK?}
      C -- no --> P[Penalty -1e6]
      C -- yes --> O["Objectives (3D)"]
      O --> PF[Pareto front]
    ```

    ## Spazio di ricerca (search space)

    Definito in `simulator/search_space.py`:

    | Parametro | Tipo | Range/step | Note |
    | --- | --- | --- | --- |
    | `lenalidomide_dose` | float | 0..40 step 2.5 | dose |
    | `bortezomib_dose` | float | 0..1.6 step 0.1 | dose |
    | `daratumumab_dose` | float | 0..16 step 1.0 | dose |
    | `time_horizon` | int | 56..224 step 28 | 2–8 cicli |
    | `interaction_strength` | float | 0..0.15 step 0.01 | interazione |

    Manopole schedule (oggi principalmente metadata):

    - `len_on_days`, `bor_weekly`, `dara_interval`

    !!! note "Importante"
        La schedule “runtime” effettiva entra nel modello tramite `dose_functions` (preset YAML).

    ## Algoritmo (Optuna)

    `simulator/optim.py` usa:

    - sampler: `TPESampler(seed=..., multivariate=True)`
    - pruner: `MedianPruner(n_startup_trials=10, ...)`
    - directions: maximize/maximize/maximize

    ## Pareto front (non-dominated)

    Una soluzione \(A\) domina \(B\) se:

    - \(A\) è **>=** su tutti gli obiettivi
    - ed è **>** su almeno uno

    ```mermaid
    flowchart LR
      A[Soluzione A] --> D{Domina B?}
      B[Soluzione B] --> D
      D -->|tutti >= e uno >| YES[Dominated]
      D -->|altrimenti| NO[Non-dominated]
    ```

    ## Interpretazione pratica

    - a parità di `tumor_reduction`, preferisci `healthy_loss` più basso
    - evita miglioramenti marginali di efficacia con grande costo in tossicità
    - l’esposizione (AUC) aiuta a scegliere regimi “più parsimoniosi”

=== "EN"
    This section explains how the Optimization Lab evaluates solutions (multi-objective / Pareto).

    ## Quick find

    | I want… | Go to… |
    | --- | --- |
    | Optuna implementation | `simulator/optim.py` |
    | objectives/constraints mapping | `simulator/objectives.py` |
    | search space | `simulator/search_space.py` |
    | visual examples | `Guides → Simulations Gallery` |

    ## Objectives & constraints

    | Name | Definition | Direction |
    | --- | --- | --- |
    | Efficacy | `tumor_reduction` | maximize |
    | Safety | `-healthy_loss` | maximize (i.e., minimize loss) |
    | Exposure | `-sum(AUC)` | maximize (i.e., minimize exposure) |

    Hard constraint:

    - `healthy_loss <= 0.25` (feasible solution)

    ```mermaid
    flowchart TD
      S[Trial params] --> R[Run model]
      R --> K[KPI summary]
      K --> C{Constraints OK?}
      C -- no --> P[Penalty -1e6]
      C -- yes --> O["Objectives (3D)"]
      O --> PF[Pareto front]
    ```

    ## Search space

    | Parameter | Type | Range/step | Notes |
    | --- | --- | --- | --- |
    | `lenalidomide_dose` | float | 0..40 step 2.5 | dose |
    | `bortezomib_dose` | float | 0..1.6 step 0.1 | dose |
    | `daratumumab_dose` | float | 0..16 step 1.0 | dose |
    | `time_horizon` | int | 56..224 step 28 | 2–8 cycles |
    | `interaction_strength` | float | 0..0.15 step 0.01 | interaction |

    ## Pareto front

    A solution \(A\) dominates \(B\) if it is not worse on all objectives and better on at least one.
