# Simulator Parameters (reference)

=== "IT"
    Questa pagina è una reference “tecnica” per capire **tutti i parametri** che entrano nella simulazione e come vengono interpretati/validati prima di arrivare al solver.

    !!! info "Fonte nel codice"
        La conversione centrale è in `simulator/models.py` → `SimulationAttempt._resolve_solver_inputs(...)`.

    ## Se cerchi…

    | Voglio… | Vai a… |
    | --- | --- |
    | teoria modelli ODE/PK/PD | `Guides → Mathematical Models` |
    | KPI e soglie | `Guides → Decision Algorithm` |
    | esempi visivi | `Guides → Simulations Gallery` |

    ## Terminologia

    - **resolved_params**: dizionario risultante da form + Twin + default; può contenere stringhe tipo `"auto"`.
    - **solver_inputs**: dizionario “sanificato” (solo numeri), pronto per `MathematicalModel`.
    - **AUTO sentinel**: valori come `"auto"`, `""`, ecc. che vengono sostituiti da default.

    ## Parametri principali (ingresso solver)

    | Key (parameters JSON) | Tipo | Default | Unità | Effetto |
    | --- | --- | --- | --- | --- |
    | `baseline_tumor_cells` | float | `1.0e9` | cells (proxy) | \(T(0)\) |
    | `baseline_healthy_cells` | float | `5.0e11` | cells (proxy) | \(H(0)\) |
    | `time_horizon` | float | `90` | days | durata simulazione |
    | `tumor_growth_rate` | float | `0.023` | 1/day | crescita tumorale |
    | `healthy_growth_rate` | float | `0.015` | 1/day | recupero “healthy” |
    | `lenalidomide_dose` | float | `25.0` | mg (UI) | dose totale usata per schedule/PK |
    | `bortezomib_dose` | float | `1.3` | mg/m² (UI) | dose totale usata per schedule/PK |
    | `daratumumab_dose` | float | `16.0` | mg/kg (UI) | dose totale usata per schedule/PK |
    | `interaction_strength` | float | `0.05` | adim. | scala l’interazione (matrice identità) |
    | `immune_compromise_index` | float | `1.0` | adim. | aumenta tossicità su \(H(t)\) |
    | `carrying_capacity_tumor` | float\|None | `None` | cells | override \(K_T\) |
    | `carrying_capacity_healthy` | float\|None | `None` | cells | override \(K_H\) |

    !!! note "Unità"
        Le unità cliniche di dose (mg, mg/m², mg/kg) sono *presentazione UI*: nel runtime vengono trattate come quantità numeriche coerenti internamente.

    ## Parametri qualitativi (guided choices)

    Questi non entrano direttamente nel solver: applicano piccoli **moltiplicatori** ai parametri già risolti.

    | Key | Valori | Dove agisce | Effetto |
    | --- | --- | --- | --- |
    | `guided_tumor_aggressiveness` | `lower/typical/higher` | `tumor_growth_rate` | ~0.85 / 1.0 / 1.15 (clamped) |
    | `guided_immune_status` | `better/typical/worse` | `immune_compromise_index` | ~0.9 / 1.0 / 1.1 (clamped) |

    ## Patient Twin (parametri derivati)

    Quando `use_twin=true` e `twin_biology_mode=auto`, la simulazione può iniettare (solo se lasciati “auto”):

    - `tumor_growth_rate`, `healthy_growth_rate`
    - `carrying_capacity_tumor`, `carrying_capacity_healthy`
    - `immune_compromise_index`

    ## Coorte / incertezza

    | Key | Tipo | Default | Effetto |
    | --- | --- | --- | --- |
    | `cohort_size` | int | `1` | se >1, calcola bande p05/p95 e statistiche |
    | `seed` | int\|None | `None` | riproducibilità della coorte |

    ## Output salvati (dove guardare)

    La simulazione salva tre JSON sul record `SimulationAttempt`:

    ### `results` (link a file / artifact)

    - `csv`: path su `MEDIA_URL` a CSV dei time-series
    - `plot`: path su `MEDIA_URL` a HTML Plotly (iframe-friendly)
    - `generated_at`: ISO timestamp
    - (opzionale) `twin_params.json`: payload twin serializzato

    ### `results_summary` (KPI e metriche)

    Contiene:

    - `tumor_reduction`, `healthy_loss`
    - `nadir`, `milestones`, `auc`, `durability_index`, `time_to_recurrence`
    - (opzionale) `cohort` con statistiche (`mean/p05/p95`, ecc.)

    ### `artifacts`

    Copie “senza timestamp” utili per UI e audit; include anche `seed` se presente.

    ## Errori comuni (e perché)

    - **ValueError “Invalid solver input …”**: un campo non numerico non è “auto” e non è convertibile a float.
    - **time_span end must be greater than start**: `time_horizon` <= 0.
    - **ODE solver failed**: parametri troppo estremi (instabilità numerica) o schedule/input non coerente.

=== "EN"
    This is a technical reference for **all simulation parameters** and how they are validated before reaching the solver.

    !!! info "Code source"
        The central conversion layer is `simulator/models.py` → `SimulationAttempt._resolve_solver_inputs(...)`.

    ## Quick find

    | I want… | Go to… |
    | --- | --- |
    | ODE/PK/PD theory | `Guides → Mathematical Models` |
    | KPIs and thresholds | `Guides → Decision Algorithm` |
    | visual examples | `Guides → Simulations Gallery` |

    ## Terminology

    - **resolved_params**: form + Twin + defaults; may contain `"auto"` strings.
    - **solver_inputs**: numeric-only payload used to instantiate `MathematicalModel`.
    - **AUTO sentinel**: values like `"auto"`/`""` that get replaced by defaults.

    ## Main parameters (solver inputs)

    | Key (parameters JSON) | Type | Default | Unit | Meaning |
    | --- | --- | --- | --- | --- |
    | `baseline_tumor_cells` | float | `1.0e9` | cells (proxy) | \(T(0)\) |
    | `baseline_healthy_cells` | float | `5.0e11` | cells (proxy) | \(H(0)\) |
    | `time_horizon` | float | `90` | days | simulation duration |
    | `tumor_growth_rate` | float | `0.023` | 1/day | tumor growth |
    | `healthy_growth_rate` | float | `0.015` | 1/day | healthy recovery |
    | `lenalidomide_dose` | float | `25.0` | UI unit | numeric dose used by PK/schedule |
    | `bortezomib_dose` | float | `1.3` | UI unit | numeric dose used by PK/schedule |
    | `daratumumab_dose` | float | `16.0` | UI unit | numeric dose used by PK/schedule |
    | `interaction_strength` | float | `0.05` | unitless | scales the interaction term |
    | `immune_compromise_index` | float | `1.0` | unitless | increases toxicity on \(H(t)\) |
    | `carrying_capacity_tumor` | float\|None | `None` | cells | override \(K_T\) |
    | `carrying_capacity_healthy` | float\|None | `None` | cells | override \(K_H\) |

    ## Qualitative knobs (guided choices)

    | Key | Values | Acts on | Effect |
    | --- | --- | --- | --- |
    | `guided_tumor_aggressiveness` | `lower/typical/higher` | `tumor_growth_rate` | ~0.85 / 1.0 / 1.15 (clamped) |
    | `guided_immune_status` | `better/typical/worse` | `immune_compromise_index` | ~0.9 / 1.0 / 1.1 (clamped) |

    ## Patient Twin (derived biology)

    When `use_twin=true` and `twin_biology_mode=auto`, the Twin may inject (only if left “auto”):

    - `tumor_growth_rate`, `healthy_growth_rate`
    - `carrying_capacity_tumor`, `carrying_capacity_healthy`
    - `immune_compromise_index`

    ## Cohort / uncertainty

    | Key | Type | Default | Meaning |
    | --- | --- | --- | --- |
    | `cohort_size` | int | `1` | if >1, computes p05/p95 bands and summary stats |
    | `seed` | int\|None | `None` | reproducibility |

    ## Stored outputs

    The attempt persists:

    - `results`: artifact links (CSV/Plotly HTML, optional twin JSON)
    - `results_summary`: KPIs and cohort stats
    - `artifacts`: UI/audit-friendly subset
