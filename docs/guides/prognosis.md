# Prognosis (PFS/OS estimation)

=== "IT"
    Questa pagina spiega la logica di stima prognostica usata dal progetto.

    Fonte principale: `simulator/prognosis.py`.

    !!! warning "Disclaimer"
        Le stime PFS/OS sono **statistiche** e **educational**. Non sostituiscono judgement clinico.

    ## Se cerchi…

    | Voglio… | Vai a… |
    | --- | --- |
    | funzione principale di stima | `simulator/prognosis.py` → `estimate_prognosis(...)` |
    | citogenetica / hazard ratios | `simulator/prognosis.py` → `CYTOGENETIC_HAZARD_RATIOS` |
    | relazione con difficulty score | `Guides → Difficulty Scoring` |
    | staging clinico (R-ISS) | `Guides → Database` (`Assessment.r_iss`) |

    ## Terminologia

    - **PFS / OS**: progression-free survival / overall survival.
    - **HR (hazard ratio)**: HR>1 peggiora la prognosi.
    - **Baseline**: valori iniziali prima dei modificatori.

    ## Obiettivo del modulo

    Il modulo produce un oggetto `PrognosisEstimate` con:

    - mediane PFS/OS (mesi)
    - probabilità PFS/OS a 12/24/36 mesi (e OS a 60 mesi)
    - categoria rischio (`standard/intermediate/high/very_high`)
    - confidenza (0–1) + lista dei modificatori applicati

    ## Pipeline logica (high-level)

    ```mermaid
    flowchart TD
      A[R-ISS baseline] --> B[Apply cytogenetic HR]
      B --> C[Apply age/ECOG modifiers]
      C --> D[Apply response + MRD modifiers]
      D --> E[Apply line-of-therapy penalty]
      E --> F[PrognosisEstimate]
    ```

    ## Baseline (R-ISS)

    La baseline è una tabella per R-ISS (I/II/III) con mediane e probabilità a timepoint.

    | R-ISS | median PFS (mo) | median OS (mo) | OS 5y | risk |
    | --- | ---: | ---: | ---: | --- |
    | I | 66 | 120 | 0.75 | standard |
    | II | 42 | 83 | 0.55 | intermediate |
    | III | 29 | 43 | 0.28 | high |

    !!! note "Nota"
        Questi valori sono “ancoraggi” (population-level) e vengono poi modificati dai fattori successivi.

    ## Citogenetica (hazard ratio)

    Nel codice, alcune lesioni aumentano l’hazard (peggiorano) PFS/OS:

    | Lesione | HR PFS | HR OS | Categoria |
    | --- | ---: | ---: | --- |
    | `del(17p)` | 2.0 | 2.5 | high |
    | `t(4;14)` | 1.8 | 1.8 | high |
    | `t(14;16)` | 2.0 | 2.2 | high |
    | `t(14;20)` | 1.8 | 2.0 | high |
    | `1q21_gain` | 1.5 | 1.6 | intermediate |
    | `1p_deletion` | 1.4 | 1.5 | intermediate |
    | `double_hit` | 3.0 | 4.0 | very_high |

    Implementazione: per ogni lesione matchata:

    \[
    \mathrm{PFS}_{mult} \leftarrow \frac{\mathrm{PFS}_{mult}}{\mathrm{HR}_{pfs}},\quad
    \mathrm{OS}_{mult} \leftarrow \frac{\mathrm{OS}_{mult}}{\mathrm{HR}_{os}}
    \]

    (HR>1 → peggiora → moltiplicatore <1).

    ## Modificatori demografici/clinici

    ### Età

    `get_age_modifier(age)` restituisce (PFS_mult, OS_mult), ad esempio:

    | Età | PFS mult | OS mult |
    | ---: | ---: | ---: |
    | <65 | 1.00 | 1.00 |
    | 65–74 | 0.95 | 0.90 |
    | 75–79 | 0.85 | 0.75 |
    | ≥80 | 0.75 | 0.60 |

    ### ECOG

    `get_ecog_modifier(ecog)`:

    | ECOG | PFS mult | OS mult |
    | ---: | ---: | ---: |
    | 0 | 1.00 | 1.00 |
    | 1 | 0.95 | 0.95 |
    | 2 | 0.85 | 0.80 |
    | 3 | 0.70 | 0.55 |
    | 4 | 0.50 | 0.30 |

    ## Modificatori di risposta e MRD

    La risposta (`sCR/CR/VGPR/PR/SD/PD`) e MRD (`negative/positive/unknown`) modulano ulteriormente la stima.

    Esempio (estratto):

    - `CR`: PFS 1.2×, OS 1.3×
    - `PD`: PFS 0.4×, OS 0.5×
    - `MRD negative`: PFS 1.5×, OS 1.4×

    ## Linea di terapia

    Il parametro `line_of_therapy` introduce una penalità (relapse → prognosi tendenzialmente peggiore). Vedi `estimate_prognosis(...)`.

    ## Output “human readable”

    Nel modulo esistono helper che producono testi brevi IT/EN (timeline a 3/6/12/24 mesi). Questo è utile per UI.

    ## Figure

    ![Baseline survival by R-ISS (toy)](../assets/images/models/prognosis_riss_baseline.svg)

    !!! tip "Interpretazione"
        La figura mostra curve illustrate (non fit su dati reali) coerenti con le baseline per R-ISS.

    ## Esempio guidato (worked example)

    Parametri:

    | Input | Valore |
    | --- | --- |
    | R-ISS | III |
    | Cytogenetics | `["del(17p)", "1q21_gain"]` |
    | Age | 77 |
    | ECOG | 2 |
    | Response | `VGPR` |
    | MRD | `unknown` |
    | Line | 2 |

    Cosa succede:

    1. prendo baseline R-ISS III
    2. applico HR citogenetici (riducendo i multiplier)
    3. applico età + ECOG
    4. applico `VGPR`
    5. applico penalità line-of-therapy

    Output: `PrognosisEstimate` con mediane + probabilità ai timepoint + lista `modifiers_applied`.

=== "EN"
    This page explains the prognosis estimation logic used by the project.

    Primary source: `simulator/prognosis.py`.

    !!! warning "Disclaimer"
        PFS/OS estimates are **statistical** and **educational**. They do not replace clinical judgement.

    ## Quick find

    | I want… | Go to… |
    | --- | --- |
    | main estimator | `simulator/prognosis.py` → `estimate_prognosis(...)` |
    | cytogenetics / hazard ratios | `simulator/prognosis.py` → `CYTOGENETIC_HAZARD_RATIOS` |
    | relation to difficulty score | `Guides → Difficulty Scoring` |
    | staging (R-ISS) | `Guides → Database` (`Assessment.r_iss`) |

    ## Terminology

    - **PFS / OS**: progression-free survival / overall survival.
    - **HR (hazard ratio)**: HR>1 means worse prognosis.
    - **Baseline**: starting values before modifiers.

    ## Module goal

    The module produces a `PrognosisEstimate` with:

    - median PFS/OS (months)
    - PFS/OS probabilities at 12/24/36 months (and OS at 60 months)
    - risk category (`standard/intermediate/high/very_high`)
    - confidence (0–1) + a list of applied modifiers

    ## Logical pipeline (high-level)

    ```mermaid
    flowchart TD
      A[R-ISS baseline] --> B[Apply cytogenetic HR]
      B --> C[Apply age/ECOG modifiers]
      C --> D[Apply response + MRD modifiers]
      D --> E[Apply line-of-therapy penalty]
      E --> F[PrognosisEstimate]
    ```

    ## Baseline (R-ISS)

    Baseline is a per-R-ISS table (I/II/III) with medians and timepoint probabilities.

    | R-ISS | median PFS (mo) | median OS (mo) | OS 5y | risk |
    | --- | ---: | ---: | ---: | --- |
    | I | 66 | 120 | 0.75 | standard |
    | II | 42 | 83 | 0.55 | intermediate |
    | III | 29 | 43 | 0.28 | high |

    ## Cytogenetics (hazard ratios)

    Some lesions increase hazard (worse) PFS/OS:

    | Lesion | HR PFS | HR OS | Category |
    | --- | ---: | ---: | --- |
    | `del(17p)` | 2.0 | 2.5 | high |
    | `t(4;14)` | 1.8 | 1.8 | high |
    | `t(14;16)` | 2.0 | 2.2 | high |
    | `t(14;20)` | 1.8 | 2.0 | high |
    | `1q21_gain` | 1.5 | 1.6 | intermediate |
    | `1p_deletion` | 1.4 | 1.5 | intermediate |
    | `double_hit` | 3.0 | 4.0 | very_high |

    Implementation: for each matched lesion:

    \[
    \mathrm{PFS}_{mult} \leftarrow \frac{\mathrm{PFS}_{mult}}{\mathrm{HR}_{pfs}},\quad
    \mathrm{OS}_{mult} \leftarrow \frac{\mathrm{OS}_{mult}}{\mathrm{HR}_{os}}
    \]

    (HR>1 → worse → multiplier <1).

    ## Demographic/clinical modifiers

    ### Age

    `get_age_modifier(age)` returns (PFS_mult, OS_mult), e.g.:

    | Age | PFS mult | OS mult |
    | ---: | ---: | ---: |
    | <65 | 1.00 | 1.00 |
    | 65–74 | 0.95 | 0.90 |
    | 75–79 | 0.85 | 0.75 |
    | ≥80 | 0.75 | 0.60 |

    ### ECOG

    `get_ecog_modifier(ecog)`:

    | ECOG | PFS mult | OS mult |
    | ---: | ---: | ---: |
    | 0 | 1.00 | 1.00 |
    | 1 | 0.95 | 0.95 |
    | 2 | 0.85 | 0.80 |
    | 3 | 0.70 | 0.55 |
    | 4 | 0.50 | 0.30 |

    ## Response and MRD modifiers

    Response (`sCR/CR/VGPR/PR/SD/PD`) and MRD (`negative/positive/unknown`) further modulate the estimate.

    Example (excerpt):

    - `CR`: PFS 1.2×, OS 1.3×
    - `PD`: PFS 0.4×, OS 0.5×
    - `MRD negative`: PFS 1.5×, OS 1.4×

    ## Line of therapy

    `line_of_therapy` adds a penalty (relapse → typically worse prognosis). See `estimate_prognosis(...)`.

    ## Human-readable output

    The module includes helpers to produce short IT/EN timeline strings (3/6/12/24 months), useful for UI.

    ## Figure

    ![Baseline survival by R-ISS (toy)](../assets/images/models/prognosis_riss_baseline.svg)

    ## Worked example

    Inputs:

    | Input | Value |
    | --- | --- |
    | R-ISS | III |
    | Cytogenetics | `["del(17p)", "1q21_gain"]` |
    | Age | 77 |
    | ECOG | 2 |
    | Response | `VGPR` |
    | MRD | `unknown` |
    | Line | 2 |

    Steps:

    1. start from R-ISS III baseline
    2. apply cytogenetic HRs (multipliers decrease)
    3. apply age + ECOG
    4. apply `VGPR`
    5. apply line-of-therapy penalty

    Output: `PrognosisEstimate` with medians + timepoint probabilities + `modifiers_applied`.
