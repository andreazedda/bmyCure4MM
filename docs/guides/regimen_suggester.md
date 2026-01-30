# Regimen Suggester (rules-based recommendations)

=== "IT"
    Questa pagina documenta il motore di suggerimento regimi (educational) usato dal progetto.

    Fonte principale: `simulator/regimen_suggester.py`.

    !!! warning "Disclaimer"
        I suggerimenti sono **educational** e basati su regole: non sono raccomandazioni cliniche.

    ## Se cerchi…

    | Voglio… | Vai a… |
    | --- | --- |
    | entrypoint suggerimenti | `simulator/regimen_suggester.py` → `suggest_regimens(...)` |
    | database regimi | `FRONTLINE_REGIMENS`, `RELAPSED_REGIMENS` |
    | dettagli di un regime per nome | `get_regimen_details(...)` |
    | come sono memorizzati i regimi nel DB | `Guides → Database` (`clinic_regimen`) |

    ## Terminologia

    - **Frontline**: prima linea (line_of_therapy=1).
    - **Relapsed/Refractory**: linee successive (line_of_therapy>=2) e refrattarietà.
    - **Contraindication**: condizione che rende un regime da evitare.

    ## Input/Output (contratto)

    ### Input principali

    | Campo | Tipo | Significato |
    | --- | --- | --- |
    | `age` | int | età |
    | `transplant_eligible` | bool/None | eleggibilità ASCT (se None viene inferita) |
    | `ecog` | int | ECOG 0–4 |
    | `r_iss` | str | I/II/III |
    | `high_risk_cytogenetics` | bool | flag high-risk |
    | `line_of_therapy` | int | 1=frontline, 2+=relapsed |
    | `prior_therapies` | list[str] | agenti già ricevuti |
    | `refractory_to` | list[str] | agenti refrattari |
    | `renal_function` | str | `normal/impaired/dialysis` |
    | `neuropathy_grade` | int | 0–4 |
    | `cardiac_issues` | bool | comorbidità cardiache |

    ### Output

    `suggest_regimens(...)` restituisce un dict con liste:

    - `preferred`
    - `alternative`
    - `consider_in_clinical_trial`
    - `avoid`
    - `patient_factors` (spiegazioni “auditabili”)
    - `disclaimer`

    Ogni item è un `RegimenSuggestion` serializzato (`to_dict()`), più un campo `rationale` aggiunto dal motore.

    ## Algoritmo (overview)

    Il motore è **rule-based**: applica regole in cascata su line-of-therapy e fattori clinici.

    ```mermaid
    flowchart TD
      A[Normalize inputs] --> B{Infer transplant eligible?}
      B --> C[Record patient_factors]
      C --> D{line_of_therapy == 1?}
      D -- yes --> E[Frontline rules]
      D -- no --> F[Relapsed rules]
      E --> G[Post-filter (contraindications)]
      F --> G
      G --> H[Return preferred/alt/avoid]
    ```

    ## Regole chiave (esempi dal codice)

    ### 1) Inferenza trapianto (se None)

    - se `age < 65` e `ecog <= 1` → `transplant_eligible=True`
    - se `age >= 75` → `False`
    - altrimenti → `None` e aggiunge nota in `patient_factors`

    ### 2) Neuropatia

    Se `neuropathy_grade >= 2`:
    - aggiunge fattore “avoid/reduce bortezomib”
    - influenza le scelte quando sono disponibili alternative (es. `KRd` vs `VRd`)

    ### 3) High-risk cytogenetics

    Se `high_risk_cytogenetics=True`:
    - suggerisce di considerare quadruplet (es. `Dara-VRd`) in frontline

    ### 4) Funzione renale

    Se `renal_function == "dialysis"`:
    - annota “adjust lenalidomide dosing, bortezomib preferred”

    ## Esempi guidati

    ### Esempio A — Frontline, trapianto eleggibile, standard risk

    ```python
    from simulator.regimen_suggester import suggest_regimens

    out = suggest_regimens(age=60, ecog=0, line_of_therapy=1, high_risk_cytogenetics=False)
    ```

    Atteso:
    - `preferred`: include tipicamente `VRd` (o `Dara-VRd` se high-risk)

    ### Esempio B — Frontline, high-risk + neuropatia

    ```python
    out = suggest_regimens(
        age=62,
        ecog=1,
        line_of_therapy=1,
        high_risk_cytogenetics=True,
        neuropathy_grade=2,
    )
    ```

    Atteso:
    - `patient_factors` include high-risk + neuropatia
    - `preferred/alternative` tende a favorire regimi con meno bortezomib (quando disponibili)

    ## Mapping al codice

    - regimi e metadata: `RegimenSuggestion`, `FRONTLINE_REGIMENS`, `RELAPSED_REGIMENS`
    - motore: `suggest_regimens(...)`
    - lookup: `get_regimen_details(...)`

=== "EN"
    This page documents the project’s regimen suggestion engine (educational, rules-based).

    Primary source: `simulator/regimen_suggester.py`.

    !!! warning "Disclaimer"
        Suggestions are **educational** and rules-based. They are not clinical recommendations.

    ## Quick find

    | I want… | Go to… |
    | --- | --- |
    | entrypoint | `simulator/regimen_suggester.py` → `suggest_regimens(...)` |
    | regimen “database” | `FRONTLINE_REGIMENS`, `RELAPSED_REGIMENS` |
    | fetch regimen details by name | `get_regimen_details(...)` |
    | how regimens are stored in DB | `Guides → Database` (`clinic_regimen`) |

    ## Terminology

    - **Frontline**: first line (line_of_therapy=1).
    - **Relapsed/Refractory**: later lines (line_of_therapy>=2) and refractoriness.
    - **Contraindication**: condition that makes a regimen unsuitable.

    ## Input/Output (contract)

    ### Key inputs

    | Field | Type | Meaning |
    | --- | --- | --- |
    | `age` | int | age |
    | `transplant_eligible` | bool/None | ASCT eligibility (if None it is inferred) |
    | `ecog` | int | ECOG 0–4 |
    | `r_iss` | str | I/II/III |
    | `high_risk_cytogenetics` | bool | high-risk flag |
    | `line_of_therapy` | int | 1=frontline, 2+=relapsed |
    | `prior_therapies` | list[str] | prior agents |
    | `refractory_to` | list[str] | refractory agents |
    | `renal_function` | str | `normal/impaired/dialysis` |
    | `neuropathy_grade` | int | 0–4 |
    | `cardiac_issues` | bool | cardiac comorbidities |

    ### Output

    `suggest_regimens(...)` returns a dict with lists:

    - `preferred`
    - `alternative`
    - `consider_in_clinical_trial`
    - `avoid`
    - `patient_factors` (auditable explanations)
    - `disclaimer`

    Each item is a serialized `RegimenSuggestion` (`to_dict()`), plus an extra `rationale` field added by the engine.

    ## Algorithm (overview)

    The engine is **rule-based**: it applies cascading rules based on line-of-therapy and patient factors.

    ```mermaid
    flowchart TD
      A[Normalize inputs] --> B{Infer transplant eligible?}
      B --> C[Record patient_factors]
      C --> D{line_of_therapy == 1?}
      D -- yes --> E[Frontline rules]
      D -- no --> F[Relapsed rules]
      E --> G[Post-filter (contraindications)]
      F --> G
      G --> H[Return preferred/alt/avoid]
    ```

    ## Key rules (examples from code)

    ### 1) Transplant inference (if None)

    - if `age < 65` and `ecog <= 1` → `transplant_eligible=True`
    - if `age >= 75` → `False`
    - else → `None` and a note is added to `patient_factors`

    ### 2) Neuropathy

    If `neuropathy_grade >= 2`:
    - adds “avoid/reduce bortezomib”
    - influences selection where alternatives exist (e.g., `KRd` vs `VRd`)

    ### 3) High-risk cytogenetics

    If `high_risk_cytogenetics=True`:
    - suggests considering quadruplet induction (e.g., `Dara-VRd`) in frontline

    ### 4) Renal function

    If `renal_function == "dialysis"`:
    - notes “adjust lenalidomide dosing, bortezomib preferred”

    ## Worked examples

    ### Example A — Frontline, transplant eligible, standard risk

    ```python
    from simulator.regimen_suggester import suggest_regimens

    out = suggest_regimens(age=60, ecog=0, line_of_therapy=1, high_risk_cytogenetics=False)
    ```

    Expected:
    - `preferred`: typically includes `VRd` (or `Dara-VRd` if high-risk)

    ### Example B — Frontline, high-risk + neuropathy

    ```python
    out = suggest_regimens(
        age=62,
        ecog=1,
        line_of_therapy=1,
        high_risk_cytogenetics=True,
        neuropathy_grade=2,
    )
    ```

    Expected:
    - `patient_factors` includes high-risk + neuropathy
    - `preferred/alternative` tends to favor less bortezomib when alternatives exist

    ## Code mapping

    - regimen definitions: `RegimenSuggestion`, `FRONTLINE_REGIMENS`, `RELAPSED_REGIMENS`
    - engine: `suggest_regimens(...)`
    - lookup: `get_regimen_details(...)`
