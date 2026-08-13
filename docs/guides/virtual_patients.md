# Virtual Patients (archetypes + parameter distributions)

> **E1 terminology notice:** this historical page uses “virtual patient” for a
> sampled synthetic model state. The distributions are `HEURISTIC`, are not a
> validated representation of a patient population, and cannot support
> patient-specific prediction.

=== "IT"
    Questa pagina descrive come il progetto genera **virtual patients** per simulazioni in-silico (cohort) in modo riproducibile e trasparente.

    Fonte principale: `simulator/virtual_patients.py`.

    ## Se cerchi…

    | Voglio… | Vai a… |
    | --- | --- |
    | elenco archetipi e distribuzioni | `simulator/virtual_patients.py` → `get_archetype_library()` |
    | come viene campionato un parametro | `ParameterDistribution.sample(...)` |
    | incertezza/cohort e bande | `Guides → Mathematical Models` + `Guides → Simulations Gallery` |
    | stato epistemico | `Governance → Epistemic Labels` |

    ## Terminologia

    - **Archetipo**: popolazione con distribuzioni parametriche definite.
    - **Cohort**: insieme di pazienti virtuali campionati (N).
    - **Seed**: seme random per riproducibilità.

    ## Concetto: perché archetipi

    Un **archetipo** rappresenta uno stato sintetico con distribuzioni euristiche definite per:

    - età/sesso
    - R-ISS
    - citogenetica
    - biologia tumorale (burden, growth rate, carrying capacity)
    - fitness (ECOG, comorbidità)
    - lab/parametri di base

    ## Pipeline (at a glance)

    ```mermaid
    flowchart TD
      A[Choose archetype] --> B[Sample parameters]
      B --> C[Create cohort (N patients)]
      C --> D[Run simulator]
      D --> E[Uncertainty bands + summary KPIs]
    ```

    ## Tipi di distribuzione supportati

    Il dataclass `ParameterDistribution` supporta:

    | Tipo | Quando usarla | Note |
    | --- | --- | --- |
    | `lognormal` | positivi continui (cell counts, growth rates) | robusta per scale log |
    | `normal` | continui “simmetrici” (lab values) | con clamp bounds |
    | `beta` | parametri (0..1) (fractions) | convertendo mean/std → (α,β) |
    | `uniform` | prior piatto | bounds richiesti |
    | `categorical` | scelte discrete | mappa {valore: probabilità} |

    ### Lognormal (conversione mean/std)

    Nel codice, `mean` e `std` sono riferiti alla **variabile reale** \(X\), mentre la lognormal è costruita da \(Z \sim \mathcal{N}(\mu,\sigma^2)\) con:

    \[
    X = e^Z
    \]

    e conversione usata:

    \[
    \mu = \ln\left(\frac{m^2}{\sqrt{m^2+s^2}}\right),\quad
    \sigma = \sqrt{\ln\left(1+\frac{s^2}{m^2}\right)}
    \]

    ### Beta (conversione mean/std)

    Per `beta`, il codice ricava \(\alpha,\beta\) dalla media e varianza:

    \[
    \mathrm{mean}=\frac{\alpha}{\alpha+\beta},\quad
    \mathrm{var}=\frac{\alpha\beta}{(\alpha+\beta)^2(\alpha+\beta+1)}
    \]

    Poi scala eventualmente a \([a,b]\).

    ## Libreria archetipi (overview)

    Gli archetipi sono definiti come enum `PatientArchetype`, con un `ArchetypeDefinition` che include prevalenza e distribuzioni.

    Esempi (nomi dal codice):

    - `NEWLY_DIAGNOSED_STANDARD_RISK`
    - `NEWLY_DIAGNOSED_HIGH_RISK`
    - `FRAIL_ELDERLY`
    - `RELAPSED_REFRACTORY`
    - `SMOLDERING_MYELOMA`
    - `TRANSPLANT_ELIGIBLE`
    - `TRANSPLANT_INELIGIBLE`

    !!! tip "Come usare in pratica"
        Quando vuoi vedere variabilità (incertezza) non usare `cohort_size=1`: aumenta `cohort_size` e fissa un `seed` per riproducibilità.

    ## Reproducibilità (seed)

    `ParameterDistribution.sample(..., seed=...)` imposta `np.random.seed(seed)`. Questo rende riproducibile il campionamento **per chiamata**.

    Se vuoi riproducibilità “end-to-end” per una cohort:
    - usa un seed globale (o un `np.random.Generator`)
    - evita di resettare il seed in più punti con lo stesso valore

    ## Esempio guidato (mini)

    Obiettivo: creare 50 pazienti “frail elderly” e campionare una variabile (es. età).

    ```python
    from simulator.virtual_patients import get_archetype_library, PatientArchetype

    lib = get_archetype_library()
    frail = lib[PatientArchetype.FRAIL_ELDERLY]

    ages = frail.age_distribution.sample(n=50, seed=123)
    ```

=== "EN"
    This page explains how the project generates **virtual patients** for in-silico simulation cohorts in a reproducible and transparent way.

    Primary source: `simulator/virtual_patients.py`.

    ## Quick find

    | I want… | Go to… |
    | --- | --- |
    | archetypes and distributions | `simulator/virtual_patients.py` → `get_archetype_library()` |
    | how a parameter is sampled | `ParameterDistribution.sample(...)` |
    | cohort uncertainty bands | `Guides → Mathematical Models` + `Guides → Simulations Gallery` |
    | epistemic status | `Governance → Epistemic Labels` |

    ## Terminology

    - **Archetype**: population with specified parameter distributions.
    - **Cohort**: sampled set of virtual patients (N).
    - **Seed**: random seed for reproducibility.

    ## Why archetypes

    An **archetype** represents a synthetic state with heuristic distributions for:

    - age/sex
    - R-ISS
    - cytogenetics
    - tumor biology (burden, growth rate, carrying capacity)
    - fitness (ECOG, comorbidities)
    - key lab/baseline parameters

    ## Pipeline (at a glance)

    ```mermaid
    flowchart TD
      A[Choose archetype] --> B[Sample parameters]
      B --> C[Create cohort (N patients)]
      C --> D[Run simulator]
      D --> E[Uncertainty bands + summary KPIs]
    ```

    ## Supported distribution types

    `ParameterDistribution` supports:

    | Type | When to use | Notes |
    | --- | --- | --- |
    | `lognormal` | positive continuous (cell counts, growth rates) | stable on log scales |
    | `normal` | “symmetric” continuous (lab values) | with bounds clamp |
    | `beta` | (0..1) parameters (fractions) | mean/std → (α,β) |
    | `uniform` | flat prior | bounds required |
    | `categorical` | discrete choices | mapping {value: probability} |

    ### Lognormal (mean/std conversion)

    In code, `mean` and `std` refer to the **real variable** \(X\), while lognormal is built from \(Z \sim \mathcal{N}(\mu,\sigma^2)\) with:

    \[
    X = e^Z
    \]

    and the conversion:

    \[
    \mu = \ln\left(\frac{m^2}{\sqrt{m^2+s^2}}\right),\quad
    \sigma = \sqrt{\ln\left(1+\frac{s^2}{m^2}\right)}
    \]

    ### Beta (mean/std conversion)

    For `beta`, the code derives \(\alpha,\beta\) from mean/variance:

    \[
    \mathrm{mean}=\frac{\alpha}{\alpha+\beta},\quad
    \mathrm{var}=\frac{\alpha\beta}{(\alpha+\beta)^2(\alpha+\beta+1)}
    \]

    Then it optionally rescales to \([a,b]\).

    ## Archetype library (overview)

    Archetypes are defined as enum `PatientArchetype`, with an `ArchetypeDefinition` including prevalence and distributions.

    Examples (names from code):

    - `NEWLY_DIAGNOSED_STANDARD_RISK`
    - `NEWLY_DIAGNOSED_HIGH_RISK`
    - `FRAIL_ELDERLY`
    - `RELAPSED_REFRACTORY`
    - `SMOLDERING_MYELOMA`
    - `TRANSPLANT_ELIGIBLE`
    - `TRANSPLANT_INELIGIBLE`

    !!! tip "Practical usage"
        If you want variability/uncertainty, avoid `cohort_size=1`: increase `cohort_size` and fix a `seed` for reproducibility.

    ## Reproducibility (seed)

    `ParameterDistribution.sample(..., seed=...)` uses `np.random.seed(seed)`, making sampling reproducible **per call**.

    For end-to-end reproducible cohorts:
    - use a single global seed (or a `np.random.Generator`)
    - avoid resetting the seed in multiple places with the same value

    ## Worked mini example

    Goal: create 50 “frail elderly” patients and sample one variable (e.g., age).

    ```python
    from simulator.virtual_patients import get_archetype_library, PatientArchetype

    lib = get_archetype_library()
    frail = lib[PatientArchetype.FRAIL_ELDERLY]

    ages = frail.age_distribution.sample(n=50, seed=123)
    ```
