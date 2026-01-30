# Optimization Theory (Optimization Lab)

=== "IT"
    Questa sezione documenta **come** l’Optimization Lab costruisce e valuta le soluzioni (multi-obiettivo / Pareto).

    !!! tip "Come usare questa pagina"
        - Se ti serve “capire al volo”: leggi **Obiettivi**, **Vincoli**, **Pareto** e guarda la figura “Pareto front”.
        - Se devi modificarlo nel codice: leggi **Mapping al codice** e **Algorithm (Optuna)**.
        - Se devi spiegarlo a un clinico/PM: leggi **Interpretazione pratica** e **Esempio guidato**.

    ## Se cerchi…

    | Voglio… | Vai a… |
    | --- | --- |
    | come funziona Optuna nel repo | `simulator/optim.py` |
    | mapping obiettivi/vincoli | `simulator/objectives.py` |
    | spazio di ricerca | `simulator/search_space.py` |
    | esempi visivi | `Guides → Simulations Gallery` |

    ## Contesto: perché multi-obiettivo

    In MM (e più in generale in terapia):

    - aumentare la **dose** può aumentare efficacia ma anche tossicità
    - ridurre esposizione può migliorare tollerabilità ma peggiorare risposta

    Quindi non esiste un “best” unico: esiste un insieme di compromessi validi.

    ## Formalizzazione (problema di ottimizzazione)

    Il sistema risolve (in forma concettuale) un problema del tipo:

    \[
    \max_{x \in \mathcal{X}} \; \mathbf{f}(x) = \left(f_1(x), f_2(x), f_3(x)\right)
    \quad \text{soggetto a} \quad g(x) \le 0
    \]

    dove:

    - \(x\) sono le **decision variables** (dosi, orizzonte, parametri di interazione, ecc.)
    - \(\mathcal{X}\) è lo **spazio di ricerca** (range + step + categorie)
    - \(\mathbf{f}(x)\) sono gli **obiettivi** (efficacia, sicurezza, esposizione)
    - \(g(x)\) è un **vincolo hard** (feasibility) costruito dai KPI (es. `healthy_loss`)

    ## Terminologia

    - **Decision variables**: parametri che l’ottimizzazione può scegliere (dosi, orizzonte, ecc.).
    - **Objective function**: funzione che mappa parametri → punteggio (qui 3 punteggi).
    - **Constraint**: condizione “hard” che rende una soluzione non valida (feasibility).
    - **Feasible**: soluzione che soddisfa tutti i vincoli.
    - **Dominance**: relazione “meglio in tutto e almeno meglio in uno”.
    - **Pareto front**: insieme delle soluzioni non dominate.

    ## Obiettivi e vincoli

    Nel codice (`simulator/objectives.py`) gli obiettivi sono:

    | Nome | Definizione | Direzione |
    | --- | --- | --- |
    | Efficacy | `tumor_reduction` | massimizza |
    | Safety | `-healthy_loss` | massimizza (equivale a minimizzare loss) |
    | Exposure | `-sum(AUC)` | massimizza (equivale a minimizzare esposizione) |

    Vincolo “hard”:

    - `healthy_loss <= 0.25` (soluzione “feasible”)

    ### Forma matematica (intuizione)

    Indicando i parametri ottimizzati con \(x\) (dosi, orizzonte, ecc.), e gli output del modello con \(y(x)\):

    - \(f_1(x) = \mathrm{tumor\_reduction}(y(x))\)  (massimizza)
    - \(f_2(x) = -\mathrm{healthy\_loss}(y(x))\)    (massimizza)
    - \(f_3(x) = -\sum_i \mathrm{AUC}_i(y(x))\)     (massimizza)

    Vincolo:

    - \(g(x) = \mathrm{healthy\_loss}(y(x)) \le 0.25\)

    ```mermaid
    flowchart TD
      S[Trial params] --> R[Run model]
      R --> K[KPI summary]
      K --> C{Constraints OK?}
      C -- no --> P[Penalty -1e6]
      C -- yes --> O["Objectives (3D)"]
      O --> PF[Pareto front]
    ```

    ### Dominanza (definizione formale)

    Con **tutti gli obiettivi massimizzati**, una soluzione \(A\) **domina** \(B\) se:

    \[
    \forall k,\; f_k(A) \ge f_k(B)\quad \text{e} \quad \exists k,\; f_k(A) > f_k(B)
    \]

    La **Pareto front** è l’insieme delle soluzioni non dominate (tipicamente filtrate anche per feasibility).

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

    ### Cosa sta facendo “davvero” Optuna (spiegazione semplice)

    - Ogni **trial** propone un set di parametri \(x\) dentro i range del search space.
    - Il sistema esegue la simulazione e produce KPI.
    - I KPI vengono trasformati in 3 obiettivi \((f_1, f_2, f_3)\).
    - Optuna usa TPE per proporre parametri “promettenti” basandosi sui risultati passati.

    ```mermaid
    sequenceDiagram
      autonumber
      participant O as Optuna sampler
      participant S as Search space
      participant M as Simulation model
      participant K as KPIs
      O->>S: sample parameters x
      S-->>O: x
      O->>M: run(x)
      M-->>K: summary (KPIs)
      K-->>O: objectives (f1,f2,f3) + constraint
    ```

    ### Pruning (perché esiste)

    Il **pruner** (es. `MedianPruner`) interrompe trial “poco promettenti” quando esiste sufficiente evidenza che non miglioreranno rispetto alla mediana dei trial precedenti.

    Beneficio:
    - riduce costo computazionale quando la simulazione è lenta (ODE + cohort)

    Rischio:
    - se i primi trial sono rumorosi o poco rappresentativi, può eliminare soluzioni buone (soprattutto con obiettivi multi-dimensionali)

    ### Penalty (perché -1e6)

    Nel codice, se una soluzione viola i vincoli, viene penalizzata sottraendo un valore enorme.

    **Perché funziona**: in uno spazio di ricerca finito, qualunque soluzione “feasible” avrà punteggi molto migliori di una penalizzata, quindi la Pareto front finale tende a contenere solo feasible.

    **Quando è un problema**: se quasi tutte le soluzioni sono infeasible, il modello non impara; in quel caso conviene:

    - restringere search space
    - rendere i vincoli più realistici
        - aggiungere prior/initial trials “safe”

    ### Alternative teoriche (non ancora implementate qui)

    Se in futuro vuoi rendere la gestione vincoli più “pulita”, le alternative classiche sono:

    - **Feasibility-first sorting**: prima feasible, poi Pareto
    - **Epsilon-constraint**: massimizzi un obiettivo e imponi soglie sugli altri (es. `healthy_loss <= ε`)
    - **Scalarization**: combini obiettivi con pesi \(w_k\) (attenzione: perdi parte della front)

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

    ### Figura: Pareto front (toy)

    ![Pareto front](../assets/images/models/pareto_front.svg)

    ## Convergenza e qualità della front (come “capire se sta funzionando”)

    In multi-obiettivo, metriche tipiche (teoriche) sono:

    - **Hypervolume (HV)**: volume dominato dalla front rispetto a un reference point (più è alto, meglio è)
    - **Spread / diversity**: quanto la front copre bene lo spazio dei compromessi

    Nel repo, oggi la validazione è principalmente:
    - confrontare visivamente la front (tabella/plot)
    - controllare che le soluzioni siano **feasible** e “clinicamente sensate”

    ## Interpretazione pratica

    - a parità di `tumor_reduction`, preferisci `healthy_loss` più basso
    - evita miglioramenti marginali di efficacia con grande costo in tossicità
    - l’esposizione (AUC) aiuta a scegliere regimi “più parsimoniosi”

    ## Esempio guidato (worked example)

    Supponiamo di confrontare tre soluzioni (output KPI già calcolati):

    | Soluzione | tumor_reduction | healthy_loss | sum(AUC) | Feasible? |
    | --- | ---:| ---:| ---:| --- |
    | A | 0.80 | 0.18 | 120 | ✅ |
    | B | 0.85 | 0.29 | 140 | ❌ (vincolo) |
    | C | 0.78 | 0.12 | 110 | ✅ |

    Interpretazione:

    - B ha efficacia più alta ma viola il vincolo → viene penalizzata e tipicamente esclusa dalla Pareto front
    - A vs C: A è più efficace, C è più sicura e meno esposta → entrambe possono essere Pareto

    ## Mapping al codice (punti chiave)

    - Obiettivi: `simulator/objectives.py`
    - Search space: `simulator/search_space.py`
    - Loop Optuna + penalty + Pareto: `simulator/optim.py`

=== "EN"
    This section explains how the Optimization Lab evaluates solutions (multi-objective / Pareto).

    !!! tip "How to read this page"
        - Quick understanding: read **Objectives**, **Constraints**, **Pareto**, and check the “Pareto front” figure.
        - If you need to change the code: read **Code mapping** and **Algorithm (Optuna)**.
        - If you must explain it to a clinician/PM: read **Practical interpretation** and **Worked example**.

    ## Quick find

    | I want… | Go to… |
    | --- | --- |
    | Optuna implementation | `simulator/optim.py` |
    | objectives/constraints mapping | `simulator/objectives.py` |
    | search space | `simulator/search_space.py` |
    | visual examples | `Guides → Simulations Gallery` |

    ## Why multi-objective

    In practice:

    - higher dose can increase efficacy but also toxicity
    - lower exposure can improve tolerability but reduce response

    So there is no single best solution—there is a set of valid trade-offs.

    ## Formal statement (optimization problem)

    Conceptually, the system solves:

    \[
    \max_{x \in \mathcal{X}} \; \mathbf{f}(x) = \left(f_1(x), f_2(x), f_3(x)\right)
    \quad \text{subject to} \quad g(x) \le 0
    \]

    where:
    - \(x\) are the **decision variables** (doses, horizon, interaction knobs, etc.)
    - \(\mathcal{X}\) is the **search space** (ranges/steps/categories)
    - \(\mathbf{f}(x)\) are the **objectives** (efficacy, safety, exposure)
    - \(g(x)\) is a **hard constraint** derived from KPIs (e.g., `healthy_loss`)

    ## Terminology

    - **Decision variables**: parameters the optimizer can choose (doses, horizon, etc.).
    - **Objective function**: maps parameters → a score (here: 3 scores).
    - **Constraint**: hard validity condition (feasibility).
    - **Feasible**: satisfies all constraints.
    - **Dominance**: “not worse in all and better in at least one”.
    - **Pareto front**: set of non-dominated feasible solutions.

    ## Objectives & constraints

    | Name | Definition | Direction |
    | --- | --- | --- |
    | Efficacy | `tumor_reduction` | maximize |
    | Safety | `-healthy_loss` | maximize (i.e., minimize loss) |
    | Exposure | `-sum(AUC)` | maximize (i.e., minimize exposure) |

    Hard constraint:

    - `healthy_loss <= 0.25` (feasible solution)

    ### Math form (intuition)

    Let \(x\) be the decision variables and \(y(x)\) the model outputs:

    - \(f_1(x) = \mathrm{tumor\_reduction}(y(x))\)
    - \(f_2(x) = -\mathrm{healthy\_loss}(y(x))\)
    - \(f_3(x) = -\sum_i \mathrm{AUC}_i(y(x))\)

    Constraint:

    - \(g(x) = \mathrm{healthy\_loss}(y(x)) \le 0.25\)

    ```mermaid
    flowchart TD
      S[Trial params] --> R[Run model]
      R --> K[KPI summary]
      K --> C{Constraints OK?}
      C -- no --> P[Penalty -1e6]
      C -- yes --> O["Objectives (3D)"]
      O --> PF[Pareto front]
    ```

    ### Dominance (formal)

    With all objectives **maximized**, a solution \(A\) dominates \(B\) if:

    \[
    \forall k,\; f_k(A) \ge f_k(B)\quad \text{and} \quad \exists k,\; f_k(A) > f_k(B)
    \]

    The **Pareto front** is the set of non-dominated (typically feasible) solutions.

    ## Search space

    | Parameter | Type | Range/step | Notes |
    | --- | --- | --- | --- |
    | `lenalidomide_dose` | float | 0..40 step 2.5 | dose |
    | `bortezomib_dose` | float | 0..1.6 step 0.1 | dose |
    | `daratumumab_dose` | float | 0..16 step 1.0 | dose |
    | `time_horizon` | int | 56..224 step 28 | 2–8 cycles |
    | `interaction_strength` | float | 0..0.15 step 0.01 | interaction |

    ## Algorithm (Optuna)

    - sampler: `TPESampler(seed=..., multivariate=True)`
    - pruner: `MedianPruner(n_startup_trials=10, ...)`
    - directions: maximize/maximize/maximize

    ### What Optuna is doing (plain language)

    - Each **trial** proposes a parameter vector \(x\) inside the search space.
    - The simulator runs and produces KPIs.
    - KPIs are mapped to the 3 objectives \((f_1,f_2,f_3)\).
    - TPE proposes “promising” samples based on past trials.

    ```mermaid
    sequenceDiagram
      autonumber
      participant O as Optuna sampler
      participant S as Search space
      participant M as Simulation model
      participant K as KPIs
      O->>S: sample parameters x
      S-->>O: x
      O->>M: run(x)
      M-->>K: summary (KPIs)
      K-->>O: objectives (f1,f2,f3) + constraint
    ```

    ### Pruning (why it exists)

    The pruner can stop trials that are unlikely to improve, saving compute when simulations are expensive.

    Trade-off:
    - too aggressive pruning can remove good-but-late solutions (especially under noise)

    ### Penalty (-1e6)

    If a solution violates hard constraints, the implementation applies a large negative penalty.

    Why it works:
    - any feasible solution will dominate a penalized one, so the final front tends to be feasible

    When it fails:
    - if most trials are infeasible, the optimizer cannot learn → narrow the search space or add safe priors

    ### Theoretical alternatives (not implemented here)

    - **Feasibility-first sorting** (feasible solutions ranked first)
    - **Epsilon-constraint** (optimize one objective while constraining the others)
    - **Scalarization** \(F(x)=\sum_k w_k f_k(x)\) (simpler but loses parts of the front)

    ## Pareto front

    A solution \(A\) dominates \(B\) if it is not worse on all objectives and better on at least one.

    ![Pareto front](../assets/images/models/pareto_front.svg)

    ## Convergence & front quality (how to tell it works)

    Common multi-objective metrics (theory):
    - **Hypervolume (HV)**: dominated volume w.r.t. a reference point (higher is better)
    - **Diversity/spread**: coverage of the trade-off space

    In this repo, the practical checks are:
    - solutions are feasible and clinically sensible
    - the Pareto set contains meaningful trade-offs (not all clustered)

    ## Practical interpretation

    - for equal `tumor_reduction`, prefer lower `healthy_loss`
    - avoid marginal efficacy gains with large toxicity costs
    - exposure (AUC) helps select “parsimonious” regimens

    ## Worked example

    | Solution | tumor_reduction | healthy_loss | sum(AUC) | Feasible? |
    | --- | ---:| ---:| ---:| --- |
    | A | 0.80 | 0.18 | 120 | ✅ |
    | B | 0.85 | 0.29 | 140 | ❌ |
    | C | 0.78 | 0.12 | 110 | ✅ |

    - B violates the constraint → gets penalized
    - A vs C: both can be Pareto since each wins on different dimensions

    ## Code mapping

    - objectives/constraints: `simulator/objectives.py`
    - search space: `simulator/search_space.py`
    - Optuna loop + penalty + Pareto: `simulator/optim.py`
