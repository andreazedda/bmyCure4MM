# Simulator KPIs (how to read)

=== "IT"

    Questa guida spiega **come leggere** l’output del Simulatore (grafici + `results_summary`) e come confrontare due run.

    !!! warning "Nota importante"
        I KPI sono **output del modello** (proxy / indicatori interni), non endpoint clinici validati.

    ## Dove li trovi

    Dopo una simulazione (attempt) trovi tipicamente:

    - un **plot** (HTML Plotly) con traiettorie e curve
    - un **CSV** con le time-series
    - un **summary** con KPI aggregati (`results_summary`)

    Per la reference tecnica su dove sono salvati e con quali chiavi, vedi:

    - `Reference → Simulator Parameters` (sezione “Output salvati”)

    ## Screenshot (esempio)

    ![Esempio output simulazione (plot)](../assets/images/screenshots/simulation_plot.png)

    ## KPI principali (cosa significano)

    I nomi possono variare leggermente con l’evoluzione del modello, ma i concetti chiave sono:

    - **Tumor reduction** (`tumor_reduction`): riduzione relativa del tumore rispetto al baseline; in generale “più alto” = risposta migliore (nel contesto del modello).
    - **Healthy loss** (`healthy_loss`): perdita relativa della componente “healthy”; in generale “più alto” = più tossicità/immunosoppressione (proxy).
    - **Nadir** (`nadir`): minimo raggiunto (es. tumore) e quando avviene; utile per capire profondità e timing della risposta.
    - **Milestones** (`milestones`): tempi a cui vengono raggiunte soglie di risposta/recupero (se previste dal modello).
    - **Time to recurrence** (`time_to_recurrence`): proxy di “durabilità”: quanto tempo passa prima che la traiettoria indichi una ripresa.
    - **Durability index** (`durability_index`): indice sintetico di stabilità della risposta (proxy, dipende dall’implementazione corrente).

    ### AUC (Area Under Curve)

    **AUC** in questo contesto è tipicamente usata come proxy di **esposizione** o “carico” sotto una curva.

    - Se applicata a curve di **concentrazione**, AUC ~ esposizione cumulativa.
    - Se applicata a curve di **stato** (tumore/healthy), è un aggregato matematico utile per confronti tra run, ma va interpretato con cautela.

    !!! note "Custom drug"
        Quando abiliti `Custom drug (experimental)`, compare una chiave extra tipo `auc[custom_*]` (esposizione del nuovo agente).

    ## Come confrontare due simulazioni (pratico)

    1) Confronta **forme e tempi** delle curve (non solo un numero): risposta precoce vs tardiva.

    2) Guarda il tradeoff:

    - aumento `tumor_reduction` con aumento `healthy_loss` = possibile “efficacia” maggiore ma più tossicità (proxy)

    3) Se usi **coorte** (`cohort_size>1`), interpreta i KPI con le bande:

    - **p05/p95** danno un’idea della variabilità del modello

    Per regole decisionali/threshold (quando presenti), vedi:

    - `Guides → Decision Algorithm`

=== "EN"

    This page explains **how to read** Simulator outputs (plots + `results_summary`) and how to compare two runs.

    !!! warning "Important note"
        KPIs are **model outputs** (internal proxies), not clinically validated endpoints.

    ## Where to find them

    Each simulation attempt typically stores:

    - a Plotly **plot** (HTML)
    - a trajectory **CSV**
    - a **summary** with aggregated KPIs (`results_summary`)

    For the technical reference (keys and persistence), see:

    - `Reference → Simulator Parameters` ("Stored outputs" section)

    ## Screenshot (example)

    ![Simulation output example (plot)](../assets/images/screenshots/simulation_plot.png)

    ## Main KPIs (meaning)

    Names may evolve over time, but the core ideas are:

    - **Tumor reduction** (`tumor_reduction`): relative reduction vs baseline; generally “higher” = better response (within the model).
    - **Healthy loss** (`healthy_loss`): relative loss of the “healthy” compartment; generally “higher” = more toxicity/immunosuppression (proxy).
    - **Nadir** (`nadir`): minimum reached (e.g. tumor) and timing.
    - **Milestones** (`milestones`): times when response thresholds are met (if defined).
    - **Time to recurrence** (`time_to_recurrence`): durability proxy.
    - **Durability index** (`durability_index`): synthetic stability index (proxy; depends on current implementation).

    ### AUC (Area Under Curve)

    Here AUC is typically an **exposure** / cumulative-load proxy.

    - For **concentration** curves, AUC ~ cumulative exposure.
    - For **state** curves (tumor/healthy), it’s a useful aggregate for comparisons but should be interpreted cautiously.

    !!! note "Custom drug"
        With `Custom drug (experimental)`, you’ll see an extra key like `auc[custom_*]` (exposure of the injected agent).

    ## Comparing two attempts (practical)

    1) Compare curve **shape and timing**, not just a single number.

    2) Look at the tradeoff:

    - higher `tumor_reduction` together with higher `healthy_loss` = “more effect” but more toxicity (proxy)

    3) With **cohorts** (`cohort_size>1`), use bands/stats:

    - **p05/p95** provide a sense of variability

    For decision rules/thresholds (when present), see:

    - `Guides → Decision Algorithm`
