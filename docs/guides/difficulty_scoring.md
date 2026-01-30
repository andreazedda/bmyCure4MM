# Difficulty Scoring (scenario difficulty + expected outcomes)

=== "IT"
    Questa pagina spiega il **modello matematico** che produce:

    - `difficulty_score` (0–100) e `difficulty_level` (Very Easy → Very Hard)
    - stime “toy” di **response**, **toxicity risk** e **survival** (educational)

    Fonte principale: `simulator/difficulty_scoring.py`.

    !!! warning "Disclaimer"
        Le stime qui sono **educational** e basate su assunzioni. Non sono raccomandazioni cliniche.

    ## Se cerchi…

    | Voglio… | Vai a… |
    | --- | --- |
    | dove viene calcolato lo score | `simulator/difficulty_scoring.py` |
    | come entra nello Scenario | `simulator/models.py` (campi `difficulty_*`) |
    | modelli ODE/PK/PD e KPI | `Guides → Mathematical Models` |
    | stime PFS/OS “separate” | `Guides → Prognosis` |

    ## Terminologia

    - **DS**: difficulty score (0–100).
    - **Component score**: contributo parziale (tumor burden, frailty, ecc.).
    - **Toy estimator**: stima illustrativa non clinicamente validata.

    ## Cosa significa “difficulty”

    L’idea: alcuni scenari sono più difficili perché combinano:

    - maggiore carico tumorale / crescita più rapida
    - citogenetica ad alto rischio (resistenza/prognosi peggiore)
    - maggiore fragilità (tollerabilità più bassa)
    - stadio (R-ISS) più sfavorevole

    Il risultato è un unico numero **interpretabile** (0–100) che può alimentare UI, sorting, didattica.

    ## Pipeline (a colpo d’occhio)

    ```mermaid
    flowchart LR
      A[Scenario + Labs] --> B[Component scores]
      B --> C[Difficulty score 0..100]
      C --> D[Difficulty level]
      C --> E[Expected response (toy)]
      C --> F[Expected toxicity (toy)]
      C --> G[Expected survival (toy)]
    ```

    ## Componenti e formule

    Lo score totale è la somma di 5 componenti (pesi “in punti”):

    \[
    DS = TB + GR + CG + PF + ST
    \]

    dove:

    - \(TB \in [0,25]\) Tumor burden
    - \(GR \in [0,20]\) Growth rate
    - \(CG \in [0,25]\) Cytogenetics
    - \(PF \in [0,15]\) Patient frailty
    - \(ST \in [5,15]\) Stage (R-ISS)

    ### 1) Tumor burden (0–25)

    \[
    TB = 25 \cdot \left(1 - e^{-\lambda \cdot (T/T_{ref})}\right)
    \]

    con default: \(T_{ref}=10^{10}\) cellule, \(\lambda=1\).

    Intuizione: crescita saturante → aumentare da \(10^9\) a \(2\cdot 10^9\) impatta meno che da \(10^8\) a \(10^9\).

    ### 2) Growth rate (0–20)

    \[
    GR = 20 \cdot \mathrm{clip}\left(\frac{r}{r_{max}}, 0, 1\right)
    \]

    con default \(r_{max}=0.06/\mathrm{day}\).

    ### 3) Cytogenetics (0–25)

    Sistema additivo (con fattori “protettivi” sottrattivi) e clamp finale su \([0,25]\).

    | Lesione | Punti |
    | --- | ---: |
    | `del(17p)` | +15 |
    | `t(4;14)` | +10 |
    | `t(14;16)` | +10 |
    | `1q21 gain` | +8 |
    | hyperdiploid | −5 |
    | `t(11;14)` | −3 |

    ### 4) Frailty (0–15)

    Il punteggio è una somma di contributi:

    - età: +1/+2/+4 (soglie 70/75/80)
    - ECOG: `ecog * 2`
    - Charlson: `min(cci * 0.5, 3)`
    - funzione renale (CrCl): +1.5 / +3 (soglie 60/30)
    - albumina: +1 / +2 (soglie 3.5 / 3.0)

    Clamp finale su \([0,15]\).

    ### 5) Stage (R-ISS) (5–15)

    | R-ISS | Punti |
    | --- | ---: |
    | I | 5 |
    | II | 10 |
    | III | 15 |

    ## Interpretazione (difficulty level)

    | DS | Level |
    | ---: | --- |
    | 0–20 | Very Easy |
    | 20–40 | Easy |
    | 40–60 | Moderate |
    | 60–80 | Hard |
    | 80–100 | Very Hard |

    ## Stime “toy” derivate dallo score

    Nel modulo esistono tre stime illustrative:

    - `estimate_response_probability(DS)`
    - `estimate_toxicity_risk(DS, frailty_score)`
    - `estimate_survival_metrics(DS)`

    !!! note "Perché “toy”"
        Servono come **segnaposto trasparente** per didattica/UI: sono formule dichiarate nel codice, non modelli clinici validati.

    ### Response probability (logistic)

    Per CR/PR viene usata una forma logistica:

    \[
    p = \sigma(\beta_0 + \beta_1 \cdot DS)
    \quad\text{con}\quad
    \sigma(z)=\frac{1}{1+e^{-z}}
    \]

    La distribuzione finale impone anche un vincolo: \(P(CR) \le 0.6\cdot P(PR)\).

    ### Toxicity risk (sigmoid)

    \[
    P(\mathrm{Grade}\ge 3)=\sigma(\alpha_0+\alpha_1 DS+\alpha_2 FS)
    \]

    con \(FS\) frailty score.

    ### Survival metrics (exponential)

    \[
    S(t)=e^{-\lambda t},\quad \lambda = \lambda_0\cdot e^{\gamma\cdot DS/100}
    \]

    e:

    \[
    t_{median}=\frac{\ln 2}{\lambda}
    \]

    ## Figure (intuizione visiva)

    ![Response probabilities vs difficulty](../assets/images/models/difficulty_response_probabilities.svg)

    ![Toxicity risk vs difficulty](../assets/images/models/difficulty_toxicity_risk.svg)

    ![Survival vs difficulty](../assets/images/models/difficulty_survival_vs_score.svg)

    !!! tip "Zoom diagrammi"
        Nel sito MkDocs puoi cliccare i diagrammi Mermaid per aprire una vista fullscreen con zoom/pan.

    ## Esempio guidato (worked example)

    Esempio (numeri illustrativi):

    | Input | Valore |
    | --- | ---: |
    | tumor cells \(T\) | \(5\cdot 10^9\) |
    | growth rate \(r\) | 0.03 / day |
    | del(17p) | sì |
    | age | 76 |
    | ECOG | 2 |
    | Charlson | 3 |
    | CrCl | 55 |
    | albumin | 3.2 |
    | R-ISS | II |

    Output (componenti):

    | Componente | Range | Output |
    | --- | ---: | ---: |
    | Tumor burden | 0–25 | ~9.8 |
    | Growth rate | 0–20 | 10.0 |
    | Cytogenetics | 0–25 | 15.0 |
    | Frailty | 0–15 | ~9.0 |
    | Stage | 5–15 | 10.0 |
    | **Total DS** | **0–100** | **~43.8 (Moderate)** |

    ## Mapping al codice

    - component scores: `TumorBurdenScore`, `GrowthRateScore`, `CytogeneticScore`, `PatientFrailtyScore`, `StageScore`
    - composizione: `DifficultyScoreCalculator`
    - stime toy: `estimate_response_probability`, `estimate_toxicity_risk`, `estimate_survival_metrics`

=== "EN"
    This page explains the **mathematical model** that produces:

    - `difficulty_score` (0–100) and `difficulty_level` (Very Easy → Very Hard)
    - “toy” estimates of **response**, **toxicity risk**, and **survival** (educational)

    Primary source: `simulator/difficulty_scoring.py`.

    !!! warning "Disclaimer"
        These estimates are **educational** and assumption-based. They are not clinical recommendations.

    ## Quick find

    | I want… | Go to… |
    | --- | --- |
    | where the score is computed | `simulator/difficulty_scoring.py` |
    | how it lands on Scenario | `simulator/models.py` (`difficulty_*` fields) |
    | ODE/PK/PD model and KPIs | `Guides → Mathematical Models` |
    | separate PFS/OS estimates | `Guides → Prognosis` |

    ## Terminology

    - **DS**: difficulty score (0–100).
    - **Component score**: partial contribution (tumor burden, frailty, etc.).
    - **Toy estimator**: illustrative estimator not clinically validated.

    ## What “difficulty” means

    Some scenarios are harder because they combine:

    - higher tumor burden / faster growth
    - high-risk cytogenetics (more resistance / worse prognosis)
    - higher frailty (lower tolerance)
    - worse stage (R-ISS)

    The result is a single **interpretable** number (0–100) used for UI, sorting, and education.

    ## Pipeline (at a glance)

    ```mermaid
    flowchart LR
      A[Scenario + Labs] --> B[Component scores]
      B --> C[Difficulty score 0..100]
      C --> D[Difficulty level]
      C --> E[Expected response (toy)]
      C --> F[Expected toxicity (toy)]
      C --> G[Expected survival (toy)]
    ```

    ## Components and formulas

    Total score is the sum of 5 components (weights are “points”):

    \[
    DS = TB + GR + CG + PF + ST
    \]

    where:

    - \(TB \in [0,25]\) tumor burden
    - \(GR \in [0,20]\) growth rate
    - \(CG \in [0,25]\) cytogenetics
    - \(PF \in [0,15]\) patient frailty
    - \(ST \in [5,15]\) stage (R-ISS)

    ### 1) Tumor burden (0–25)

    \[
    TB = 25 \cdot \left(1 - e^{-\lambda \cdot (T/T_{ref})}\right)
    \]

    defaults: \(T_{ref}=10^{10}\) cells, \(\lambda=1\).

    ### 2) Growth rate (0–20)

    \[
    GR = 20 \cdot \mathrm{clip}\left(\frac{r}{r_{max}}, 0, 1\right)
    \]

    default \(r_{max}=0.06/\mathrm{day}\).

    ### 3) Cytogenetics (0–25)

    Additive points (with protective subtractive factors), clamped to \([0,25]\).

    | Lesion | Points |
    | --- | ---: |
    | `del(17p)` | +15 |
    | `t(4;14)` | +10 |
    | `t(14;16)` | +10 |
    | `1q21 gain` | +8 |
    | hyperdiploid | −5 |
    | `t(11;14)` | −3 |

    ### 4) Frailty (0–15)

    Sum of contributions:

    - age: +1/+2/+4 (thresholds 70/75/80)
    - ECOG: `ecog * 2`
    - Charlson: `min(cci * 0.5, 3)`
    - renal function (CrCl): +1.5 / +3 (thresholds 60/30)
    - albumin: +1 / +2 (thresholds 3.5 / 3.0)

    Final clamp to \([0,15]\).

    ### 5) Stage (R-ISS) (5–15)

    | R-ISS | Points |
    | --- | ---: |
    | I | 5 |
    | II | 10 |
    | III | 15 |

    ## Difficulty level

    | DS | Level |
    | ---: | --- |
    | 0–20 | Very Easy |
    | 20–40 | Easy |
    | 40–60 | Moderate |
    | 60–80 | Hard |
    | 80–100 | Very Hard |

    ## “Toy” estimates derived from the score

    The module includes three illustrative estimators:

    - `estimate_response_probability(DS)`
    - `estimate_toxicity_risk(DS, frailty_score)`
    - `estimate_survival_metrics(DS)`

    !!! note "Why “toy”"
        They are a **transparent placeholder** for UI/education: formulas are explicit in code but not clinically validated.

    ### Response probability (logistic)

    \[
    p = \sigma(\beta_0 + \beta_1 \cdot DS),\quad
    \sigma(z)=\frac{1}{1+e^{-z}}
    \]

    The final distribution also enforces \(P(CR) \le 0.6\cdot P(PR)\).

    ### Toxicity risk (sigmoid)

    \[
    P(\mathrm{Grade}\ge 3)=\sigma(\alpha_0+\alpha_1 DS+\alpha_2 FS)
    \]

    where \(FS\) is the frailty score.

    ### Survival metrics (exponential)

    \[
    S(t)=e^{-\lambda t},\quad \lambda = \lambda_0\cdot e^{\gamma\cdot DS/100}
    \]

    and:

    \[
    t_{median}=\frac{\ln 2}{\lambda}
    \]

    ## Figures (visual intuition)

    ![Response probabilities vs difficulty](../assets/images/models/difficulty_response_probabilities.svg)

    ![Toxicity risk vs difficulty](../assets/images/models/difficulty_toxicity_risk.svg)

    ![Survival vs difficulty](../assets/images/models/difficulty_survival_vs_score.svg)

    !!! tip "Zoom diagrams"
        On the MkDocs site you can click Mermaid diagrams to open a fullscreen view with zoom/pan.

    ## Worked example

    Example (illustrative numbers):

    | Input | Value |
    | --- | ---: |
    | tumor cells \(T\) | \(5\cdot 10^9\) |
    | growth rate \(r\) | 0.03 / day |
    | del(17p) | yes |
    | age | 76 |
    | ECOG | 2 |
    | Charlson | 3 |
    | CrCl | 55 |
    | albumin | 3.2 |
    | R-ISS | II |

    Output (components):

    | Component | Range | Output |
    | --- | ---: | ---: |
    | Tumor burden | 0–25 | ~9.8 |
    | Growth rate | 0–20 | 10.0 |
    | Cytogenetics | 0–25 | 15.0 |
    | Frailty | 0–15 | ~9.0 |
    | Stage | 5–15 | 10.0 |
    | **Total DS** | **0–100** | **~43.8 (Moderate)** |

    ## Code mapping

    - component scores: `TumorBurdenScore`, `GrowthRateScore`, `CytogeneticScore`, `PatientFrailtyScore`, `StageScore`
    - composition: `DifficultyScoreCalculator`
    - toy estimators: `estimate_response_probability`, `estimate_toxicity_risk`, `estimate_survival_metrics`
