---
title: Current Model Registry and Status
status: CURRENT_VERIFIED
owner: Andrea Zedda
audience: model developers, reviewers and research users
last_verified_at: 2026-08-24
last_verified_git_sha: bf097810b337dc6b766cda04497005670cd96513
machine_readable_source: twin_engine/model_registry.json
---

# Current model registry and status

The machine-readable registry is authoritative for IDs, versions, schemas and entry points. This page adds formulas, examples and interpretation boundaries.

## Status summary

| Model ID | Version | Entry point | Status |
|---|---|---|---|
| `patient_twin_state_model` | `research-twin-v1` | `twin_engine.simulation_bridge.run_patient_simulation` | CURRENT_HEURISTIC_PROTOTYPE |
| `observation_model` | `observation-model-v1` | `twin_engine.observation_model.predict_biomarkers` | CURRENT_HEURISTIC_PROTOTYPE |
| `lenalidomide_exposure_model` | `exposure-bridge-v1` | `twin_engine.exposure_bridge.build_exposure_dose_function` | CURRENT_EXECUTABLE |
| `hepatic_signal_model` | `toxicity-prototype-v1` | `twin_engine.toxicity_dynamics.compute_toxicity_dynamics` | CURRENT_HEURISTIC_PROTOTYPE |
| `neutropenia_signal_model` | `toxicity-prototype-v1` | `twin_engine.toxicity_dynamics.compute_toxicity_dynamics` | CURRENT_HEURISTIC_PROTOTYPE |
| `counterfactual_model` | `research-twin-v1` | `twin_engine.counterfactual.run_counterfactual` | CURRENT_HEURISTIC_PROTOTYPE |
| `therapy_design_toy_model` | `therapy-design-toy-v1` | `simulator.design.reporting.run_design_report` | CURRENT_HYPOTHETICAL_PROTOTYPE |

## 1. Patient Twin state / mechanistic runtime

### Question

How does the current abstract tumour/healthy system evolve under a configured therapy schedule?

### Core equations

$$
\frac{dT}{dt} = r_T T\left(1-\frac{T}{K_T}\right)-\left(\sum_i E_i(C_i)\right)T
$$

$$
\frac{dH}{dt} = r_H H\left(1-\frac{H}{K_H}\right)-I\,\overline{E(C)}\,H
$$

$$
\frac{dC_i}{dt}=-k_{elim,i}C_i+u_i(t)
$$

$$
E_i(C_i)=E_{max,i}\frac{C_i}{EC50_i+C_i}
$$

Variables:

- $T$: abstract tumour-cell burden, cells;
- $H$: abstract healthy-cell compartment, cells;
- $C_i$: model concentration for drug $i$, model units;
- $r_T,r_H$: growth rates, day$^{-1}$;
- $K_T,K_H$: carrying capacities, cells;
- $I$: dimensionless immune-compromise/toxicity multiplier;
- $u_i(t)$: schedule-resolved model input;
- $k_{elim,i}=\ln(2)/t_{1/2,i}$, day$^{-1}$;
- $E_{max,i}$: maximum model effect;
- $EC50_i$: model concentration at half maximum effect.

Parameters in the current Twin registry are primarily `HEURISTIC`. The equations are mechanistic research abstractions, not a validated MM natural-history model.

### Worked example

With $T=10^9$, $K_T=10^{12}$ and $r_T=0.02$ day$^{-1}$, untreated logistic growth is approximately:

$$
0.02\times10^9\left(1-10^{-3}\right)\approx1.998\times10^7\ \text{cells/day}
$$

If the total model drug effect is $0.03$ day$^{-1}$, the kill term is $3.0\times10^7$ cells/day and the instantaneous derivative becomes negative. This illustrates model behavior; it does not predict a patient's tumour-cell count.

### Falsification

Reject or revise the current model for a task when it fails declared numerical invariants, cannot improve on simple temporal baselines, produces non-identifiable parameters or does not explain held-out observations within the validation protocol.

## 2. Observation model

### Question

How are latent $T$ and $H$ mapped to three currently supported observed biomarkers?

$$
\widehat{M}=\alpha_M T+\beta_M
$$

$$
\widehat{F}=\alpha_F\left(\frac{T}{T_{ref}}\right)^{\gamma_F}
$$

$$
\widehat{Hb}=Hb_0\left(\frac{H}{H_{ref}}\right)^{\eta_H}
$$

- $\widehat{M}$: predicted M-protein, g/dL;
- $\widehat{F}$: predicted FLC ratio, dimensionless;
- $\widehat{Hb}$: predicted hemoglobin, g/dL;
- $\alpha_M,\beta_M,\alpha_F,\gamma_F,Hb_0,T_{ref},H_{ref},\eta_H$: heuristic or calibration parameters.

Default values are defined in `twin_engine/observation_model.py` and can be initialized from a baseline assessment.

### Worked example

If $\alpha_M=10^{-9}$ g/dL/cell, $\beta_M=0$ and $T=8\times10^8$ cells, then:

$$
\widehat{M}=10^{-9}\times8\times10^8=0.8\ \text{g/dL}
$$

This result is model-relative. FLC ratio is not a universal direct tumour-burden scalar, and hemoglobin is affected by factors not represented by $H$ alone.

### Falsification

The mapping is rejected for patient-specific prediction if residuals show systematic bias, parameters are non-identifiable or a simple baseline performs as well or better on future observations.

## 3. Exposure bridge

### Question

Does the system preserve the temporal identity of a treatment schedule?

For a horizon of $N$ days with administered dose $d_j$ on day $j$:

$$
D_{cum}(n)=\sum_{j=0}^{n} d_j
$$

$$
\bar d=\frac{1}{N}\sum_{j=0}^{N-1}d_j
$$

The profile stores the daily vector, cumulative vector, peak dose, interruption days, schedule type and a canonical profile hash.

### Worked example

```text
5 mg daily:         [5, 5, 5, 5]
10 mg alternate:   [10, 0, 10, 0]
```

Both average 5 mg/day, but mean absolute temporal distance is:

$$
\frac{|5-10|+|5-0|+|5-10|+|5-0|}{4}=5\ \text{mg/day}
$$

Therefore `same_average_exposure=true` and `different_temporal_profile=true`.

### Limitation

This is administered-dose exposure identity, not patient-specific plasma concentration or PK validation.

## 4. Hepatic and neutropenia signal models

### Question

Can exposure and prior history generate normalized prototype risk signals for stress testing?

The current implementation uses clipped heuristic combinations. A simplified day-$t$ representation is:

$$
S_{AST,t}=clip_{[0,1]}(b_{AST}+\alpha_{AST}d_t^*+\beta_s S+\gamma_{AST}P_{AST})
$$

$$
S_{ALT,t}=clip_{[0,1]}(b_{ALT}+\alpha_{ALT}d_t^*+\beta_s S+\gamma_{ALT}P_{ALT})
$$

$$
S_{liver,t}=\max(S_{AST,t},S_{ALT,t})
$$

$$
S_{neu,t}=clip_{[0,1]}(b_{neu}+\alpha_{neu}D_{cum,t}^*+\beta_{prior}N_{prior})
$$

where normalized daily and cumulative dose are referenced to 10 mg and 100 mg, respectively. Coefficients are defaults or limited single-patient estimates.

### Worked example

Using $b_{neu}=0.15$, $\alpha_{neu}=0.45$, normalized cumulative exposure $0.6$, $\beta_{prior}=0.25$ and prior-neutropenia flag $1$:

$$
S_{neu}=clip(0.15+0.45\times0.6+0.25)=0.67
$$

The value `0.67` is a normalized prototype signal, not a predicted ANC or probability of clinical neutropenia.

### Falsification

Do not promote this model to attribution or prediction unless it is source-verified, compared against competing explanations, calibrated and validated prospectively or temporally under a governed protocol.

## 5. Counterfactual model

### Question

How do a baseline and alternative mechanistic configuration differ under the same governed Twin state?

The current utility is heuristic:

$$
U_1=E_{tumour}+(1-L_{healthy})+D_{durability}-P_{constraint}
$$

A second version adds prototype toxicity penalties:

$$
U_2=U_1-\lambda_L S_{liver}-\lambda_N S_{neu}
$$

Current default weights are $\lambda_L=\lambda_N=0.5$.

### Worked example

If tumour reduction is $0.70$, healthy loss $0.10$, durability $0.80$, constraint penalty $0.15$, liver signal $0.30$ and neutropenia signal $0.40$:

$$
U_1=0.70+0.90+0.80-0.15=2.25
$$

$$
U_2=2.25-0.5(0.30)-0.5(0.40)=1.90
$$

This score ranks model configurations only. It does not identify clinical benefit or a causal effect.

## 6. Therapy-design toy model

This experimental model exercises target, modality, logic, bulk/reservoir/subclone and Pareto software abstractions. Its antigens, coefficients, escape and relapse proxies are hypothetical/heuristic. It is useful for architecture tests and falsifiable research design, not for selecting or claiming efficacy of ADC, CAR-T, bispecific or logic-gated therapies.

The target M4/M5 models require evidence-linked parameterization, competing hypotheses, external comparators, sensitivity analysis and validation before their status can be upgraded.

## Version and invalidation rule

Any equation, endpoint, unit, parameter default or dependency change that alters numerical output requires a model-version decision and an output-diff report. Prior runs are never silently reinterpreted.
