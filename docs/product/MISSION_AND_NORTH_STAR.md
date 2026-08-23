---
title: Mission and North Star
status: CANONICAL_STRATEGY
owner: Andrea Zedda
audience: researchers, clinicians, biologists, pharmacologists, bioengineers, data scientists, contributors
last_verified_at: 2026-08-24
last_verified_git_sha: bf097810b337dc6b766cda04497005670cd96513
source_of_truth: docs/product/MISSION_AND_NORTH_STAR.md
---

# Mission and North Star

## Mission

bmyCure4MM is intended to become a **virtual and reproducible Multiple Myeloma research laboratory**.

Its purpose is to make Multiple Myeloma research questions computable, inspectable and falsifiable. The platform should allow interdisciplinary users to:

1. acquire structured knowledge about disease biology, measurement, pharmacology, current treatments, innovative modalities, uncertainty and evidence;
2. construct lineage-bound virtual patients and synthetic cohorts;
3. define experiments using established regimens, alternative schedules, combinations, sequences, novel molecules and experimental modalities;
4. execute mechanistic simulations under explicit assumptions;
5. compare efficacy signals, toxicity, residual disease, resistance, antigen escape, uncertainty and model failure;
6. reject weak or non-identifiable hypotheses before committing to more expensive experimental work;
7. translate computational results into source-linked, reproducible and testable research protocols.

The platform is not currently a treatment-selection system. It does not provide validated patient-specific benefit prediction, identified causal effects or clinical decision support.

## Why this architecture is required

The architecture follows from the research loop:

```text
Knowledge
→ source assertions
→ versioned data
→ virtual patient
→ experiment definition
→ simulation
→ comparison
→ falsification
→ research protocol
→ new evidence
```

A database alone cannot execute this loop. A simulator without provenance cannot make results reproducible. A Patient Twin without observation, calibration and validation contracts cannot support credible patient-specific research. An optimizer without uncertainty and evidence penalties can rank invalid hypotheses. A knowledge portal without links to executable models cannot close the loop from literature to experiment.

## Research objective

The general decision problem is:

$$
a^* = \arg\max_{a \in \mathcal{A}} U(a \mid D, M, C)
$$

where:

- $a$ is a research action: regimen, schedule, sequence, combination, molecule, modality or experiment;
- $\mathcal{A}$ is the set of admissible research actions;
- $D$ is the versioned evidence and dataset available for the experiment;
- $M$ is the explicit model set and its assumptions;
- $C$ is the set of biological, privacy, security, feasibility and intended-use constraints;
- $U$ is a declared research utility, not a clinical recommendation.

A representative future utility may combine:

$$
U(a) =
w_B E_{bulk}(a)
+ w_R E_{reservoir}(a)
+ w_I E_{immune}(a)
+ w_O R_{organ}(a)
- w_T T_{tox}(a)
- w_L P_{relapse}(a)
- w_U U_{uncertainty}(a)
- w_F C_{infeasible}(a)
$$

Every term must have a versioned definition, measurable proxy, uncertainty statement and falsification condition before it can be used as evidence.

### Example

Suppose two hypothetical schedules have the same average dose. Schedule A administers 5 mg every day; schedule B administers 10 mg every other day. Their average daily dose may be equal, but their time-resolved exposure vectors differ. The platform must therefore preserve schedule identity, simulate both profiles and report whether the result depends on temporal exposure rather than collapse them into one scalar.

## Knowledge Commons

The Knowledge Commons is a first-class target capability. It should connect:

```text
claim
↔ source
↔ evidence level
↔ biological entity
↔ model parameter
↔ executable experiment
↔ result
↔ limitation
```

It must support progressive entry into the field without lowering epistemic standards. A new contributor should be able to learn the relevant terminology, inspect the evidence behind a claim, see how that claim affects a model, and understand which conclusions remain forbidden.

## Primary users

- **Hematology researchers:** formulate and inspect longitudinal and mechanistic hypotheses.
- **Clinicians in research roles:** examine data readiness, assumptions, model limitations and candidate experiments without receiving prescriptive output.
- **Biologists and pharmacologists:** connect targets, pathways, resistance and toxicity mechanisms to experiments.
- **Bioengineers and quantitative scientists:** implement state estimation, PK/PD, uncertainty, control and validation methods.
- **Data scientists:** build governed adapters, benchmarks, predictive comparators and diagnostics.
- **Students and new contributors:** follow a structured learning path into Multiple Myeloma research.

## Scientific success criteria

A successful result is not necessarily a positive result. Valid outcomes include:

```text
INSUFFICIENT_DATA
CALIBRATION_FAILED
MODEL_NOT_IDENTIFIABLE
NO_OUT_OF_SAMPLE_GAIN
BACKTEST_UNACCEPTABLE
HYPOTHESIS_NOT_SUPPORTED
```

The platform creates value when it reduces unsupported confidence, exposes missing data, distinguishes competing mechanisms and produces reproducible next experiments.

## Non-goals for the current roadmap

- automated treatment choice;
- patient-specific dosing instructions;
- claims of clinical superiority;
- causal treatment-effect claims without identification;
- regulated or SaMD use;
- representation of a target model as implemented merely because it appears in a design document.
