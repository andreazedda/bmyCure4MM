---
title: Real-Patient Research Loop v0.1
status: TARGET_APPROVED
owner: Andrea Zedda
audience: M1 contributors and reviewers
last_verified_at: 2026-08-24
last_verified_git_sha: bf097810b337dc6b766cda04497005670cd96513
canonical_tracker: https://github.com/andreazedda/bmyCure4MM/issues/28
---

# Real-Patient Research Loop v0.1

> **TARGET M1 workflow.** The repository has many prerequisites, but this complete loop has not yet been demonstrated.

## Flow

```text
source-verified evidence
→ immutable dataset
→ frozen computational input
→ baseline Twin
→ pre-calibration residuals
→ calibration accepted or rejected
→ exposure-aware scenarios
→ toxicity attribution with competing explanations
→ temporal backtest against baselines
→ uncertainty and identifiability report
→ immutable research report and manifest
```

## Completed prerequisites

- Twin and computational-input lineage;
- E1 intended use and epistemic output language;
- deterministic dependency graph;
- immutable run manifests, model registry, artifact integrity and comparability.

## Remaining milestone prerequisite

M0-R must close through documentation, protected CI and minimum QA issues `#14 → #13 → #23`.

Private read-only source verification may prepare in parallel. New canonical calibration and patient-specific what-if claims cannot be promoted before M0-R closes.

## Scientific work order

```text
#32 source assertions and dataset
→ #33 observation model
→ #34 day-resolved lenalidomide exposure
→ #35 toxicity attribution
→ #36 calibration, identifiability and temporal validation
→ #18 immutable end-to-end loop
```

## Acceptance is not positivity

Valid outcomes include:

```text
CALIBRATION_ACCEPTED
INSUFFICIENT_DATA
CALIBRATION_FAILED
MODEL_NOT_IDENTIFIABLE
NO_OUT_OF_SAMPLE_GAIN
BACKTEST_UNACCEPTABLE
```

## Forbidden conclusions

- an alternative regimen would have produced a superior clinical outcome;
- a mechanistic trajectory is an identified causal effect;
- optimizer completion proves validity;
- unresolved source ambiguity may be silently imputed;
- in-sample RMSE improvement alone establishes patient-specific prediction.
