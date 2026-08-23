---
title: Twin and Computational-Input Lineage
status: CURRENT_VERIFIED
owner: Andrea Zedda
audience: model developers and research reviewers
last_verified_at: 2026-08-24
last_verified_git_sha: bf097810b337dc6b766cda04497005670cd96513
---

# Twin and computational-input lineage

## Purpose

A Twin state is interpretable only when its data subset, model, configuration, code identity and parent state are explicit.

## Computational-input identity

The current contract binds a governed structured subset to model/configuration/code identity and produces content-addressed input metadata.

A change to a governed input must produce a new computational-input identity. Historical states without the current lineage contract remain explicitly legacy/unbound; they are not silently reinterpreted.

## Twin lineage

A lineage-bound state records at least:

```text
lineage schema/version
computational-input identity
model version
configuration hash
code/Git identity
dataset identity/subset hash
parent Twin state when applicable
calibration relationship
```

Calibration creates a child state; it does not mutate the parent state into a different scientific object.

## Readiness invariant

Residuals or calibration artifacts from another lineage cannot satisfy readiness for the current Twin.

```text
same patient
≠ same computational input
≠ same Twin lineage
```

## Historical data

The lineage schema migration preserved historical states with an empty/legacy lineage representation. Such states remain readable but cannot be treated as current lineage-bound evidence without an explicit new state.

## Scientific boundary

Valid lineage proves identity and traceability. It does not prove:

- that the state variables are biologically correct;
- that parameters are identifiable;
- that calibration is accepted;
- that future prediction improves over a baseline;
- that a counterfactual is causal.
