---
title: Run Comparability and Invalidation
status: CURRENT_VERIFIED
owner: Andrea Zedda
audience: research users and model maintainers
last_verified_at: 2026-08-24
last_verified_git_sha: bf097810b337dc6b766cda04497005670cd96513
---

# Run comparability and invalidation

## Immutable runs

A completed research run is not edited. A change in data, code, model, configuration, evidence, endpoint definition or seed creates a new run and a new manifest.

## Version vector

A run identity contains the governed application/API/schema, dataset and subset identity, model and schema versions, configuration hash, evidence/validation/report versions, Git SHA, container identity where available, random seed, intended-use level and epistemic label.

## Direct comparability

Direct comparison is computed, not assumed. At minimum the system evaluates compatibility of:

```text
model major version
input/output schema
endpoint definitions
units
intended-use class
data compatibility
```

When these conditions fail, results are marked:

```text
NOT_DIRECTLY_COMPARABLE_WITHOUT_RERUN
```

## Example

Two runs with identical data and code but different definitions of `toxicity_score` are not directly comparable, even if both values are between zero and one. The endpoint-semantic change requires a model/endpoint version decision and usually a rerun.

## Invalidation

Corrections or changes that may invalidate prior runs include:

- a corrected clinical value or unit;
- a changed dataset subset;
- a model equation or parameter default;
- changed endpoint semantics;
- a numerical dependency upgrade that alters output;
- an artifact hash mismatch;
- an intended-use change.

Invalidation is append-only and records the reason; prior artifacts are not silently overwritten.

## Boundary

Comparability and integrity prove whether outputs may be contrasted under the contract. They do not prove that either model is scientifically valid.
