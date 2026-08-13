# Research run manifests and model registry

## Purpose

Every persisted scientific execution is bound to an immutable
`ResearchRunManifest`. A successful process exit is not a scientific result;
the manifest establishes identity and integrity only. Calibration,
identifiability, validation, and clinical meaning remain separate gates.

The implementation is split across:

- `twin_engine/run_manifest.py`: identity construction and artifact checks;
- `twin_engine/model_registry.json`: governed model registrations;
- `twin_engine/contracts.py` and `twin_engine/schemas/`: versioned scientific
  JSON contracts;
- `twin_engine/comparability.py`: direct-comparison and invalidation rules;
- `ResearchRunManifest`, `ResearchRunInvalidation`, and
  `SimulationRunMetadata`: persisted control-plane records.

## Mandatory version vector

Each manifest binds:

| Identity | Meaning |
| --- | --- |
| `app_version`, `api_version`, `db_schema_version` | application and schema identity |
| `dataset_id`, `dataset_version`, `dataset_sha256` | immutable dataset binding |
| `record_subset_sha256` | exact computational subset binding |
| `model_id`, `model_version`, `model_card_version` | governed model identity |
| `configuration_sha256` | solver/configuration identity |
| `evidence_graph_version`, `validation_protocol_version` | evidence and validation contracts |
| `report_template_version` | interpretation-template identity |
| `git_sha`, `dependency_lock_sha256`, `model_registry_sha256` | code and dependency identity |
| `container_digest` | OCI image identity, or explicit `UNAVAILABLE_LOCAL` |
| `random_seed`, `created_at` | stochastic and temporal identity |
| `intended_use_level`, `epistemic_label` | governed interpretation boundary |

Missing dataset provenance is recorded as `UNBOUND`; it is never replaced by
an inferred or fabricated dataset identifier. The record-subset hash remains
mandatory. Local runs without an OCI identity record `UNAVAILABLE_LOCAL`.

## Artifact integrity

Inputs and outputs always receive canonical SHA-256 references. Persisted JSON
reports and trajectories additionally receive `media://` references and file
hashes. The matching manifest artifact is written under
`MEDIA_ROOT/research_manifests/`.

`verify_manifest_artifacts()` fails closed when:

- the database manifest and manifest artifact differ;
- a file is missing or its SHA-256 differs;
- a path escapes `MEDIA_ROOT`;
- a hash-only reference is malformed.

Report views return `ARTIFACT_INTEGRITY_FAILED` instead of rendering a
manifest-bound artifact that fails verification. Legacy unversioned reports
remain readable but are labelled as requiring rerun for direct comparison.

Manifest artifacts contain hashes, version identifiers, and relative artifact
references only. They must never contain names, MRNs, dates of birth, free-text
clinical source excerpts, or private dataset payloads.

## Immutability and diagnostics

Manifests and manifest-bound metadata cannot be updated or deleted through the
model API. Completed counterfactual records cannot be mutated. Uncertainty,
sensitivity, and robustness diagnostics therefore create new manifests linked
to the source run; they do not append fields to a completed result.

Historical rows created before the contract are explicitly tagged
`legacy-unversioned-v0` by migration `0004`. They are not silently promoted to
the current schemas.

## Direct comparability

`compare_manifests()` returns one of:

- `DIRECTLY_COMPARABLE`;
- `NOT_DIRECTLY_COMPARABLE_WITHOUT_RERUN`.

Direct comparison requires all of the following:

1. neither run has an invalidation record;
2. model IDs and declared compatible major versions match;
3. the model-registry identity matches, binding schemas, endpoint definitions,
   parameter definitions, and units;
4. dataset identity is bound (not `UNBOUND`/`UNAVAILABLE`) and dataset version,
   physical dataset hash, and record-subset hash match;
5. intended-use levels match.

A user may still inspect runs that fail this gate, but must not present their
numerical difference as a direct like-for-like comparison.

## Invalidation

Corrections never rewrite old manifests. `invalidate_manifest()` appends an
immutable `ResearchRunInvalidation` containing the change type, previous and
replacement identities, reason, and a canonical change hash.

Governed change types cover clinical values, datasets, unit mappings, models,
parameter defaults, dependencies, and configuration. The corrected execution
must produce a new run manifest. The old artifact remains available for audit
and rollback.

## Model registry change control

Each registered model declares its model/version identifiers, input/output
schema versions, endpoint-definition version, parameter names, units, evidence
status, implementation entry point, model card, compatible major versions, and
explicitly invalidated prior versions.

Changes to equations, endpoint meaning, units, input/output semantics, or
parameter meaning require scientific review and may require a `model!:` change
and major `model_version` increment. A registry edit changes
`model_registry_sha256`, so prior and new runs are not treated as directly
comparable without explicit compatibility evidence and rerun.

## Rollback

Code rollback is a PR revert plus reverse migration rehearsal on a verified
database copy. Real database migration requires a checksum-verified backup and
cardinality inventory first. Run artifacts and invalidation records are retained
for audit; rollback does not rewrite scientific history.
