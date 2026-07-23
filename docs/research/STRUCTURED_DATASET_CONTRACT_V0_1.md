# bmyCure4MM Structured Research Dataset Contract v0.1

## Status

Schema version: `0.1.0`

Epistemic level: research prototype.

This contract defines the private structured longitudinal dataset consumed by
bmyCure4MM. It does not define a clinical recommendation format.

## Core invariant

Clinical evidence must flow through:

document/source
→ atomic source assertion
→ structured research record
→ validated database import
→ model/run provenance

A value without traceable source metadata is incomplete research data.

## Privacy boundary

Real or potentially identifying source material belongs outside Git.

Git may contain:

- schema;
- generic validation logic;
- synthetic test fixtures;
- pseudonymous metadata that cannot identify a patient.

Git must not contain:

- patient names;
- MRNs derived from real clinical systems;
- birth dates tied to real persons;
- clinical source PDFs;
- private clinical excerpts;
- private dataset payloads.

Private datasets should live below `local_private/`.

## Dataset identity

Every private dataset has:

- `schema_version`;
- `dataset_id`;
- `dataset_version`;
- `case_ref`;
- `created_at`;
- source-document registry;
- atomic records.

`dataset_id + dataset_version` identifies an immutable dataset snapshot.

A correction to a clinical value, unit, date, interpretation, source mapping,
or validation status produces a new dataset version.

Previously used dataset versions must remain reproducible.

## Record identity

The database identity of a fact is distinct from the identity of the evidence
supporting that fact.

Examples:

### Longitudinal laboratory result

Semantic identity:

`case_ref + date + analyte`

### Adverse event

Semantic identity:

`case_ref + date + event_type`

### Therapy interruption

Semantic identity:

`case_ref + drug + start_date + end_date`

### Therapy course

Semantic identity:

`case_ref + regimen_ref + start_date`

### Assessment

Semantic identity:

`case_ref + date`

These keys prevent duplicate semantic records.

They do NOT replace provenance.

## Atomic provenance

Every structured record must contain one or more provenance assertions.

Each assertion contains:

- `assertion_id`;
- `source_id`;
- `source_locator`;
- `source_sha256`;
- `extraction_method`;
- `technical_validation_status`;
- `clinical_validation_status`;
- `source_quality`;
- `epistemic_class`;
- optional notes.

`source_locator` may be a page, section, table, row, or another stable
location inside the original source.

`source_sha256` identifies the exact source artifact without storing that
artifact in Git.

## Extraction states

Allowed extraction methods:

- `manual_verified`
- `ai_suggested_human_verified`
- `legacy_manual_extraction`

AI-only extraction is not sufficient for a validated clinical observation.

## Technical validation

Allowed statuses:

- `validated`
- `pending`
- `failed`
- `legacy_unverified`

Technical validation answers:

"Was the structured value faithfully transferred from the source?"

## Clinical validation

Allowed statuses:

- `not_required`
- `pending`
- `validated`
- `disputed`

Clinical validation answers questions requiring medical interpretation.

A raw laboratory value generally does not need clinical interpretation to
exist as an observation.

A causal attribution or interpretation may.

## Epistemic classes

Allowed classes:

- `observed`
- `derived`
- `inferred`
- `simulated`
- `hypothetical`

Clinical source observations should normally be `observed`.

Derived values must never silently become observed values.

## Import idempotency

For an unchanged dataset version:

first import:
`created >= 0`

second import:
`created = 0`
`changed = 0`
`unchanged = N`

A rerun must not rewrite unchanged semantic records merely because an upsert
operation was invoked.

If a semantic record differs from the stored version, the importer must report
the difference explicitly.

## Conflict handling

Two sources may support the same semantic observation.

They must not automatically create duplicate clinical facts.

Their provenance assertions should coexist.

If two sources disagree on value, unit, date, or interpretation, the importer
must classify the record as a conflict and must not silently overwrite the
previous observation.

## Dataset hash

Every immutable private dataset snapshot must have a SHA-256 hash computed from
its canonical serialized representation.

Research runs must eventually record at least:

- dataset ID;
- dataset version;
- dataset hash;
- Git commit;
- model version;
- configuration hash;
- random seed where relevant.

## Research Loop v0.1 gate

The dataset layer is ready only when all are true:

1. schema validation passes;
2. no direct identifiers are exported;
3. required units and dates are valid;
4. each clinical fact has atomic provenance;
5. validation states are explicit;
6. import is semantically idempotent;
7. conflicts are surfaced rather than overwritten;
8. the dataset has an immutable version and hash.
