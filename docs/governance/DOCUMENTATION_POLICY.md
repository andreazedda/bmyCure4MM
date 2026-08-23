---
title: Documentation Policy
status: CANONICAL
owner: Andrea Zedda
audience: document authors and reviewers
last_verified_at: 2026-08-24
last_verified_git_sha: bf097810b337dc6b766cda04497005670cd96513
---

# Documentation policy

## Required document statuses

Every governed document must use one status:

```text
CURRENT_VERIFIED
CURRENT_PARTIAL
CURRENT_HEURISTIC_PROTOTYPE
CURRENT_HYPOTHETICAL_PROTOTYPE
TARGET_APPROVED
CANDIDATE_COMPARATOR
RESEARCH_HYPOTHESIS
DEFERRED
REJECTED
HISTORICAL_ARCHIVE
```

## Required metadata

Current and target documents must declare:

```yaml
title:
status:
owner:
audience:
last_verified_at:
last_verified_git_sha:
source_of_truth:
supersedes:
superseded_by:
```

Fields that do not apply may be omitted, but `status`, `owner`, `last_verified_at` and `last_verified_git_sha` are mandatory for current documents.

## Current versus target

- `CURRENT_*` content must map to an existing code/schema/registry entry point.
- `TARGET_APPROVED` content must state prerequisites and acceptance evidence.
- `RESEARCH_HYPOTHESIS` content must include competing hypotheses and a falsification path.
- historical content must begin with an archive banner and must not appear in current navigation.

## Formula standard

Every formula that affects project interpretation must include:

1. scientific or computational question;
2. equation;
3. definition and unit of every variable;
4. parameter provenance/evidence status;
5. source-code entry point and model version;
6. worked numerical example;
7. valid domain and assumptions;
8. sensitivity/uncertainty note;
9. falsification or rejection criterion;
10. allowed and forbidden conclusions.

A formula in a DOCX must use native Word equation objects. The canonical repository representation remains Markdown/LaTeX plus code mapping.

## Diagram standard

- label every diagram `CURRENT` or `TARGET`;
- version the source;
- prefer several small diagrams over one dense diagram;
- use no more than approximately eight primary nodes per operational view;
- avoid crossing connectors where practical;
- use landscape/wider layout for complex flows;
- keep labels readable at normal page zoom;
- render Mermaid during validation;
- provide a caption, scope and interpretation note;
- do not place floating overlays over text in exported office documents.

## Claims standard

Current documents inherit `INTENDED_USE.md` and `CLAIMS_POLICY.md`. Archived claims cannot be reintroduced through examples, screenshots, translations or generated API documentation.

## Review cadence

- verify current-state pages whenever a release or model/schema change affects them;
- verify at least every 90 days while the project is active;
- update the verified Git SHA after review;
- create an issue when a document cannot be verified rather than asserting currency.

## Inventory and navigation

`docs/DOCUMENTATION_INVENTORY.yaml` classifies the documentation corpus. Only `CANONICAL`, `CURRENT_SUPPORTING` and explicitly selected `TARGET_SUPPORTING` documents may appear in current MkDocs navigation.

Unmatched files resolve to `UNCLASSIFIED_REQUIRES_REVIEW` and must stay outside current navigation until reviewed.
