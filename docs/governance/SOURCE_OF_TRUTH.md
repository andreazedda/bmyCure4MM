---
title: Source of Truth
status: CANONICAL
owner: Andrea Zedda
audience: all contributors
last_verified_at: 2026-08-24
last_verified_git_sha: a33418fb8ae9cb9fd05832dd9bc1cb0778e08533
---

# Source of truth

## Hierarchy

When two records disagree, use this precedence order:

1. **Code, schemas, migrations, registries and machine-readable governance** at the canonical Git commit.
2. **Current GitHub documentation** explicitly verified against that commit.
3. **GitHub issues, pull requests and validation evidence**.
4. **Notion current project dashboard and decision records**.
5. **Architecture DOCX and other external design proposals**.
6. **Historical archives and implementation logs**.

A lower layer must be reconciled; it must not silently override a higher layer.

## GitHub boundary

GitHub owns:

- code and migrations;
- model/data/input/run schemas;
- model registry and versions;
- current technical documentation;
- issues, milestones, PRs and acceptance evidence;
- numerical and migration evidence;
- release gates and artifacts.

## Notion boundary

Notion owns:

- current project summary;
- strategic rationale;
- decision records;
- cross-project planning;
- human-readable status and risk summaries.

Notion must link to canonical GitHub issues and a verified Git SHA. It must not duplicate mutable file-level acceptance criteria or retain a superseded SHA as current state.

## DOCX boundary

The architecture DOCX is a **target architecture proposal**. Every revision must state:

```yaml
document_version: <version>
document_status: TARGET_ARCHITECTURE_PROPOSAL
verified_against_git_sha: <sha>
current_intended_use: E1_research_prototype
target_intended_use: E2_reproducible_research
```

Every major component/model must be marked `CURRENT`, `PARTIAL`, `TARGET`, `HYPOTHESIS` or `DEFERRED`.

## Conflict protocol

1. Freeze the disputed claim.
2. Identify the highest-precedence evidence.
3. Record the discrepancy in issue `#14` or the relevant issue.
4. Update lower-precedence records.
5. Re-run applicable documentation, semantic, migration or numerical checks.
6. Record the new verified SHA/date.

## Current canonical baseline

```text
repository = andreazedda/bmyCure4MM
branch = master
verified_head = a33418fb8ae9cb9fd05832dd9bc1cb0778e08533
```
