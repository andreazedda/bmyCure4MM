---
title: Documentation Inventory
status: CANONICAL
owner: Andrea Zedda
audience: maintainers and documentation reviewers
last_verified_at: 2026-08-24
last_verified_git_sha: bf097810b337dc6b766cda04497005670cd96513
machine_readable_source: ../DOCUMENTATION_INVENTORY.yaml
---

# Documentation inventory

The machine-readable inventory is [`docs/DOCUMENTATION_INVENTORY.yaml`](../DOCUMENTATION_INVENTORY.yaml).

It classifies documentation by explicit path and ordered glob rules. The first matching rule wins. Any unmatched file is:

```text
UNCLASSIFIED_REQUIRES_REVIEW
```

and must remain outside current navigation until reviewed.

## Inventory classifications

| Classification | Meaning |
|---|---|
| `CANONICAL` | Current authority for project state, governance, architecture or contracts |
| `CURRENT_SUPPORTING` | Current support material governed by a canonical parent |
| `TARGET_SUPPORTING` | Approved target design, explicitly not current implementation |
| `HISTORICAL_ARCHIVE` | Retained history; not a current product claim |
| `DUPLICATE` | Redundant with a canonical document |
| `DELETE_CANDIDATE` | No continuing value after review |
| `UNCLASSIFIED_REQUIRES_REVIEW` | Currency and ownership not yet proven |

## Review rule

A document may be moved into current navigation only when:

1. its status and owner are explicit;
2. its claims match the canonical intended use;
3. current/target/hypothesis content is separated;
4. referenced paths, routes, commands, formulas and model versions are verified;
5. its `last_verified_git_sha` is recorded.
