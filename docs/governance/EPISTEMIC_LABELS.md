# Epistemic Labels

Status: **current and authoritative**

| Label | Meaning | Does not mean |
|---|---|---|
| `OBSERVED` | Direct structured observation retained with provenance. | Complete, correct, or causal. |
| `DERIVED` | Deterministic transformation of declared inputs. | Independently observed. |
| `USER_PROVIDED` | Value entered by a user and not independently verified. | Source-verified. |
| `SIMULATED` | Output of a versioned computational model. | Patient outcome or causal effect. |
| `HEURISTIC` | Rule or threshold used for research orientation. | Clinically validated rule. |
| `LITERATURE_INFORMED` | Traceable literature informed the statement or parameter. | Externally validated in this implementation. |
| `HYPOTHETICAL` | Explicit assumption or scenario introduced for exploration. | Observed fact. |
| `VALIDATED_EXTERNAL` | Passed a separately governed external validation protocol. | Universal validity or clinical authority. |

Labels describe epistemic status, not quality scores. Multiple labels may be
needed when an output combines observed inputs, derived features, and simulated
results. Missing provenance must fail closed to `HEURISTIC`, `UNKNOWN`, or
`NEEDS_EVIDENCE` rather than being promoted to literature support.
