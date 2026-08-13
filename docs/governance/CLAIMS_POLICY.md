# Claims Policy

Status: **current and authoritative**
Version: `claims-policy-v1`

## Allowed claims

- “mechanistic research simulation”;
- “exploratory scenario comparison under explicit assumptions”;
- “model-relative diagnostic flag”;
- “data-completeness or model-readiness assessment”;
- “hypothesis generation or falsification result”;
- “reproducible research artifact.”

Every such claim must include an epistemic label and model/version context.

## Forbidden current claims

Runtime interfaces, APIs, reports, README content, and current documentation
must not state or imply:

- a recommended or best treatment;
- a patient-specific dose, schedule, hold, restart, or regimen instruction;
- predicted patient benefit;
- a clinically superior regimen;
- an identified causal treatment effect;
- a validated patient prognosis without a separately governed validation
  record.

A disclaimer does not make a prescriptive instruction acceptable. The
instruction itself must be absent.

## Evidence discipline

When biological or parameter provenance is missing, use `HEURISTIC`,
`UNKNOWN`, or `NEEDS_EVIDENCE`. `LITERATURE_INFORMED` means a traceable source
informed the content; it does not imply external validation or patient-level
applicability. `VALIDATED_EXTERNAL` may be used only when a governed external
validation record exists.

Historical material that documents superseded claims may remain only under
`docs/archive/` with an explicit non-current status notice.

## Invalidated pre-policy outputs

Chemistry efficacy scores, survival estimates, toxicity/risk-benefit values,
and patient-prognosis outputs generated before `claims-policy-v1` are
invalidated for scientific interpretation and comparison. They were not bound
to validated estimators or governed evidence records and must not be treated as
reproducible research evidence. This policy change does not alter the current
ODE equations or numerical model parameters; it removes unsupported output
claims and records their epistemic invalidation.
