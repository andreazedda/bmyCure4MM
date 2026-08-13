# Exploratory regimen catalog context

The legacy catalog in `simulator/regimen_suggester.py` contains literature
metadata and heuristic buckets. Current runtime views do not expose its raw
prescriptive rationale strings.

The governed HTML and JSON surfaces expose:

- the structured or derived inputs used by the bucket rules;
- `HEURISTIC` and `LITERATURE_INFORMED` labels;
- population-level evidence metadata when present;
- constraint flags without patient-specific actions;
- missing-data and uncertainty limitations.

Catalog membership does not establish clinical applicability, safety,
comparative benefit, or a causal treatment effect. Regimen names are references
for research comparison, not treatment selections.
