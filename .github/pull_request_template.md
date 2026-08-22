## Linked issue

Closes #

## Scope and source of truth

Describe the bounded change, its source-of-truth documents or contracts, and
explicit non-goals.

## Classification and risk

Change classification (select all that apply):

- [ ] Code
- [ ] Scientific model or equation
- [ ] Database schema or migration
- [ ] Dataset or importer
- [ ] Evidence, claim, or epistemic label
- [ ] Security, privacy, or authorization
- [ ] Documentation or operations

Overall risk: `low / medium / high / critical`

| Area | Impact, or `none` with justification |
| --- | --- |
| Model/version | |
| Schema/migrations | |
| Data/provenance | |
| Evidence/claims | |
| Intended use | |
| Authorization | |

## Scientific and numerical impact

State whether equations, coefficients, solver behavior, tolerances, parameter
defaults, units, observation semantics, comparability, or numerical outputs
change. For an approved change, attach the model-version decision and a
privacy-safe output diff.

- [ ] No numerical output changed, or model-version/output-diff evidence is attached.
- [ ] No intended-use expansion occurred.

## Migrations

Migration status: `none / required and included / deferred (explain)`

- [ ] Migration status is declared and `makemigrations --check --dry-run` is clean.

## Privacy, PHI, and security

Describe the privacy, PHI, secret-handling, artifact, and authorization impact.
Do not attach patient data, source clinical documents, credentials, databases,
logs, or private payloads.

- [ ] No private data or direct identifiers were added.
- [ ] The research safety gate was executed on the candidate diff.
- [ ] Screenshots, if any, contain synthetic data only.

## Tests and evidence

List commands, exact outcomes, CI run links, numerical-baseline comparison, and
any privacy-safe failure evidence. Passing software checks do not establish
model validation, external validation, or clinical utility.

## Rollback

Describe the reversible rollback path, including dependency, migration, model,
dataset, and artifact implications where applicable.

## Known limitations

List residual risks, deferred work, and conclusions this change does not support.
