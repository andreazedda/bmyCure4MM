# Research Cockpit Test Strategy

## Django View Tests

- Staff and owners can render `/research/patient/<patient_id>/cockpit/`.
- Non-owners cannot access another patient's cockpit.
- Staff users can render `/research/developer/`; non-staff users receive 403.
- The guided initialize route renders on GET and avoids raw validation-error pages.

## Service Tests

- Assessment recommendation chooses the earliest highest-completeness usable assessment.
- Cockpit context uses structured labs before assessment fallback.
- Latest completed scenario per label is selected and sorted by heuristic utility.
- Schedule-collapse detection reports duplicate trajectory fingerprints.
- Developer checks return status/detail/object/next-action fields.

## Privacy And Safety Tests

- The direct-identifier scanner catches explicit patient names and common PHI markers.
- Local feedback writes to `local_private` only.
- The pre-push script is executable and contains the required Django check/test commands.

## Regression Tests

- The patient detail page includes `Patient Journey Tools` with explicit contrast classes.
- The patient detail page includes the `Twin Inputs (over time)` empty state when LDH, Beta-2M, and R-ISS are unavailable.
- Raw JSON is available only inside developer details blocks on research result pages.

## Content QA Tests

- The cockpit includes a workflow map and the explanation schema blocks: Question answered, Inputs used, Method / computation, Output meaning, What you can do next, Limitations, and Developer details.
- The cockpit includes definitions for RMSE, MAE, residual, research utility, toxicity constraint, mechanistic simulation, causal effect, provenance, and schedule collapse.
- The cockpit displays formulas for residuals, research utility, mechanistic simulation, and causal effect notation.
- The initialize page explains what initialization is and what it is not.
- The patient detail page explains source-data workspace versus research cockpit and shows a recommended next step.
- The developer console explains check groups and keeps raw details collapsed.
- The glossary route renders core terms.
- Forbidden clinical-claim phrases are absent from source templates and docs, excluding content-standard guidance that intentionally lists forbidden language.