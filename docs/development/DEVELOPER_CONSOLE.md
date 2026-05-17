# Developer Console

The developer console is available to staff users at:

- `/research/developer/`

It runs offline checks and displays actionable rows with status, object ids, and next actions.

The page now follows the same explanation schema as the cockpit:

- Question answered
- Inputs used
- Method / computation
- Output meaning
- What you can do next
- Limitations

Raw details are collapsed by default under row-level developer details.

## Check Groups

- Data: orphan records, missing structured schedules, duplicate labs, unit completeness, key biomarker availability.
- Model: current twin state status, pre/post residuals, RMSE direction, underdetermined calibration, trajectory artifacts, predicted biomarkers, utility fields, schedule-collapse warnings.
- Causal: assumption-set presence and unsupported causal-effect flags.
- Scientific: model-reference entries from `twin_engine/model_references.json`, including placeholder status.
- Privacy: ignored local/private paths, staged sensitive paths, staged direct-identifier scan.
- Artifact: report and trajectory artifacts, provenance metadata, output hashes.

## Status Meaning

- `pass`: no issue detected by the offline check.
- `warn`: review required before sharing artifacts or using the output in a demo.
- `fail`: blocker before push, sharing, or research review.

Passing checks are audit signals only. They do not establish clinical validity, causal identification, or external model validation.

## Local Feedback

The feedback form writes JSONL records to `local_private/feedback/research_cockpit_feedback.jsonl`. That path is ignored by git and is intended for local developer notes only.

## Access Control

The route is login-required and staff-only. Non-staff users receive a permission error.