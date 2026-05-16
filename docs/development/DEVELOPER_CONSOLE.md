# Developer Console

The developer console is available to staff users at:

- `/research/developer/`

It runs offline checks and displays actionable rows with status, object ids, and next actions.

## Check Groups

- Data: orphan records, missing structured schedules, duplicate labs, unit completeness, key biomarker availability.
- Model: current twin state status, pre/post residuals, RMSE direction, underdetermined calibration, trajectory artifacts, predicted biomarkers, utility fields, schedule-collapse warnings.
- Causal: assumption-set presence and unsupported causal-effect flags.
- Scientific: model-reference entries from `twin_engine/model_references.json`, including placeholder status.
- Privacy: ignored local/private paths, staged sensitive paths, staged direct-identifier scan.
- Artifact: report and trajectory artifacts, provenance metadata, output hashes.

## Local Feedback

The feedback form writes JSONL records to `local_private/feedback/research_cockpit_feedback.jsonl`. That path is ignored by git and is intended for local developer notes only.

## Access Control

The route is login-required and staff-only. Non-staff users receive a permission error.