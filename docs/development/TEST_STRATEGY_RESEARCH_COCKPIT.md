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