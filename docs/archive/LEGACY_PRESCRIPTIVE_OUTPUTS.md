# Archived legacy prescriptive-output design

> **ARCHIVED — NOT CURRENT BEHAVIOR.** This record documents a superseded
> pre-governance implementation. It must not be used for research
> interpretation, treatment selection, or clinical action.

Before `claims-policy-v1`, the repository contained a threshold-based UI called
“decision support” and documentation that translated simulation metrics into
clinical-sounding actions. The implementation was heuristic, did not establish
patient-specific validity, and did not identify causal effects.

Issue #12 removed those runtime instructions and replaced them with
model-relative diagnostic codes such as `constraint_flag`,
`simulated_low_impact_zone`, and `model_regrowth_signal`. Historical Git commits
remain the audit source for the exact superseded text and code.

Current authority: [Canonical Intended Use](../governance/INTENDED_USE.md).

Other superseded narratives removed from the current documentation tree in the
same governance change include the former efficacy-estimation, survival and
toxicity, form-enhancement, implementation-log, and learning-guide pages. Git
history preserves their exact content.
