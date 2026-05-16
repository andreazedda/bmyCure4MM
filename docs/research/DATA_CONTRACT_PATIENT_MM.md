# Data Contract: Patient MM Research Workflow

## Scope

This contract describes the minimum structured data expected by the real-patient research workflow.

The workflow is for research and modeling only.

## Identifier Policy

- demographics may exist in `clinic.Patient` for the operational app
- research artifacts must use pseudonymized references only
- do not include `first_name`, `last_name`, `mrn`, or free-text notes in research exports by default

## Required Core Fields

### Demographics pseudonymized

- patient internal id
- owner or staff access context
- sex
- birth date stored operationally, but not exported by default

### Diagnosis date

- `clinic.Patient.diagnosis_date`

### Staging

- `clinic.Assessment.r_iss`
- optional supporting labs: beta-2 microglobulin, albumin, LDH

### Cytogenetics

- `clinic.PatientCytogenetics.abnormality`
- `clinic.PatientCytogenetics.detected_on`
- `clinic.PatientCytogenetics.method`

### Longitudinal labs

At least one of the following per research-relevant assessment:

- `m_protein_g_dl`
- `flc_ratio`
- `hemoglobin_g_dl`

Additional supported fields:

- `kappa_mg_l`
- `lambda_mg_l`
- `calcium_mg_dl`
- `creatinine_mg_dl`
- `beta2m_mg_l`
- `albumin_g_dl`
- `ldH_u_l`
- `response`

### Therapy timeline

- `clinic.PatientTherapy.start_date`
- `clinic.PatientTherapy.end_date`
- `clinic.PatientTherapy.regimen`
- `clinic.PatientTherapy.outcome`
- `clinic.PatientTherapy.adverse_events`

### Structured dose/schedule

Required for transparent schedule conversion:

- `clinic.PatientTherapy.doses`
- `clinic.PatientTherapy.cycle_length_days`
- `clinic.PatientTherapy.days_on`

Optional but recommended:

- `clinic.PatientTherapy.schedule_notes`

### Response

- `clinic.Assessment.response`
- longitudinal deltas across assessments

### Toxicity

- `clinic.PatientTherapy.adverse_events`
- `clinic.models_symptoms.SymptomAssessment` neuropathy, infection, ECOG, and CRAB flags

### Symptoms

- bone pain
- fatigue
- neuropathy
- infection status
- appetite and weight loss fields when present

### Imaging, genomics, MRD placeholders

Not yet modeled as first-class tables in the research layer.

Expected placeholder categories:

- imaging summary
- genomics summary
- MRD summary
- source and collection date

### Source quality

Each `clinic.PatientTherapy` can now store:

- `unknown`
- `inferred`
- `patient_reported`
- `clinical_record`
- `curated_research`

### Provenance

Each structured therapy row can now store `provenance` JSON.

Research runs additionally persist:

- code commit hash
- config hash
- preset hashes
- input hash
- output hash
- solver parameters
- model version

## Minimal Research Readiness Rules

A patient is minimally ready when:

- at least one assessment exists
- the latest assessment has a date and at least one of `m_protein_g_dl`, `flc_ratio`, or `hemoglobin_g_dl`
- a research twin state can be initialized from an assessment

Calibration-ready usually means:

- at least two dated assessments
- a coherent time interval

Transparent therapy schedule conversion requires:

- structured dose entries in `doses`
- cycle length
- explicit `days_on`

If those fields are missing, the system must report missingness explicitly instead of inventing a schedule.

## Export Rules

Allowed by default:

- pseudonymized patient reference
- research twin state id
- summary metrics
- residual metrics
- provenance hashes

Blocked by default:

- names
- MRN
- free-text notes
- birth date
- schedule notes
