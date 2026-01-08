# In‑App Guidance Improvements (Twin + Simulator + Regimens + Therapies)

Date: 2026‑01‑08

## Goal

Make the platform *self‑explanatory inside the UI* (no Admin redirects, no external docs needed) by:

- Clarifying what really affects the Patient Twin (structured assessments/labs only).
- Replacing “type numbers” with guided choices where possible.
- Making regimen and therapy workflows feel guided and non‑ambiguous.
- Keeping guidance aligned with actual code behavior (no “fake” UX).

## Key UX principles adopted

- **Transparency:** clearly show what inputs drive the model.
- **Guided entry:** prefer dropdowns, presets, and step‑by‑step hints over free text.
- **No magic notes:** free‑text notes are for humans; they do not change the Twin.
- **Follow the data:** Twin changes when new *assessments* (labs) are added.

## What was implemented

### 1) Twin: “what affects what” made explicit

**User problem:** “Non capisco quali scelte influenzano quali parametri del digital twin” and confusion about notes.

**Solution:** in‑app explanations that map:

- `Assessment` labs (R‑ISS, LDH, β2M, FLC ratio)
- ⟶ derived `risk_score`
- ⟶ biology parameters used by the simulator (growth rates, carrying capacities, immune index)

**Also made explicit:** notes do not auto‑change the Twin; only structured labs do.

#### Primary code/UX touchpoints

- `simulator/twin.py`: Twin building logic (`build_patient_twin(assessment)`)
- `simulator/templates/simulator/interactive_tutorial.html`: mapping table + explanation
- `clinic/templates/clinic/patient_detail.html`: short tips near clinical notes / assessments

---

### 2) Simulator: guided choices (no numeric typing)

**User problem:** “Perché non metti la possibilità di fare delle scelte guidate?”

**Solution:** add guided dropdowns that adjust biology parameters without requiring manual numbers.

- Guided tumor aggressiveness (lower / typical / higher)
- Guided immune status (better / typical / worse)

These apply multiplicative adjustments (with clamping) to the resolved parameter payload after Twin auto‑merge.

#### Primary code touchpoints (Simulator)

- `simulator/forms.py`: added guided fields
- `simulator/templates/simulator/_simulation_form.html`: rendered guided controls
- `simulator/models.py`: `_apply_guided_twin_choices(params)` applied in `SimulationAttempt.run_model()`

---

### 3) Regimens: strongly guided creation (presets + wizard + drug builder)

**User problem:** regimen creation felt like guesswork (“troppo poco guidato”).

**Solution:** in‑page guidance on the regimen form:

- 1‑click preset templates (fill common regimens)
- Step‑by‑step wizard (line + family → autofill Name/Line/Components)
- Drug checkbox builder that fills the components string
- Server‑side normalization of components (trim, lowercase, de‑duplicate)

**Plus:** added an in‑page section explaining **drug → health parameter mapping** used by simulator safety rules, and where those parameters live in the UI.

#### Primary code touchpoints (Regimens)

- `clinic/forms.py`: `RegimenForm.clean_components()` normalization and improved placeholders
- `clinic/templates/clinic/regimen_form.html`: presets + wizard + checkbox builder + mapping section
- `clinic/tests/test_patient_crud.py`: regression test for components normalization

---

### 4) Therapies: make entry guided (patient detail page)

**User problem:** “anche per la sezione therapies non è abbastanza guidata” — unclear what to do.

**Solution:** guide the “Add Therapy” block in the patient detail page:

- “Quick guide (4 steps)” at the top of the form
- Extra hints under date fields (start date meaning; leave end date blank if ongoing)
- Outcome suggestions via `datalist` (CR/VGPR/PR/SD/PD…)
- Better examples for adverse events
- More informative empty-state and a button that opens “Add Therapy” directly

Also added a small note clarifying how **Observed effect** is computed (needs an assessment before and after therapy start; otherwise “–”).

#### Primary code touchpoints (Therapies)

- `clinic/forms.py`: `PatientTherapyForm` outcome datalist + clearer placeholders
- `clinic/templates/clinic/patient_detail.html`: quick guide + outcome datalist + empty-state improvements

## Reality check: what affects the Twin

- **Therapy entries**: record regimen + dates + outcome/AEs for tracking; do **not** directly change the Twin.
- **Assessments (labs)**: are the inputs that update the Twin (and therefore affect simulator parameters).

## Tests executed

- `python manage.py test clinic.tests.test_patient_crud -v 2` (29 tests) ✅
- `python manage.py test clinic.tests.test_patient_crud.RegimenInAppCRUDTests -v 2` ✅

## Files changed (high level)

- `clinic/forms.py`
- `clinic/templates/clinic/patient_detail.html`
- `clinic/templates/clinic/regimen_form.html`
- `clinic/tests/test_patient_crud.py`
- `simulator/forms.py`
- `simulator/models.py`
- `simulator/templates/simulator/_simulation_form.html`
- `simulator/templates/simulator/interactive_tutorial.html`
- `simulator/twin.py`

## Notes / next optional improvements

These are optional and were discussed as potential follow‑ups:

- Add per‑drug badges next to regimen drug checkboxes (so mapping is visible without expanding details).
- Expand drug → parameter mapping as simulator safety rules grow.
- Add clickable “quick picks” for common adverse events in therapy entry.
