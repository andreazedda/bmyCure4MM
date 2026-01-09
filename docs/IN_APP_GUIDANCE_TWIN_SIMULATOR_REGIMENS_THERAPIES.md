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

### 2b) Simulator: results that feel like a “clinical response”

**User expectation:** “Lancio la terapia simulata e voglio un responso: 1 mese, 3 mesi, tempo stimato a recidiva, informazioni significative… altrimenti non è una simulazione.”

**What we implemented:** the results panel now surfaces **readable, timepoint‑based outputs** instead of only a plot/CSV.

#### Output now shown in the Results card

- **Key timepoints** (~1 month/30d, ~2 months/60d, ~3 months/90d, end of horizon)
  - tumor reduction (proxy of response depth)
  - healthy loss (proxy of toxicity / immunosuppression)
- **Best response (nadir)**: day of minimum tumor burden + reduction at nadir
- **Durability index**: fraction of the horizon where tumor stays below baseline
- **Time to recurrence** (if reached within horizon): first time tumor rises above a threshold after nadir

#### “Average / mean” outputs via Virtual Cohort Size

The form already exposes **Virtual cohort size** (1/10/50/200). When set > 1, the simulator now runs an internal cohort (no extra DB rows) and displays **mean + P05–P95 bands** for:

- tumor reduction
- healthy loss
- durability
- time‑to‑recurrence stats (when recurrence happens in the cohort)
- milestone bands at ~1 month, ~3 months, end

#### Important safety note about “overall survival / life expectancy”

The current mathematical model does **not** compute a clinically validated overall survival (“tempo medio di vita”). We therefore show **time‑to‑event proxies** (time to recurrence + durability) and uncertainty bands, clearly labeled as model outputs.

#### Primary code/UX touchpoints (Results)

- `simulator/models.py`: richer `results_summary` (milestones, nadir, durability, cohort bands)
- `simulator/templates/simulator/_simulation_results.html`: renders timepoints + cohort tables
- `simulator/explain.py`: updated tooltip copy for durability

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

- Attribution / Attribuzione:
  - EN: Panel concept by Andrea Zedda, PhD Bioengineer — BMI (Bioingegneria Modulare Italiana) project.
  - IT: Idea del pannello di Andrea Zedda, PhD Bioengineer — progetto BMI (Bioingegneria Modulare Italiana).

- Add per‑drug badges next to regimen drug checkboxes (so mapping is visible without expanding details).
- Expand drug → parameter mapping as simulator safety rules grow.
- Add clickable “quick picks” for common adverse events in therapy entry.
