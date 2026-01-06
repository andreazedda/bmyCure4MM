# Tutorial: First Complete Simulation (Zero to Result)

**Goal:** Simulate the effect of a treatment regimen on a Multiple Myeloma patient, see tumor curves and patient impact.

---

## The Complete Flow (5 Steps)

```
┌─────────────────┐
│ 1. Create Patient│  ← Basic demographic record
└────────┬────────┘
         │
┌────────▼──────────┐
│ 2. Add Labs       │  ← Lab values (R-ISS, LDH, β2M, FLC)
└────────┬──────────┘     This is the "snapshot" for the Twin
         │
┌────────▼────────────┐
│ 3. Go to Simulator  │  ← Choose a scenario/regimen
└────────┬────────────┘
         │
┌────────▼──────────────────┐
│ 4. Configure & Enable Twin│  ← Twin uses labs to compute biological params
└────────┬──────────────────┘
         │
┌────────▼──────┐
│ 5. Run & Review│  ← Results: curves, toxicity, derived parameters
└─────────────────┘
```

---

## STEP 1: Create a Patient

1. **Navigate to:** `http://127.0.0.1:8000/patients/` (or click **Clinic** → **Patients**)
2. **Click:** green **"+ New Patient"** button at the top
3. **Fill minimal fields:**
   - First Name: `John`
   - Last Name: `Doe`
   - Date of Birth: `1960-01-15`
   - Sex: `M`
4. **Save**

✅ Patient created. Now we need a lab "snapshot" for this patient.

---

## STEP 2: Add an Assessment (Lab Snapshot)

An **Assessment** = a lab/clinical snapshot **at a specific date**.

1. **After saving the patient**, you're on the patient detail page
2. **Find the "Assessments" section** (below demographics)
3. **Click:** button **"Add Assessment"**
4. **Fill lab fields (example high-risk patient):**
   - **Date:** today (e.g., `2026-01-05`)
   - **R-ISS Stage:** `III` (high risk)
   - **LDH (U/L):** `600` (very high, normal is <280)
   - **Beta-2 Microglobulin (mg/L):** `12.0` (high, normal <3.5)
   - **FLC Ratio (κ/λ):** `8.5` (very imbalanced, normal ~1.0)
   - **M-protein (g/dL):** `4.5` (optional)
   - **Response:** leave empty or `PD` (progressive disease)
5. **Save**

✅ Lab snapshot created. Now the **Patient Twin** can compute biological parameters from these values.

---

## STEP 3: Go to the Simulator

1. **Go back to home** (click logo at top or `/`)
2. **Click:** **"Simulator"** section → **"Browse Scenarios"**
3. **Choose a scenario** (e.g., scenario ID `2` - typically a triplet regimen VRd: Bortezomib + Lenalidomide + Dexamethasone)
4. **Click on the scenario title** to open the detail page

You'll now see the simulation form panel.

---

## STEP 4: Configure Simulation & Enable Twin

### 4.1 Choose a Preset (optional but recommended)

- In the **"Preset"** field, select `standard_triplet` or similar
- This auto-populates "safe" doses

### 4.2 Enable Patient Twin

1. **Scroll down** to section **"4. Advanced & Twin"**
2. **Toggle the switch** "Use Patient Twin"
3. **Select the Assessment:**
   - In the **"Twin Assessment"** dropdown, you'll see the assessment you just created
   - Select the one with today's date (e.g., `John Doe - 2026-01-05`)
4. **Choose mode:**
   - **Twin Biology Mode:** select **"Auto"**
   - This makes the Twin compute `tumor_growth_rate`, `carrying_capacity_tumor`, etc. automatically from labs

### 4.3 Set Horizon and Cohort (optional)

- **Time Horizon (days):** `180` (6 months simulation)
- **Cohort Size:** `1` (keep low for speed; later try 10 or 50 to see uncertainty bands)

### 4.4 (Optional) Check Doses

If you chose a preset, doses are pre-filled. Otherwise:
- **Lenalidomide:** 25 mg
- **Bortezomib:** 1.3 mg/m²
- **Daratumumab:** 16 mg/kg

---

## STEP 5: Run Simulation and View Results

1. **Click:** blue **"Run Simulation"** button at bottom
2. **Wait** 10-30 seconds (you'll see a loading indicator)
3. **Results appear below the form:**

### What You'll See:

#### A) **Time Plot (Interactive Chart)**
   - **X-axis:** days (0 → 180)
   - **Y-axis:** cells
   - **Red curve:** tumor cells (should drop!)
   - **Green curve:** healthy cells (try to keep high)
   - If cohort size > 1, you'll see variability bands

#### B) **KPIs (Key Performance Indicators)**
   - **Tumor Reduction (%):** e.g., `92.5%` → excellent response!
   - **Healthy Loss (%):** e.g., `15.2%` → manageable toxicity
   - **Time to Recurrence (days):** e.g., `210` → time before tumor regrows
   - **AUC (Area Under Curve):** total drug exposure

#### C) **Downloadable Files (at bottom of results)**
   - **`twin_params.json`** ← Parameters derived by Twin (risk_score, tumor_growth_rate, etc.)
   - **`simulation_data.csv`** ← Raw day-by-day data
   - **`plot.html`** ← Standalone interactive chart

#### D) **What to Look For in `twin_params.json`**

Download and open with a text editor. You'll see something like:

```json
{
  "risk_score": 0.82,
  "tumor_growth_rate": 0.0356,
  "healthy_growth_rate": 0.0136,
  "carrying_capacity_tumor": 2100000000.0,
  "carrying_capacity_healthy": 280000000000.0,
  "immune_compromise_index": 1.215,
  "assessment_id": 3
}
```

**Interpretation:**
- `risk_score: 0.82` → high-risk patient (from R-ISS=III, high LDH, high β2M)
- `tumor_growth_rate: 0.0356` → tumor grows fast (near max 0.04)
- `healthy_growth_rate: 0.0136` → slow healthy cell recovery (near min)
- These values were **automatically computed** from the labs you entered in the Assessment

---

## Quick FAQs

### Q: The "Twin Assessment" dropdown is empty, why?

**A:** Two possible reasons:
1. **You didn't create an Assessment** → go back to STEP 2
2. **Permissions:** if you're not staff/editor, you only see Assessments for patients with `owner = your user`. Quick fix for testing: log in as admin/staff.

### Q: What happens if I DON'T enable the Twin?

**A:** Simulation uses default biological values (fixed, not personalized to patient labs). Still works, but not "patient-specific".

### Q: How do I see the effect of different doses?

**A:**
1. Modify dose fields (e.g., reduce Lenalidomide from 25 to 15 mg)
2. Run simulation again
3. Compare KPIs: lower tumor reduction? lower healthy loss? Efficacy/toxicity trade-off.

### Q: Can I simulate multiple patients at once?

**A:** Increase **Cohort Size** (e.g., 50). You'll see variability bands in the chart. But each run uses **one same** patient (same Twin/Assessment): variability represents stochastic model uncertainty, not different patients.

### Q: Where's the scientific explanation of the Twin?

**A:** Detailed file: `docs/development/PATIENT_TWIN_AS_BUILT.md`

---

## Next Steps (After First Test)

1. **Try patients with different risk:**
   - Create a second Assessment with R-ISS=I, LDH=200, β2M=2.5 → you'll see much lower `risk_score` and different biological parameters
2. **Compare different regimens:**
   - Scenario A: triplet VRd
   - Scenario B: doublet Rd (without Bortezomib)
   - Which gives better tumor reduction? Which has less toxicity?
3. **Optimize doses:** use the "Optimize Regimen" button (if you're an editor) to let the system find optimal doses

---

## Ultra-Quick Glossary

| Term | What It Means |
|------|---------------|
| **Assessment** | Lab/clinical snapshot (R-ISS, LDH, etc.) at a specific date |
| **Patient Twin** | Computation that transforms Assessment into biological parameters (growth rate, etc.) |
| **Scenario** | A "clinical case" with an associated drug regimen |
| **Preset** | Pre-configured dose set for a standard regimen |
| **Cohort Size** | Number of simulation "replicas" to estimate variability |
| **Time Horizon** | Simulation duration in days |
| **Auto (Twin Mode)** | Twin **overwrites** base biological parameters with lab-derived ones |
| **Manual (Twin Mode)** | Twin **doesn't** overwrite; you use manually entered values |

---

## Common Troubleshooting

### Error: "You have 1 unapplied migration(s)"

**Fix:**
```bash
python manage.py migrate
```

### Error: "no such column: clinic_patient.owner_id"

**Fix:** you merged without migrating. Run:
```bash
python manage.py migrate
```

### Form not showing / broken layout

**Fix:** check Django is serving static files:
```bash
python manage.py collectstatic --noinput
```

### No Assessments in dropdown

**Possible causes:**
1. You didn't create one → see STEP 2
2. Patient has `owner != your user` and you're not staff → log in as admin or assign owner
3. Refresh the page (Cmd+R or Ctrl+R)

---

**End of Tutorial!** 🎉

Questions? Check `PATIENT_TWIN_AS_BUILT.md` for technical details or open a GitHub issue.
