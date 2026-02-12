# PK/PD Presets

## YAML Schema
```yaml
name: lenalidomide
pk:
  half_life: 9.0   # hours
  Vd: 60.0         # L
pd:
  Emax: 0.85
  EC50: 8.0
unit: mg
schema: 21on/7off
schedule:
  type: cycle
  cycle_length: 28
  days_on: 21
  administration_window_days: 1
```
- `type`: cycle | weekly | interval | pulsed | continuous.
- `dose_range.min/max` enforce form guardrails.

## Tutorial
1. Tweak `days_on` to 14 in `lenalidomide.yaml`.
2. Reload form (select preset) – range note updates.
3. Compare AUC in results vs previous 21-day on schedule; expect lower exposure.

## Are drug algorithms tested?

Yes, at a practical “does it run + does it change outputs” level:

- The full suite is runnable via `python manage.py test`.
- PK/PD YAML loading + schedule building is covered by `simulator/tests/test_pkpd_yaml.py`.
- Core model behavior is covered by `simulator/tests/test_simulation.py` (dose-response sanity) and `simulator/tests/test_numerical_sanity.py` (outputs contain concentrations, AUC, and stable KPIs).
- Dose guardrails/warnings are covered by `simulator/tests/test_parameter_validation.py`.

What’s *not* claimed: clinical validity for every drug/regimen listed elsewhere in the UI.

## Simulate a “new drug” digitally (current limitation + workaround)

Today the web Simulator UI (and `SimulationAttempt._resolve_solver_inputs`) uses **three fixed drug channels**:

- `lenalidomide`
- `bortezomib`
- `daratumumab`

So a truly new drug name (e.g. `carfilzomib`) is not yet a first-class simulation input without code changes.

### Workaround: use a proxy slot (no UI/code changes)

If you want to prototype a new agent digitally and compare it to current therapies, you can **reuse one slot** (commonly `daratumumab`) as a proxy:

1. Pick the proxy channel you will “reinterpret” (example: `daratumumab`).
2. Set the other channels you don’t want to use to `0` in the simulator form.
3. Edit the corresponding YAML preset under `simulator/presets/drugs/`.

Important: keep `name` equal to the proxy key (e.g. `daratumumab`), otherwise the simulator won’t apply the profile.

Example (proxying `daratumumab` as a new candidate):

```yaml
name: daratumumab
display_name: "NEW_CANDIDATE (proxy slot)"
unit: mg/kg
dose_range:
  min: 1
  max: 20
pk:
  half_life: 120.0
  Vd: 40.0
pd:
  Emax: 0.70
  EC50: 15.0
schedule:
  type: interval
  interval_days: 14
  administration_window_days: 1
```

### How to compare against existing therapies

To compare “baseline therapy” vs “new digital drug” in a reproducible way:

- Use the **same Scenario**.
- Use the same `seed` and `cohort_size`.
- Run two attempts:
  - A) baseline YAML (or `PREDLAB_V2=0` to force legacy defaults)
  - B) modified YAML + `PREDLAB_V2=1`

Minimal CLI example:

```bash
PREDLAB_V2=1 python manage.py shell -c "\
from django.contrib.auth import get_user_model;\
from simulator import models;\
User=get_user_model(); u=User.objects.get(username='admin');\
s=models.Scenario.objects.first();\
p={'baseline_tumor_cells':1e9,'baseline_healthy_cells':5e11,\
   'lenalidomide_dose':25.0,'bortezomib_dose':1.3,'daratumumab_dose':16.0,\
   'time_horizon':84.0,'tumor_growth_rate':0.023,'healthy_growth_rate':0.015,\
   'interaction_strength':0.05};\
a=models.SimulationAttempt.objects.create(scenario=s,user=u,parameters=p,cohort_size=30,seed=2026);\
print(a.run_model()['auc'])"
```

Compare KPIs like `tumor_reduction`, `healthy_loss`, `time_to_recurrence`, and per-drug `auc`.

### If you need a true 4th drug channel

That requires code changes (new form field(s) + updating `SimulationAttempt._resolve_solver_inputs` + templates + defaults/tests). If you want, I can implement this next in a minimal way.
