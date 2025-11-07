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
