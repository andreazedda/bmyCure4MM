# Preset PK/PD

## Schema YAML
```yaml
name: lenalidomide
pk:
  half_life: 9.0
  Vd: 60.0
pd:
  Emax: 0.85
  EC50: 8.0
unit: mg
schedule:
  type: cycle
  cycle_length: 28
  days_on: 21
```
- `type`: cycle | weekly | interval | pulsed | continuous.
- `dose_range` definisce i limiti mostrati nel form.

## Tutorial
1. Modifica `days_on` da 21 a 14.
2. Aggiorna il preset nel form: il range mostrato cambia.
3. Riesegui la simulazione e confronta l'AUC con il setup originale.
