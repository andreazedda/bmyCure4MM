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

## Gli algoritmi dei farmaci sono testati?

Sì, nel senso pratico “il modello gira, produce KPI coerenti e i preset YAML cambiano davvero l’output”:

- La suite completa è eseguibile con `python manage.py test`.
- Caricamento/validazione preset PK/PD + schedule da YAML: `simulator/tests/test_pkpd_yaml.py`.
- Sanity sul modello (dose-response + persistenza output): `simulator/tests/test_simulation.py`.
- Sanity numerica (concentrazioni, AUC, KPI stabili): `simulator/tests/test_numerical_sanity.py`.
- Guardrail e warning sulle dosi: `simulator/tests/test_parameter_validation.py`.

Nota: questo non significa “validazione clinica” per ogni farmaco/regime; è test di robustezza e coerenza del modello.

## Simulare un “nuovo farmaco” digitalmente (limite attuale + workaround)

Ad oggi il Simulatore web (e `SimulationAttempt._resolve_solver_inputs`) usa **3 canali farmaco fissi**:

- `lenalidomide`
- `bortezomib`
- `daratumumab`

Quindi un vero “quarto farmaco” (es. `carfilzomib`) non è ancora un input di prima classe senza modifiche al codice.

### Workaround: usare uno slot come proxy (senza cambiare UI/codice)

Per fare un prototipo digitale e confrontarlo con le terapie attuali puoi “reinterpretare” uno slot (tipicamente `daratumumab`):

1. Scegli lo slot proxy (esempio: `daratumumab`).
2. Imposta a `0` le dosi degli slot che non vuoi usare.
3. Modifica il relativo preset YAML in `simulator/presets/drugs/`.

Importante: lascia `name` uguale alla key dello slot (es. `daratumumab`), altrimenti il simulatore non applica quel profilo.

Esempio (proxy su `daratumumab`):

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

### Come confrontarlo con le terapie attuali

Per un confronto riproducibile:

- Usa lo stesso Scenario.
- Imposta lo stesso `seed` e `cohort_size`.
- Esegui due run:
  - A) baseline (YAML originale, o `PREDLAB_V2=0` per usare i default legacy)
  - B) YAML modificato + `PREDLAB_V2=1`

Esempio minimo da CLI:

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

Confronta KPI tipo `tumor_reduction`, `healthy_loss`, `time_to_recurrence` e `auc` per farmaco.

### Se ti serve un vero 4° canale farmaco

Serve una piccola estensione al codice (nuovi field nel form + update di `SimulationAttempt._resolve_solver_inputs` + template + default/test). Se vuoi posso implementarla in modo minimale.
