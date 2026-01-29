# App `simulator`

## Scopo

Implementa:

- scenari clinici “training” (`Scenario`)
- esecuzione del modello matematico + reportistica (`SimulationAttempt.run_model`)
- regole/utility per “patient twin”
- preset PK/PD (`simulator/presets/`)
- API e UX auditing (varie `api_*.py`, `api_ux.py`)

## Modelli principali

Vedi:

- `simulator/models.py`: `Scenario`, `SimulationAttempt`
- `simulator/models_help.py`: `HelpArticle`

## Dati M2M

`Scenario.recommended_regimens` crea tabella di join:
- `simulator_scenario_recommended_regimens`

## Preset PK/PD

File YAML in `simulator/presets/drugs/` (es. `lenalidomide.yaml`).

!!! tip "Validazione"
    Ci sono test dedicati ai preset: vedi `simulator/tests/test_pkpd_yaml.py`.

