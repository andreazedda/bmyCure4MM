# Gemello Paziente

## Input
- Campi Assessment: R-ISS, LDH, β2M, rapporto FLC.

## Mapping
- Punteggio di rischio pesato → `_lerp` per crescita e capacità di carico.
- YAML `twin_risk.yaml` contiene pesi e range modificabili.

## Tutorial
1. Assessment basso rischio: R-ISS I, LDH normale.
2. Assessment alto rischio: R-ISS III, LDH 600, β2M 12.
3. Esegui simulazioni con `use_twin=True` e confronta KPI.
4. Apri `twin_params.json` per verificare i parametri derivati.
