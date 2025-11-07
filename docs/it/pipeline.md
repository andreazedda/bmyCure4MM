# Pipeline

## Componenti
- `pipelines/processes` – script batch.
- `pipelines/configs` – YAML di configurazione.
- `pipelines/data` – percorsi, log, output.

## Flusso
1. Genera percorsi (`paths_generator.py`).
2. Combina setting (`settings_generator.py`).
3. Lancia il processo (es. `binding_visualizer.py`).
4. Raccogli log/artifact per l’interfaccia web.

## Suggerimenti
- Versiona i preset.
- Pulisci regolarmente `outputs/` per evitare file voluminosi.
