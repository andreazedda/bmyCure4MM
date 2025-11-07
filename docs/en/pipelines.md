# Pipelines

## Layout
- `pipelines/processes/*.py` – orchestrators (binding visualizer, paths generator, drug evaluator).
- `pipelines/configs/*.yaml` – presets and general settings.
- `pipelines/data/paths.json` – resolved input/output directories.
- `pipelines/data/logs` / `outputs` – artifacts from batch runs.

## Flow
1. Generate paths (`paths_generator.py`) → writes JSON.
2. Load settings via `settings_generator.py` to combine general + task-specific YAML.
3. Execute process (e.g., `binding_visualizer.py`) – logs go to `pipelines/data/logs`.
4. Review outputs and sync with web UI through artifacts list.

## Tips
- Version-control YAML configs alongside code reviews.
- Keep `pipelines/data/outputs` small; archive heavy runs externally.
