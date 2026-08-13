# API Overview

## Endpoints
| Method | Path | Description |
| --- | --- | --- |
| GET | `/api/glossary/` | Governed bilingual model-KPI glossary. |
| GET | `/api/drugs/` | Governed hypothetical PK/PD model-input profiles. |
| GET | `/api/simulations/<id>` | (planned) Simulation summary + artifacts. |
| POST | `/api/simulations/run` | (planned) Launch run with payload. |
| GET | `/api/optim/pareto` | (planned) Stream Pareto front data. |

## Usage

- Intended use: `E1_research_prototype`.
- Scientific responses include `governance` with an epistemic label and model version.
- Model-input profiles are not clinical dose limits or treatment selections.
- Auth: session or token (future).
- Rate limit recommended for automation.
- Example:
```bash
curl -H "Accept: application/json" https://portal/api/glossary/
```
