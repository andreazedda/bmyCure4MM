# API Overview

## Endpoints
| Method | Path | Description |
| --- | --- | --- |
| GET | `/api/glossary/` | Bilingual KPI glossary for frontends. |
| GET | `/api/drugs/` | PK/PD preset profiles with ranges + schedules. |
| GET | `/api/simulations/<id>` | (planned) Simulation summary + artifacts. |
| POST | `/api/simulations/run` | (planned) Launch run with payload. |
| GET | `/api/optim/pareto` | (planned) Stream Pareto front data. |

## Usage
- Auth: session or token (future).
- Rate limit recommended for automation.
- Example:
```bash
curl -H "Accept: application/json" https://portal/api/glossary/
```
