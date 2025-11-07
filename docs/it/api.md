# API

## Endpoint
| Metodo | Percorso | Descrizione |
| --- | --- | --- |
| GET | `/api/glossary/` | Glossario KPI bilingue. |
| GET | `/api/drugs/` | Profili PK/PD con range e schedule. |
| GET | `/api/simulations/<id>` | (pianificato) Riepilogo run + artifact. |
| POST | `/api/simulations/run` | (pianificato) Avvia simulazioni via payload. |
| GET | `/api/optim/pareto` | (pianificato) Dati Pareto. |

## Uso
- Autenticazione: sessione o token (future release).
- Consigliato rate-limit per script.
- Esempio: `curl -H "Accept: application/json" https://portal/api/glossary/`.
