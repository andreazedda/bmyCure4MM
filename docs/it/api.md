# API

## Endpoint
| Metodo | Percorso | Descrizione |
| --- | --- | --- |
| GET | `/api/glossary/` | Glossario bilingue governato dei KPI del modello. |
| GET | `/api/drugs/` | Profili ipotetici e governati degli input PK/PD. |
| GET | `/api/simulations/<id>` | (pianificato) Riepilogo run + artifact. |
| POST | `/api/simulations/run` | (pianificato) Avvia simulazioni via payload. |
| GET | `/api/optim/pareto` | (pianificato) Dati Pareto. |

## Uso

- Uso previsto: `E1_research_prototype`.
- Le risposte scientifiche includono `governance`, etichetta epistemica e versione del modello.
- I profili di input non sono limiti clinici di dose né selezioni terapeutiche.
- Autenticazione: sessione o token (future release).
- Consigliato rate-limit per script.
- Esempio: `curl -H "Accept: application/json" https://portal/api/glossary/`.
