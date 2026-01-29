# Glossary (terminologia)

Questa pagina raccoglie la terminologia più usata nella piattaforma, con definizioni **leggibili da umani** e riferimenti a dove si vede in UI/codice.

## Clinica (MM)

- **MM (Multiple Myeloma)**: neoplasia delle plasmacellule; nel simulatore è modellata come carico tumorale \(T(t)\).
- **CRAB**: *Calcium, Renal, Anemia, Bone* (manifestazioni tipiche).
- **R-ISS**: staging prognostico; nella piattaforma compare in `Assessment` e negli scenari.
- **FLC ratio**: rapporto catene leggere; usato come indicatore di attività di malattia (Twin/score).
- **LDH**: marker di aggressività/turnover; entra in risk/score.
- **β2-microglobulin**: marker prognostico; entra in R-ISS e rischio.

## Modelli matematici

- **PK (Pharmacokinetics)**: come \(C(t)\) cambia per input ed eliminazione.
- **PD (Pharmacodynamics)**: come \(C(t)\) si traduce in effetto \(E(C)\).
- **Emax**: modello saturante dell’effetto; vedi `Guides → Mathematical Models`.
- **EC50**: concentrazione che produce 50% dell’effetto massimo.
- **K (carrying capacity)**: “limite superiore” di crescita in modelli logistici.
- **Half-life**: tempo di dimezzamento; determina \(k_{elim}\).
- **AUC**: area sotto la curva di concentrazione; proxy di esposizione complessiva.
- **KPI**: metriche riassuntive (es. `tumor_reduction`, `healthy_loss`).

## Ottimizzazione

- **Multi-objective optimization**: ottimizzare più obiettivi contemporaneamente (efficacia/sicurezza/esposizione).
- **Pareto front**: insieme di soluzioni non dominate (nessuna è “meglio in tutto” delle altre).
- **Trial**: una proposta di parametri (dosi/horizon) valutata dal modello.
- **Constraint**: vincolo “hard” (es. `healthy_loss <= 0.25`).

## Software / dati

- **Django model**: classi Python ORM (es. `clinic.Patient`) che mappano tabelle DB.
- **Migration**: versione dello schema DB (Django migrations).
- **Artifact**: file generato (CSV/HTML/JSON) salvato sotto `MEDIA_ROOT`.
- **JSONField**: campo strutturato; in SQLite è serializzato in `TEXT` con check JSON.

## Dove trovare “la cosa giusta”

- Oggetti DB: `Reference → Database Objects`
- Parametri simulatore: `Reference → Simulator Parameters`
- Modelli e KPI: `Guides → Mathematical Models`
- Ottimizzazione: `Guides → Optimization Theory`

