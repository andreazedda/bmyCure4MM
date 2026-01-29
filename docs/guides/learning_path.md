# Percorso guidato (trovare info in 30 secondi)

Questa pagina è il “**GPS**” della documentazione: cosa leggere, in che ordine, e dove trovare rapidamente le informazioni.

## Se cerchi… vai qui

| Voglio… | Vai a… |
| --- | --- |
| capire cosa fa il sistema | `Guides → Architecture` |
| vedere **tabelle DB** e relazioni | `Guides → Database` |
| trovare **un oggetto DB** (Model ↔ tabella ↔ file) | `Reference → Database Objects` |
| capire il **modello ODE / PK/PD** e KPI | `Guides → Mathematical Models` |
| capire l’**Optimization Lab** (Pareto, obiettivi) | `Guides → Optimization Theory` |
| vedere regole e soglie “auditabili” | `Guides → Decision Algorithm` |
| avviare local/dev/docker | `Guides → Operations` + `Guides → Development` |
| consultare documenti lunghi/tecnici | `Deep Dives` |

## Percorso “Clinico” (30–45 min)

1. `IT → Simulatore` (workflow + KPI)
2. `Guides → Mathematical Models` (solo sezioni “KPI” e “Patient Twin”)
3. `Guides → Decision Algorithm` (soglie + regole principali)
4. `IT → Optimization Lab` (come leggere la Pareto table)

Se vuoi un percorso “hands-on” con dati demo: `Deep Dives → Learning (IT) → Guida apprendimento`.

## Percorso “Developer” (45–90 min)

1. `Guides → Architecture`
2. `Reference → Configuration`
3. `Reference → Database Objects`
4. `Guides → Database` (DDL e relazioni)
5. `Guides → Mathematical Models` (sezioni “pipeline/solver/outputs”)

## Percorso “Ricerca/Modeling” (60–120 min)

1. `Guides → Mathematical Models` (tutto)
2. `Guides → Optimization Theory`
3. `Deep Dives → Mathematical Models` (documento lungo già presente nel repo)

## Come cercare bene (MkDocs)

La search di MkDocs Material è la via più veloce:

- cerca per `PatientTherapy`, `clinic_patienttherapy`, `SimulationAttempt`, `simulator_simulationattempt`
- cerca per KPI: `tumor_reduction`, `healthy_loss`, `AUC`, `time_to_recurrence`

!!! tip "Se vedi ‘Syntax error in text’"
    È un errore di parsing Mermaid in una pagina. Vai su `Home` e verifica che *diagramma* e *formula* si vedano; se si vedono lì ma non altrove, il problema è localizzato a un blocco Mermaid specifico.
