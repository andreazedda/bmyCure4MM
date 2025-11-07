# Avvio rapido Predictive Lab

## Obiettivo e vincoli
- Permettere a non-esperti di simulare regimi MM in sicurezza.
- Rispettare i range dei preset e il limite di tossicità (Perdita sani ≤ 0.25).

## Passi guidati
1. **Seleziona lo scenario** – leggi profilo e R-ISS.
2. **Attiva il Gemello Paziente** – sincronizza biologia con LDH/β2M/FLC.
3. **Scegli il preset** – VRd, Dara-Rd, KRd con bounds protetti.
4. **Regola dosi/schedule** – resta nei range e pattern mostrati.
5. **Imposta orizzonte + coorte** – default 84gg / 1 paziente; aumenta a 10/50/200 per bande.
6. **Esegui la simulazione** – spinner + toast, poi controlla i KPI.
7. **Leggi i KPI** – riduzione tumorale ≥90% e perdita sani ≤25% = profilo forte.
8. **Apri l’Optimization Lab** – cerca Pareto con vincoli e ri-simula i migliori.
9. **Itera** – esporta CSV/grafici per audit.

## FAQ
- **AUC vs tossicità?** AUC alta = maggiore esposizione → più rischio di immunosoppressione.
- **Soglia perdita sani?** Guardrail a 0.25, avvisi gialli a 0.2, rossi a 0.3.
- **Quando aumentare la coorte?** 10 per sensibilità, 50+ per report con bande di incertezza.
