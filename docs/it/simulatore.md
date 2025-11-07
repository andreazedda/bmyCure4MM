# Guida al Simulatore

## Panoramica
- Modello PK/PD + crescita logistica risolto via `MathematicalModel`.
- Form HTMX in quattro gruppi con suggerimenti bilingue.

## Flusso
1. Controlla laboratori (ANC, piastrine, gravidanza) – blocco duro se fuori soglia.
2. Conferma il preset; usa **Reset al preset** per tornare ai default.
3. Inserisci opzionalmente il **Seed** per ripetibilità.
4. Scegli orizzonte/coorte e invia.
5. Scarica CSV/Plot/Twin JSON dagli artifact.

## Output
- KPI con tooltip spiegati (riduzione tumorale, perdita sani, AUC, tempo alla recidiva).
- Badge “Safe/Borderline/Unsafe” per l’interazione.
- Note cliniche per audit.
