# Simulatore di ricerca

Stato: documento corrente secondo `claims-policy-v1`.

Il simulatore confronta scenari meccanicistici sotto assunzioni esplicite. Gli
output relativi a tumore, cellule sane, esposizione, tossicità e recidiva sono
`SIMULATED`. Le soglie associate sono segnali diagnostici `HEURISTIC` relativi
al modello.

## Cosa può mostrare una simulazione

- traiettorie prodotte da modello, parametri e calendario selezionati;
- superamento di una soglia relativa al modello;
- sensibilità a calendari, parametri e assunzioni di incertezza;
- artefatti riproducibili quando identità di modello, input e versione sono
  registrate.

## Cosa non può mostrare

- quale terapia debba ricevere un paziente;
- se una dose o un calendario debbano cambiare;
- beneficio o prognosi individuale validati;
- un effetto causale del trattamento.

Usare `/learn/` per tutorial sintetici e `/research/` per flussi di ricerca con
lineage. Vedere [Uso previsto canonico](../governance/INTENDED_USE.md) e
[Linguaggio degli output](../governance/MODEL_OUTPUT_LANGUAGE.md).
