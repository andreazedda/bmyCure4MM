# MM Efficacy Estimation — Visual Guide

Questa pagina mostra **come** presentiamo il punteggio “MM efficacy estimation” in modo moderno (leggibile, scansionabile, con progress bar) e **come interpretarlo**.

## Terminologia (cosa stai guardando)

- **Score (0–100)**: indicatore sintetico di “evidenza” (non è una probabilità clinica).
- **Confidence level**: livello qualitativo in base a fonti/completezza dati.
- **Evidence sources**: categorie di segnali (trial, letteratura, meccanismo/target).
- **Limitations**: cosa manca (es. nessun meccanismo ChEMBL, pochi trial, ecc.).

## Pipeline (rendering/valutazione)

```mermaid
flowchart TD
  A[Input: compound / target] --> B[Fetch evidence sources]
  B --> C[Score components]
  C --> D[Aggregate 0..100]
  D --> E[Confidence level]
  E --> F[UI card + progress bar + findings]
```

## Esempi (come appare in UI)

### Example 1 — High efficacy (es. bortezomib)

<div class="bmy-card">
  <h3>Multiple Myeloma Efficacy Estimation <span class="bmy-badge bmy-badge--good">95 / 100</span></h3>
  <progress class="bmy-progress bmy-progress--good" max="100" value="95"></progress>
  <p><strong>Interpretation:</strong> Strong evidence for MM efficacy</p>
  <p><strong>Confidence:</strong> High</p>
  <ul>
    <li>Evidence sources: Clinical trials, PubMed literature, Target analysis</li>
    <li>Clinical status: FDA approved / many MM trials</li>
    <li>Key findings: strong target relevance + mechanism + trials</li>
  </ul>
</div>

### Example 2 — Medium efficacy (investigational)

<div class="bmy-card">
  <h3>Multiple Myeloma Efficacy Estimation <span class="bmy-badge bmy-badge--mid">58 / 100</span></h3>
  <progress class="bmy-progress bmy-progress--mid" max="100" value="58"></progress>
  <p><strong>Interpretation:</strong> Moderate evidence for MM efficacy</p>
  <p><strong>Confidence:</strong> Medium</p>
  <ul>
    <li>Evidence sources: Clinical trials, PubMed literature</li>
    <li>Clinical status: Phase II trials (limited MM signal)</li>
    <li>Limitations: missing or partial mechanism data</li>
  </ul>
</div>

### Example 3 — Low/no efficacy

<div class="bmy-card">
  <h3>Multiple Myeloma Efficacy Estimation <span class="bmy-badge bmy-badge--low">8 / 100</span></h3>
  <progress class="bmy-progress bmy-progress--low" max="100" value="8"></progress>
  <p><strong>Interpretation:</strong> Insufficient evidence for MM efficacy</p>
  <p><strong>Confidence:</strong> Insufficient data</p>
  <ul>
    <li>Evidence sources: (none)</li>
    <li>Limitations: no MM trials + limited literature + missing mechanism</li>
  </ul>
</div>

## Guida interpretazione score

| Score range | Color | Interpretazione | Azione tipica |
| --- | --- | --- | --- |
| 90–100 | green | “Excellent evidence” | farmaco consolidato / alta fiducia |
| 70–89 | green | “Strong evidence” | late-stage / molte fonti |
| 50–69 | yellow | “Moderate evidence” | mid-stage / segnale parziale |
| 30–49 | blue/neutral | “Limited evidence” | esplorativo |
| 0–29 | gray | “Insufficient data” | nessuna evidenza MM |

## Note di sicurezza e limiti

!!! warning
    Questo score è un indicatore di *evidenza* (informativo). Non sostituisce:
    - linee guida cliniche
    - judgment del clinico
    - validazione regolatoria

