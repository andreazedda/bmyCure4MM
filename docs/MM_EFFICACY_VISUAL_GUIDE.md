# MM Efficacy Estimation - Visual Guide

## What Users Will See

### Example 1: Bortezomib (High Efficacy)

```
┌─────────────────────────────────────────────────────────────────────┐
│ 🔴 Multiple Myeloma Efficacy Estimation                    95/100   │
│                                                                      │
│ ████████████████████████████████████████████████ 95%                │
│ Strong evidence for MM efficacy                                     │
│                                                                      │
│ Confidence Level: High              Evidence Sources:               │
│                                     • Clinical Trials                │
│                                     • PubMed Literature              │
│                                     • Target Analysis                │
│                                                                      │
│ ℹ Clinical Status: FDA approved / 12 MM trials, 3 active            │
│                                                                      │
│ Key Findings:                                                        │
│ • Targets known MM-relevant protein (score: 30/30)                  │
│ • Mechanism relevant to MM pathology (score: 25/25)                 │
│ • MM clinical trials: FDA approved (score: 35/35)                   │
│ • Strong MM literature support (score: 10/10)                       │
│                                                                      │
│ ⚠ Note: This estimation is based on target relevance, mechanism    │
│ of action, clinical trials, and published literature. It does not   │
│ replace clinical judgment or regulatory approval status.            │
└─────────────────────────────────────────────────────────────────────┘
```

**Color**: Green background, green progress bar
**Badge**: Green with "95/100"
**Interpretation**: "Strong evidence for MM efficacy"

---

### Example 2: Investigational Drug (Medium Efficacy)

```
┌─────────────────────────────────────────────────────────────────────┐
│ 🔴 Multiple Myeloma Efficacy Estimation                    58/100   │
│                                                                      │
│ █████████████████████████░░░░░░░░░░░░░░░░░░░ 58%                   │
│ Moderate evidence for MM efficacy                                   │
│                                                                      │
│ Confidence Level: Medium            Evidence Sources:               │
│                                     • Clinical Trials                │
│                                     • PubMed Literature              │
│                                                                      │
│ ℹ Clinical Status: Phase II trials (2 MM trials), 1 active          │
│                                                                      │
│ Key Findings:                                                        │
│ • Targets kinase pathway (score: 18/30)                            │
│ • MM clinical trials: Phase II (score: 20/35)                       │
│ • Moderate MM literature support (score: 6/10)                      │
│                                                                      │
│ ⚠ Limitations: No ChEMBL mechanism data available                   │
│                                                                      │
│ ⚠ Note: This estimation is based on target relevance, mechanism    │
│ of action, clinical trials, and published literature.               │
└─────────────────────────────────────────────────────────────────────┘
```

**Color**: Yellow background, yellow progress bar
**Badge**: Yellow with "58/100"
**Interpretation**: "Moderate evidence for MM efficacy"

---

### Example 3: Non-MM Drug (Low/No Efficacy)

```
┌─────────────────────────────────────────────────────────────────────┐
│ 🔴 Multiple Myeloma Efficacy Estimation                     8/100   │
│                                                                      │
│ ███░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░░ 8%                  │
│ Insufficient evidence for MM efficacy                               │
│                                                                      │
│ Confidence Level: Insufficient Data Evidence Sources:               │
│                                     • (None)                         │
│                                                                      │
│ ⚠ Limitations:                                                      │
│ • No ChEMBL mechanism data available                                │
│ • No MM clinical trial data found                                   │
│ • Limited literature evidence                                       │
│                                                                      │
│ ⚠ Note: This compound does not appear to have evidence supporting  │
│ use in Multiple Myeloma treatment.                                  │
└─────────────────────────────────────────────────────────────────────┘
```

**Color**: Gray background, gray progress bar
**Badge**: Gray with "8/100"
**Interpretation**: "Insufficient evidence for MM efficacy"

---

## Score Interpretation Guide

### Visual Indicators

| Score Range | Progress Bar Color | Badge Color | Interpretation | User Action |
|-------------|-------------------|-------------|----------------|-------------|
| **90-100** | 🟢 Dark Green | 🟢 Green | "Excellent MM efficacy evidence" | FDA-approved or equivalent |
| **70-89** | 🟢 Green | 🟢 Green | "Strong MM efficacy evidence" | Late-stage development |
| **50-69** | 🟡 Yellow | 🟡 Yellow | "Moderate MM efficacy evidence" | Mid-stage development |
| **30-49** | 🔵 Blue | 🔵 Blue | "Limited MM efficacy evidence" | Early investigation |
| **0-29** | ⚪ Gray | ⚪ Gray | "Insufficient MM efficacy evidence" | No clear MM application |

### Confidence Badges

| Level | Color | Meaning |
|-------|-------|---------|
| **High** | 🟢 Green | ≥3 evidence sources + score ≥70 |
| **Medium** | 🟡 Yellow | ≥2 evidence sources + score ≥50 |
| **Low** | 🔵 Blue | Score ≥30 but limited sources |
| **Insufficient** | ⚪ Gray | Score <30 or no evidence |

---

## Real Examples

### Bortezomib (Velcade)
- **Score**: 95/100
- **Target**: PSMB5 (proteasome)
- **Mechanism**: Proteasome inhibitor
- **Clinical**: FDA approved, 12+ trials
- **Literature**: 100+ papers
- **Display**: Green, "Strong evidence"

### Lenalidomide (Revlimid)
- **Score**: 97/100
- **Target**: CRBN (cereblon)
- **Mechanism**: IMiD
- **Clinical**: FDA approved, 15+ trials
- **Literature**: 120+ papers
- **Display**: Green, "Strong evidence"

### Daratumumab (Darzalex)
- **Score**: 88/100
- **Target**: CD38
- **Mechanism**: Monoclonal antibody
- **Clinical**: FDA approved, 8+ trials
- **Literature**: 50+ papers
- **Display**: Green, "Strong evidence"

### Selinexor (Xpovio)
- **Score**: 75/100
- **Target**: XPO1
- **Mechanism**: Export inhibitor
- **Clinical**: FDA approved, 5+ trials
- **Literature**: 30+ papers
- **Display**: Green, "Strong evidence"

### Experimental Kinase Inhibitor
- **Score**: 52/100
- **Target**: BTK kinase
- **Mechanism**: Kinase inhibitor
- **Clinical**: Phase II, 2 trials
- **Literature**: 8 papers
- **Display**: Yellow, "Moderate evidence"

### Aspirin (Control - Non-MM Drug)
- **Score**: 2/100
- **Target**: COX enzymes
- **Mechanism**: Anti-inflammatory
- **Clinical**: No MM trials
- **Literature**: No MM papers
- **Display**: Gray, "Insufficient evidence"

---

## Mobile View

On smaller screens, the layout stacks vertically:

```
┌─────────────────────────┐
│ MM Efficacy      95/100 │
│ ████████████████ 95%    │
│ Strong evidence         │
│                         │
│ Confidence: High        │
│                         │
│ Evidence:               │
│ • Clinical Trials       │
│ • Literature            │
│                         │
│ Clinical Status:        │
│ FDA approved, 12 trials │
│                         │
│ [Show Details ▼]        │
└─────────────────────────┘
```

---

## Interaction

### Hover Effects
- Progress bar scales slightly
- Badges show tooltip with breakdown
- Evidence sources are clickable (links to source)

### Click Actions
- "Show Details" expands full breakdown
- Evidence sources link to:
  - Clinical Trials → ClinicalTrials.gov search
  - PubMed Literature → PubMed search
  - ChEMBL → ChEMBL molecule page

---

## Accessibility

### Screen Reader Support
- Progress bar has `aria-label` with score
- Color is not the only indicator (text + icons)
- Semantic HTML with proper headings
- WCAG 2.1 AA compliant

### Keyboard Navigation
- All interactive elements are focusable
- Tab order follows visual order
- Enter/Space activate links

---

## Localization

### English Version
```
🔴 Multiple Myeloma Efficacy Estimation
Confidence Level: High
Evidence Sources: Clinical Trials, PubMed Literature
Note: This estimation is based on target relevance...
```

### Italian Version
```
🔴 Stima dell'Efficacia nel Mieloma Multiplo
Livello di Confidenza: Alto
Fonti di Evidenza: Trial Clinici, Letteratura PubMed
Nota: Questa stima si basa sulla rilevanza del target...
```

Both versions display simultaneously with language toggle.

---

## Performance

### Load Time
- **Additional rendering**: ~50ms
- **No API calls**: Uses cached data
- **Minimal DOM**: ~200 lines HTML
- **Total impact**: Negligible

### Caching
- Efficacy profile cached with PDB metadata
- Recalculated only when data changes
- Session storage for repeat views

---

## Summary

The MM Efficacy Estimation feature provides:

✅ **Visual clarity** - Color-coded scoring with progress bars
✅ **Evidence transparency** - Shows all contributing factors
✅ **Clinical context** - Links to trials and literature
✅ **User guidance** - Clear interpretation of results
✅ **Accessibility** - Screen reader support, keyboard nav
✅ **Bilingual** - English and Italian

Users can quickly assess MM relevance while understanding the evidence basis for the score.
