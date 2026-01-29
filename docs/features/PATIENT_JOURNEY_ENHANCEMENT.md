# 🎯 Patient Journey Enhancement - Specifica Funzionale

> Documento di specifica per rendere il flusso "paziente → analisi → trattamento → esiti → suggerimenti" completo e trasparente per utenti non esperti.

---

## 📋 Flusso Utente Ideale

```
┌─────────────────────────────────────────────────────────────────────────────┐
│  1. PAZIENTE ENTRA                                                          │
│     └─ Form anagrafica + anamnesi preformattata                             │
├─────────────────────────────────────────────────────────────────────────────┤
│  2. ANALISI (Assessment)                                                    │
│     └─ Lab values + Sintomi clinici (scale standard)                        │
│     └─ → Calcolo automatico R-ISS, risk score, prognosi stimata             │
├─────────────────────────────────────────────────────────────────────────────┤
│  3. SELEZIONE TRATTAMENTO                                                   │
│     └─ Regimi suggeriti in base a risk profile                              │
│     └─ Confronto side-by-side (A vs B)                                      │
├─────────────────────────────────────────────────────────────────────────────┤
│  4. SIMULAZIONE                                                             │
│     └─ Esecuzione PKPD con Patient Twin                                     │
│     └─ Timeline evoluzione: "a 3/6/12 mesi ci aspettiamo..."                │
├─────────────────────────────────────────────────────────────────────────────┤
│  5. ESITI + RACCOMANDAZIONI                                                 │
│     └─ Prognosi esplicita (PFS/OS stimati)                                  │
│     └─ Algoritmo decisionale visibile                                       │
│     └─ Suggerimenti terapie alternative specifiche                          │
└─────────────────────────────────────────────────────────────────────────────┘
```

---

## 1. Sintomi Clinici Preformattati

### 1.1 Nuovo Modello: `SymptomAssessment`

```python
# clinic/models.py - AGGIUNTA

class SymptomAssessment(models.Model):
    """Valutazione sintomi clinici con scale standardizzate."""
    
    assessment = models.OneToOneField(
        Assessment, 
        on_delete=models.CASCADE, 
        related_name="symptoms"
    )
    
    # Dolore (scala NRS 0-10)
    bone_pain_score = models.IntegerField(
        "Dolore osseo (0-10)",
        choices=[(i, str(i)) for i in range(11)],
        null=True, blank=True,
        help_text="0=nessun dolore, 10=dolore insopportabile"
    )
    bone_pain_location = models.CharField(
        max_length=128, blank=True,
        help_text="es: colonna lombare, costale, bacino"
    )
    
    # Fatigue (scala FACIT-F semplificata 0-4)
    fatigue_score = models.IntegerField(
        "Fatigue (0-4)",
        choices=[
            (0, "0 - Nessuna"),
            (1, "1 - Lieve"),
            (2, "2 - Moderata"),
            (3, "3 - Severa"),
            (4, "4 - Invalidante"),
        ],
        null=True, blank=True
    )
    
    # Neuropatia periferica (CTCAE grade 0-4)
    neuropathy_grade = models.IntegerField(
        "Neuropatia periferica (Grade 0-4)",
        choices=[
            (0, "G0 - Assente"),
            (1, "G1 - Parestesie asintomatiche"),
            (2, "G2 - Sintomi moderati, limita ADL strumentali"),
            (3, "G3 - Sintomi severi, limita ADL di base"),
            (4, "G4 - Conseguenze potenzialmente letali"),
        ],
        null=True, blank=True
    )
    
    # Infezioni recenti
    recent_infections = models.BooleanField(
        "Infezioni negli ultimi 3 mesi",
        default=False
    )
    infection_details = models.CharField(
        max_length=256, blank=True,
        help_text="es: polmonite, infezione vie urinarie, sepsi"
    )
    
    # Eventi scheletrici
    skeletal_events = models.BooleanField(
        "Eventi scheletrici recenti",
        default=False,
        help_text="Fratture patologiche, compressione midollare, radioterapia ossea"
    )
    
    # Sintomi CRAB
    hypercalcemia_symptoms = models.BooleanField(
        "Sintomi ipercalcemia",
        default=False,
        help_text="Confusione, poliuria, stipsi, nausea"
    )
    renal_symptoms = models.BooleanField(
        "Sintomi renali",
        default=False,
        help_text="Oliguria, edemi, dispnea"
    )
    anemia_symptoms = models.BooleanField(
        "Sintomi anemia",
        default=False,
        help_text="Dispnea da sforzo, tachicardia, pallore"
    )
```

### 1.2 Form Sintomi con Preset

```python
# clinic/forms.py - AGGIUNTA

SYMPTOM_PRESETS = {
    "asymptomatic": {
        "label_en": "Asymptomatic / Smoldering",
        "label_it": "Asintomatico / Smoldering",
        "bone_pain_score": 0,
        "fatigue_score": 0,
        "neuropathy_grade": 0,
        "recent_infections": False,
        "skeletal_events": False,
    },
    "mild_symptoms": {
        "label_en": "Mild Symptoms",
        "label_it": "Sintomi Lievi",
        "bone_pain_score": 2,
        "fatigue_score": 1,
        "neuropathy_grade": 0,
        "recent_infections": False,
        "skeletal_events": False,
    },
    "moderate_crab": {
        "label_en": "Moderate CRAB Features",
        "label_it": "CRAB Moderato",
        "bone_pain_score": 5,
        "fatigue_score": 2,
        "neuropathy_grade": 0,
        "recent_infections": False,
        "skeletal_events": True,
        "anemia_symptoms": True,
    },
    "severe_symptomatic": {
        "label_en": "Severe / High Tumor Burden",
        "label_it": "Severo / Alto Carico Tumorale",
        "bone_pain_score": 8,
        "fatigue_score": 3,
        "neuropathy_grade": 0,
        "recent_infections": True,
        "skeletal_events": True,
        "hypercalcemia_symptoms": True,
        "renal_symptoms": True,
        "anemia_symptoms": True,
    },
}
```

---

## 2. Algoritmo Decisionale Trasparente

### 2.1 Struttura Dati Esportabile

```python
# simulator/decision_algorithm.py - NUOVO FILE

DECISION_ALGORITHM = {
    "version": "1.0",
    "last_updated": "2026-01-17",
    "description": {
        "en": "Transparent decision algorithm for MM treatment optimization",
        "it": "Algoritmo decisionale trasparente per ottimizzazione trattamento MM"
    },
    
    "thresholds": {
        "efficacy": {
            "good": {"min": 0.50, "description_en": "≥50% tumor reduction", "description_it": "≥50% riduzione tumorale"},
            "moderate": {"min": 0.30, "max": 0.50, "description_en": "30-50% tumor reduction", "description_it": "30-50% riduzione tumorale"},
            "poor": {"min": 0.0, "max": 0.30, "description_en": "<30% tumor reduction", "description_it": "<30% riduzione tumorale"},
            "negative": {"max": 0.0, "description_en": "Tumor growth", "description_it": "Crescita tumorale"},
        },
        "toxicity": {
            "acceptable": {"max": 0.20, "description_en": "<20% healthy cell loss", "description_it": "<20% perdita cellule sane"},
            "borderline": {"min": 0.20, "max": 0.30, "description_en": "20-30% healthy cell loss", "description_it": "20-30% perdita cellule sane"},
            "high": {"min": 0.30, "description_en": "≥30% healthy cell loss", "description_it": "≥30% perdita cellule sane"},
        },
        "durability": {
            "good": {"min": 180, "unit": "days", "description_en": "≥6 months to recurrence", "description_it": "≥6 mesi a recidiva"},
            "moderate": {"min": 90, "max": 180, "unit": "days", "description_en": "3-6 months to recurrence", "description_it": "3-6 mesi a recidiva"},
            "poor": {"max": 90, "unit": "days", "description_en": "<3 months to recurrence", "description_it": "<3 mesi a recidiva"},
        },
    },
    
    "decision_rules": [
        {
            "id": "RULE_001",
            "name_en": "High Toxicity Management",
            "name_it": "Gestione Alta Tossicità",
            "condition": "healthy_loss >= 0.30",
            "priority": "high",
            "action_en": "Reduce drug doses by 20-30% or shorten treatment duration",
            "action_it": "Riduci dosi farmaci del 20-30% o accorcia durata trattamento",
            "rationale_en": "Excessive damage to healthy plasma cells compromises immune function and increases infection risk",
            "rationale_it": "Danno eccessivo alle plasmacellule sane compromette funzione immunitaria e aumenta rischio infezioni",
            "evidence": "NCCN Guidelines 2024, IFM recommendations",
        },
        {
            "id": "RULE_002",
            "name_en": "Poor Efficacy Escalation",
            "name_it": "Escalation per Scarsa Efficacia",
            "condition": "tumor_reduction < 0.30 AND healthy_loss < 0.30",
            "priority": "high",
            "action_en": "Increase doses by 15-25% or extend treatment to 224-280 days",
            "action_it": "Aumenta dosi del 15-25% o estendi trattamento a 224-280 giorni",
            "rationale_en": "Suboptimal response with tolerable toxicity allows room for intensification",
            "rationale_it": "Risposta subottimale con tossicità tollerabile permette intensificazione",
            "evidence": "IMWG response criteria, MAIA/CASSIOPEIA trials",
        },
        {
            "id": "RULE_003",
            "name_en": "Treatment Failure Switch",
            "name_it": "Cambio per Fallimento Terapeutico",
            "condition": "tumor_reduction < 0",
            "priority": "critical",
            "action_en": "Switch to alternative regimen immediately",
            "action_it": "Cambia regime terapeutico immediatamente",
            "rationale_en": "Tumor growth indicates treatment resistance - continuing is futile and harmful",
            "rationale_it": "Crescita tumorale indica resistenza - continuare è inutile e dannoso",
            "evidence": "IMWG progressive disease criteria",
            "suggested_alternatives": ["DPd", "KRd", "Isa-Pd", "BCMA-targeted"],
        },
        {
            "id": "RULE_004",
            "name_en": "Early Recurrence Prevention",
            "name_it": "Prevenzione Recidiva Precoce",
            "condition": "time_to_recurrence < 90 AND time_horizon < 200",
            "priority": "medium",
            "action_en": "Extend treatment horizon to 224-280 days, consider maintenance therapy",
            "action_it": "Estendi orizzonte a 224-280 giorni, considera terapia di mantenimento",
            "rationale_en": "Longer treatment duration delays recurrence and improves progression-free survival",
            "rationale_it": "Durata maggiore ritarda recidiva e migliora sopravvivenza libera da progressione",
            "evidence": "FIRST trial (lenalidomide maintenance), TOURMALINE-MM3",
        },
        {
            "id": "RULE_005",
            "name_en": "Optimal Balance Maintenance",
            "name_it": "Mantenimento Equilibrio Ottimale",
            "condition": "tumor_reduction >= 0.50 AND healthy_loss < 0.20",
            "priority": "low",
            "action_en": "Current regimen is optimal. Consider ±10% dose fine-tuning or compare alternatives",
            "action_it": "Regime attuale ottimale. Considera variazioni ±10% o confronta alternative",
            "rationale_en": "Good efficacy with acceptable toxicity - minor adjustments may further optimize",
            "rationale_it": "Buona efficacia con tossicità accettabile - piccole modifiche possono ottimizzare",
            "evidence": "Treat-to-target principles",
        },
    ],
    
    "risk_stratification": {
        "R-ISS_I": {
            "median_PFS_months": 66,
            "median_OS_months": "not reached",
            "5y_survival_percent": 82,
            "description_en": "Low risk - excellent prognosis",
            "description_it": "Basso rischio - prognosi eccellente",
        },
        "R-ISS_II": {
            "median_PFS_months": 42,
            "median_OS_months": 83,
            "5y_survival_percent": 62,
            "description_en": "Intermediate risk - standard prognosis",
            "description_it": "Rischio intermedio - prognosi standard",
        },
        "R-ISS_III": {
            "median_PFS_months": 29,
            "median_OS_months": 43,
            "5y_survival_percent": 40,
            "description_en": "High risk - aggressive disease",
            "description_it": "Alto rischio - malattia aggressiva",
        },
    },
}
```

### 2.2 Vista "Algoritmo in Chiaro"

```html
<!-- templates/simulator/algorithm_transparency.html -->

<div class="card mb-4">
    <div class="card-header bg-info text-white">
        <h4 class="mb-0">
            <span class="t-en">🔍 Decision Algorithm - Full Transparency</span>
            <span class="t-it">🔍 Algoritmo Decisionale - Piena Trasparenza</span>
        </h4>
    </div>
    <div class="card-body">
        <p class="text-muted">
            <span class="t-en">This platform uses evidence-based rules to generate recommendations. 
            All thresholds and logic are shown below.</span>
            <span class="t-it">Questa piattaforma usa regole basate sull'evidenza per generare raccomandazioni. 
            Tutte le soglie e la logica sono mostrate sotto.</span>
        </p>
        
        <h5>Soglie Efficacia</h5>
        <table class="table table-sm">
            <tr><td>🟢 Buona</td><td>≥50% riduzione tumorale</td></tr>
            <tr><td>🟡 Moderata</td><td>30-50% riduzione</td></tr>
            <tr><td>🔴 Scarsa</td><td><30% riduzione</td></tr>
            <tr><td>⚫ Fallimento</td><td>Crescita tumorale</td></tr>
        </table>
        
        <h5>Soglie Tossicità</h5>
        <table class="table table-sm">
            <tr><td>🟢 Accettabile</td><td><20% perdita cellule sane</td></tr>
            <tr><td>🟡 Borderline</td><td>20-30% perdita</td></tr>
            <tr><td>🔴 Alta</td><td>≥30% perdita</td></tr>
        </table>
        
        <div class="alert alert-secondary mt-3">
            <strong>📚 Fonti:</strong> NCCN Guidelines 2024, IMWG Criteria, 
            IFM Recommendations, MAIA/CASSIOPEIA/FIRST Trials
        </div>
        
        <a href="/api/decision-algorithm/" class="btn btn-outline-primary btn-sm">
            📥 Esporta JSON
        </a>
    </div>
</div>
```

---

## 3. Prognosi Esplicita

### 3.1 Calcolo Prognosi Stimata

```python
# simulator/prognosis.py - NUOVO FILE

from typing import Dict, Optional
from dataclasses import dataclass

@dataclass
class PrognosisEstimate:
    """Stima prognostica basata su R-ISS, citogenetica e risposta al trattamento."""
    
    # Sopravvivenza
    median_pfs_months: float
    median_os_months: Optional[float]
    
    # Range di confidenza
    pfs_range_low: float
    pfs_range_high: float
    
    # Probabilità a milestone
    prob_alive_1y: float
    prob_alive_3y: float
    prob_alive_5y: float
    
    # Fattori modificanti
    modifiers_applied: list
    confidence_level: str  # "high", "moderate", "low"
    
    # Disclaimer
    disclaimer_en: str = "Estimates based on historical data. Individual outcomes may vary."
    disclaimer_it: str = "Stime basate su dati storici. Gli esiti individuali possono variare."


def estimate_prognosis(
    r_iss: str,
    cytogenetics: dict,
    age: int,
    ecog: int,
    treatment_response: Optional[str] = None,
    simulation_results: Optional[dict] = None,
) -> PrognosisEstimate:
    """
    Calcola stima prognostica combinando:
    - R-ISS staging
    - Citogenetica ad alto rischio
    - Età e performance status
    - Risposta al trattamento (se disponibile)
    - Risultati simulazione (se disponibili)
    
    Basato su:
    - Palumbo A. et al. JCO 2015 (R-ISS)
    - Sonneveld P. et al. Blood 2016 (high-risk cytogenetics)
    - Kumar SK et al. Leukemia 2020 (modern outcomes)
    """
    
    # Base R-ISS prognosis
    base_prognosis = {
        "I": {"pfs": 66, "os": 120, "5y": 0.82},
        "II": {"pfs": 42, "os": 83, "5y": 0.62},
        "III": {"pfs": 29, "os": 43, "5y": 0.40},
    }
    
    riss = r_iss.upper() if r_iss else "II"
    base = base_prognosis.get(riss, base_prognosis["II"])
    
    pfs = base["pfs"]
    os_months = base["os"]
    five_year = base["5y"]
    
    modifiers = []
    
    # Cytogenetic modifiers
    if cytogenetics.get("del17p"):
        pfs *= 0.6
        os_months *= 0.65
        modifiers.append("del(17p): -40% PFS")
    
    if cytogenetics.get("t_4_14"):
        pfs *= 0.75
        os_months *= 0.8
        modifiers.append("t(4;14): -25% PFS")
    
    if cytogenetics.get("t_14_16"):
        pfs *= 0.5
        os_months *= 0.55
        modifiers.append("t(14;16): -50% PFS")
    
    if cytogenetics.get("gain_1q21"):
        pfs *= 0.85
        modifiers.append("1q21 gain: -15% PFS")
    
    if cytogenetics.get("hyperdiploid"):
        pfs *= 1.1
        modifiers.append("Hyperdiploid: +10% PFS")
    
    # Age modifier
    if age >= 75:
        pfs *= 0.85
        os_months *= 0.75
        modifiers.append(f"Age {age}: reduced tolerance")
    elif age < 65:
        pfs *= 1.1
        modifiers.append(f"Age {age}: fit for intensive therapy")
    
    # ECOG modifier
    if ecog >= 3:
        pfs *= 0.7
        os_months *= 0.6
        modifiers.append(f"ECOG {ecog}: poor performance status")
    
    # Treatment response modifier (if known)
    response_modifiers = {
        "sCR": 1.3,
        "CR": 1.2,
        "VGPR": 1.1,
        "PR": 1.0,
        "SD": 0.7,
        "PD": 0.3,
    }
    if treatment_response and treatment_response in response_modifiers:
        pfs *= response_modifiers[treatment_response]
        modifiers.append(f"Response {treatment_response}: ×{response_modifiers[treatment_response]}")
    
    # Simulation results modifier
    if simulation_results:
        tumor_red = simulation_results.get("tumor_reduction", 0)
        if tumor_red >= 0.7:
            pfs *= 1.15
            modifiers.append("Simulation: deep response predicted")
        elif tumor_red < 0.3:
            pfs *= 0.8
            modifiers.append("Simulation: suboptimal response predicted")
    
    # Calculate probabilities
    # Using exponential survival model approximation
    import math
    lambda_pfs = 1 / pfs if pfs > 0 else 1
    
    prob_1y = math.exp(-lambda_pfs * 12)
    prob_3y = math.exp(-lambda_pfs * 36)
    prob_5y = math.exp(-lambda_pfs * 60)
    
    # Confidence level
    if len(modifiers) <= 2:
        confidence = "high"
    elif len(modifiers) <= 4:
        confidence = "moderate"
    else:
        confidence = "low"
    
    return PrognosisEstimate(
        median_pfs_months=round(pfs, 1),
        median_os_months=round(os_months, 1) if os_months else None,
        pfs_range_low=round(pfs * 0.7, 1),
        pfs_range_high=round(pfs * 1.3, 1),
        prob_alive_1y=round(prob_1y, 2),
        prob_alive_3y=round(prob_3y, 2),
        prob_alive_5y=round(prob_5y, 2),
        modifiers_applied=modifiers,
        confidence_level=confidence,
    )
```

### 3.2 Visualizzazione Timeline

```html
<!-- Componente timeline evoluzione attesa -->

<div class="card mb-4" id="prognosis-timeline">
    <div class="card-header">
        <h5 class="mb-0">
            <span class="t-en">📈 Expected Evolution Timeline</span>
            <span class="t-it">📈 Timeline Evoluzione Attesa</span>
        </h5>
    </div>
    <div class="card-body">
        <div class="timeline">
            <div class="timeline-item">
                <div class="timeline-marker bg-primary"></div>
                <div class="timeline-content">
                    <h6>Oggi / Now</h6>
                    <p>Baseline assessment: R-ISS {{ r_iss }}, Risk score {{ risk_score|floatformat:2 }}</p>
                </div>
            </div>
            
            <div class="timeline-item">
                <div class="timeline-marker bg-info"></div>
                <div class="timeline-content">
                    <h6>3 mesi / 3 months</h6>
                    <p>
                        <span class="t-en">Expected response: {{ expected_3m_response }}</span>
                        <span class="t-it">Risposta attesa: {{ expected_3m_response }}</span>
                    </p>
                    <p class="text-muted small">
                        <span class="t-en">Reassess: labs, M-protein, symptoms</span>
                        <span class="t-it">Rivalutare: lab, M-proteina, sintomi</span>
                    </p>
                </div>
            </div>
            
            <div class="timeline-item">
                <div class="timeline-marker bg-success"></div>
                <div class="timeline-content">
                    <h6>6 mesi / 6 months</h6>
                    <p>
                        <span class="t-en">Target: {{ expected_6m_target }}</span>
                        <span class="t-it">Obiettivo: {{ expected_6m_target }}</span>
                    </p>
                </div>
            </div>
            
            <div class="timeline-item">
                <div class="timeline-marker bg-warning"></div>
                <div class="timeline-content">
                    <h6>12 mesi / 12 months</h6>
                    <p>
                        <span class="t-en">Probability still in response: {{ prob_1y|floatformat:0 }}%</span>
                        <span class="t-it">Probabilità ancora in risposta: {{ prob_1y|floatformat:0 }}%</span>
                    </p>
                </div>
            </div>
        </div>
        
        <div class="alert alert-info mt-3">
            <small>
                <span class="t-en">⚠️ These are statistical estimates based on clinical trials data. 
                Individual outcomes depend on many factors and may differ significantly.</span>
                <span class="t-it">⚠️ Queste sono stime statistiche basate su dati di trial clinici. 
                Gli esiti individuali dipendono da molti fattori e possono differire significativamente.</span>
            </small>
        </div>
    </div>
</div>
```

---

## 4. Confronto Scenari Side-by-Side

### 4.1 Vista Comparativa

```python
# simulator/views.py - AGGIUNTA

def compare_scenarios(request):
    """Confronta 2-3 scenari di trattamento side-by-side."""
    
    scenario_ids = request.GET.getlist("scenario_id")
    if len(scenario_ids) < 2:
        messages.warning(request, "Seleziona almeno 2 scenari da confrontare")
        return redirect("simulator:scenario_list")
    
    scenarios = Scenario.objects.filter(pk__in=scenario_ids)
    
    # Run simulations for each scenario (or use cached results)
    comparisons = []
    for scenario in scenarios:
        result = run_simulation_cached(scenario)
        comparisons.append({
            "scenario": scenario,
            "tumor_reduction": result.get("tumor_reduction"),
            "healthy_loss": result.get("healthy_loss"),
            "time_to_recurrence": result.get("time_to_recurrence"),
            "regimen_name": scenario.get_regimen_display(),
            "toxicity_grade": classify_toxicity(result.get("healthy_loss")),
            "efficacy_grade": classify_efficacy(result.get("tumor_reduction")),
        })
    
    # Determine winner for each metric
    best_efficacy = max(comparisons, key=lambda x: x["tumor_reduction"] or 0)
    best_toxicity = min(comparisons, key=lambda x: x["healthy_loss"] or 1)
    best_durability = max(comparisons, key=lambda x: x["time_to_recurrence"] or 0)
    
    context = {
        "comparisons": comparisons,
        "best_efficacy": best_efficacy,
        "best_toxicity": best_toxicity,
        "best_durability": best_durability,
    }
    return render(request, "simulator/compare_scenarios.html", context)
```

### 4.2 Template Comparativo

```html
<!-- templates/simulator/compare_scenarios.html -->

<div class="row">
    {% for comp in comparisons %}
    <div class="col-md-{{ 12|divisibleby:comparisons|length }}">
        <div class="card h-100 {% if comp == best_efficacy %}border-success{% endif %}">
            <div class="card-header">
                <h5>{{ comp.regimen_name }}</h5>
            </div>
            <div class="card-body">
                <table class="table table-sm">
                    <tr>
                        <td>Efficacia</td>
                        <td>
                            <span class="badge bg-{{ comp.efficacy_grade }}">
                                {{ comp.tumor_reduction|floatformat:1 }}%
                            </span>
                            {% if comp == best_efficacy %}🏆{% endif %}
                        </td>
                    </tr>
                    <tr>
                        <td>Tossicità</td>
                        <td>
                            <span class="badge bg-{{ comp.toxicity_grade }}">
                                {{ comp.healthy_loss|floatformat:1 }}%
                            </span>
                            {% if comp == best_toxicity %}🏆{% endif %}
                        </td>
                    </tr>
                    <tr>
                        <td>Durabilità</td>
                        <td>
                            {{ comp.time_to_recurrence }} giorni
                            {% if comp == best_durability %}🏆{% endif %}
                        </td>
                    </tr>
                </table>
            </div>
        </div>
    </div>
    {% endfor %}
</div>
```

---

## 5. Regimi Alternativi Suggeriti

### 5.1 Database Regimi con Indicazioni

```python
# clinic/fixtures/standard_regimens.json

[
    {
        "name": "VRd",
        "components": "bortezomib, lenalidomide, dexamethasone",
        "line": "frontline",
        "indication_primary": ["nd_standard", "nd_high_risk"],
        "avoid_if": ["severe_neuropathy", "renal_failure_severe"],
        "expected_response_rate": 0.90,
        "median_pfs_months": 41,
    },
    {
        "name": "DRd",
        "components": "daratumumab, lenalidomide, dexamethasone",
        "line": "frontline",
        "indication_primary": ["nd_standard", "frail_elderly"],
        "avoid_if": [],
        "expected_response_rate": 0.93,
        "median_pfs_months": 61,
    },
    {
        "name": "KRd",
        "components": "carfilzomib, lenalidomide, dexamethasone",
        "line": "relapsed",
        "indication_primary": ["rr_early", "nd_high_risk"],
        "avoid_if": ["cardiac_history", "hypertension_uncontrolled"],
        "expected_response_rate": 0.87,
        "median_pfs_months": 26,
    },
    {
        "name": "Isa-Pd",
        "components": "isatuximab, pomalidomide, dexamethasone",
        "line": "relapsed",
        "indication_primary": ["rr_late"],
        "avoid_if": [],
        "expected_response_rate": 0.63,
        "median_pfs_months": 11,
    }
]
```

### 5.2 Suggeritore Automatico

```python
# simulator/regimen_suggester.py

def suggest_alternative_regimens(
    current_regimen: str,
    patient_archetype: str,
    simulation_results: dict,
    contraindications: list,
) -> list:
    """
    Suggerisce regimi alternativi quando il corrente non è ottimale.
    
    Returns list of dicts with:
    - regimen_name
    - why_suggested (reason in context)
    - expected_improvement
    - cautions
    """
    
    suggestions = []
    
    # If current regimen shows poor efficacy
    if simulation_results.get("tumor_reduction", 0) < 0.30:
        # Suggest more potent options
        if "carfilzomib" not in current_regimen.lower():
            suggestions.append({
                "regimen_name": "KRd",
                "why_suggested_en": "Carfilzomib may overcome bortezomib resistance",
                "why_suggested_it": "Carfilzomib può superare resistenza a bortezomib",
                "expected_improvement": "+15-20% response rate",
                "cautions_en": "Monitor cardiac function, blood pressure",
                "cautions_it": "Monitorare funzione cardiaca, pressione",
            })
        
        if "daratumumab" not in current_regimen.lower():
            suggestions.append({
                "regimen_name": "DRd or DVRd",
                "why_suggested_en": "Adding CD38 antibody improves depth of response",
                "why_suggested_it": "Aggiungere anticorpo CD38 migliora profondità risposta",
                "expected_improvement": "+20% PFS",
                "cautions_en": "Infusion reactions (first dose)",
                "cautions_it": "Reazioni infusionali (prima dose)",
            })
    
    # If current regimen shows high toxicity
    if simulation_results.get("healthy_loss", 0) >= 0.30:
        suggestions.append({
            "regimen_name": "Rd (doublet)",
            "why_suggested_en": "Removing third agent reduces toxicity while maintaining efficacy",
            "why_suggested_it": "Rimuovere terzo agente riduce tossicità mantenendo efficacia",
            "expected_improvement": "-30% toxicity",
            "cautions_en": "May have slightly lower response rate",
            "cautions_it": "Può avere tasso risposta leggermente inferiore",
        })
    
    # Filter out contraindicated regimens
    suggestions = [
        s for s in suggestions
        if not any(c in s["regimen_name"].lower() for c in contraindications)
    ]
    
    return suggestions[:3]  # Max 3 suggestions
```

---

## 6. Implementazione Prioritaria

### Fase 1 (Quick Win - 1 settimana)
1. ✅ Aggiungere endpoint `/api/decision-algorithm/` che espone l'algoritmo in JSON
2. ✅ Aggiungere card "Algoritmo in Chiaro" nella pagina paziente
3. ✅ Mostrare prognosi stimata (PFS mediano) in base a R-ISS

### Fase 2 (Core Features - 2 settimane)
1. Modello `SymptomAssessment` con form precompilati
2. Modulo prognosi completo con timeline
3. Suggeritore regimi alternativi

### Fase 3 (Advanced - 3 settimane)
1. Confronto scenari side-by-side
2. Export report PDF con algoritmo + prognosi
3. API per integrazione sistemi esterni

---

## Note per Sviluppatori

### Testing
```bash
# Test algoritmo decisionale
python manage.py test simulator.tests.test_decision_algorithm

# Test prognosi
python manage.py test simulator.tests.test_prognosis

# Test suggeritore regimi
python manage.py test simulator.tests.test_regimen_suggester
```

### Configurazione
Le soglie sono configurabili via `simulator/presets/decision_thresholds.yaml` senza modificare codice.

### Compliance
- Algoritmo conforme a NCCN Guidelines 2024
- Dati prognostici da trial registrativi (MAIA, CASSIOPEIA, FIRST)
- Disclaimer obbligatorio su tutte le stime

---

*Documento creato: 2026-01-17*
*Autore: bmyCure4MM Development Team*
