# bmyCure4MM — Patient-Specific Digital Twin: Audit, Target Architecture & Roadmap

> **Versione**: 1.0.0  
> **Data**: 2026-03-19  
> **Tipo**: Documento progettuale — audit + architettura target + backlog  
> **Scope**: Evoluzione della piattaforma da simulatore educativo a sistema digital twin paziente-specifico per il Mieloma Multiplo

---

## Indice

1. [Executive Summary](#executive-summary)
2. [Audit architetturale corrente](#audit-architetturale-corrente)
   - 2.1 [Layer dati paziente reale](#21-layer-dati-paziente-reale)
   - 2.2 [Layer follow-up longitudinale](#22-layer-follow-up-longitudinale)
   - 2.3 [Layer twin corrente](#23-layer-twin-corrente)
   - 2.4 [Layer simulazione matematica](#24-layer-simulazione-matematica)
   - 2.5 [Layer scenario/educativo](#25-layer-scenarioeducativo)
   - 2.6 [Layer decision support](#26-layer-decision-support)
   - 2.7 [Provenance e artifact](#27-provenance-e-artifact)
3. [Feature Verification Matrix](#feature-verification-matrix)
4. [Current vs Target Twin Assessment](#current-vs-target-twin-assessment)
5. [Architettura target](#architettura-target)
   - 5.1 [Layer 1 — Real Patient](#layer-1--real-patient)
   - 5.2 [Layer 2 — Twin State](#layer-2--twin-state)
   - 5.3 [Layer 3 — Observation](#layer-3--observation)
   - 5.4 [Layer 4 — Simulation](#layer-4--simulation)
   - 5.5 [Layer 5 — Counterfactual](#layer-5--counterfactual)
   - 5.6 [Layer 6 — Decision Support](#layer-6--decision-support)
   - 5.7 [Layer 7 — Reporting & Provenance](#layer-7--reporting--provenance)
6. [Minimal Viable Real-Patient Digital Twin](#minimal-viable-real-patient-digital-twin)
7. [Modelli Django proposti](#modelli-django-proposti)
8. [Architettura moduli/servizi proposta](#architettura-moduliservizi-proposta)
9. [Piano di integrazione con il repo esistente](#piano-di-integrazione-con-il-repo-esistente)
10. [Backlog GitHub Issues](#backlog-github-issues)
11. [Documentation Roadmap](#documentation-roadmap)
12. [Conclusione tecnica finale](#conclusione-tecnica-finale)

---

## Executive Summary

**bmyCure4MM** è una piattaforma Django per il Mieloma Multiplo composta da quattro blocchi funzionali:

1. **Layer dati paziente reale** (`clinic/models.py`): modelli `Patient`, `Assessment`, `PatientTherapy`, `PatientCytogenetics`, `SymptomAssessment`. CRUD completo, timeline terapeutica, biomarker IMWG, citogenetica, sintomi con scale cliniche validate (NRS, FACIT-F, CTCAE).

2. **"Twin" leggero** (`simulator/twin.py`): funzione stateless `build_patient_twin(assessment)` che mappa 4 biomarker (R-ISS, LDH, β2M, FLC ratio) in un `risk_score` scalare, poi interpola linearmente in parametri del simulatore (growth rate, carrying capacity, immune compromise index). **Non è un digital twin dinamico**: è un layer di derivazione parametrica one-shot senza stato, senza persistenza, senza aggiornamento temporale, senza modello di osservazione.

3. **Simulatore ODE** (`simulator/models_simulation.py`, `simulator/models.py`): sistema accoppiato PK/PD + dinamica tumorale logistica, risolto con `scipy.integrate.solve_ivp`, supporto per 3 farmaci (lenalidomide, bortezomib, daratumumab) con schedule posologici configurabili via YAML, matrice d'interazione, replicazione stocastica per coorte.

4. **Decision support euristico** (`clinic/views.py`, `simulator/decision_algorithm.py`): interpretazione a soglie fisse dei KPI di simulazione (tumor reduction ≥0.50, healthy loss <0.20, TTR ≥180 giorni) con raccomandazioni bilingui (EN/IT). Interamente a regole, nessun ragionamento model-based paziente-specifico.

**Diagnosi**: La piattaforma è un **workbench di simulazione educativa con layer opzionale di risk-mapping biomarker→parametri**. NON è un sistema digital twin paziente-specifico. Il gap verso un vero digital twin è strutturale: mancano stato latente, modello di osservazione, loop di aggiornamento, calibrazione, confronto predetto-vs-osservato, e branching controfattuale.

---

## Audit architetturale corrente

### 2.1 Layer dati paziente reale

**Status: Implementato.**

Il modulo `clinic` fornisce un layer CRUD completo per dati paziente reali.

#### Modelli implementati

| Modello | File | Scopo | Campi chiave |
|---------|------|-------|--------------|
| `Patient` | `clinic/models.py:9` | Demografiche paziente | `mrn` (unique), `first_name`, `last_name`, `birth_date`, `sex`, `diagnosis_date`, `owner` (FK → User), `notes` |
| `Assessment` | `clinic/models.py:83` | Valori di laboratorio IMWG | `patient` (FK), `date`, `m_protein_g_dl`, `kappa_mg_l`, `lambda_mg_l`, `flc_ratio`, `hemoglobin_g_dl`, `calcium_mg_dl`, `creatinine_mg_dl`, `beta2m_mg_l`, `albumin_g_dl`, `ldH_u_l`, `r_iss`, `response` (IMWG), `notes` |
| `CytogeneticAbnormality` | `clinic/models.py:52` | Catalogo anomalie citogenetiche | `code` (unique), `description` |
| `PatientCytogenetics` | `clinic/models.py:62` | Anomalie per paziente | `patient` (FK), `abnormality` (FK), `detected_on`, `method` |
| `Regimen` | `clinic/models.py:121` | Regimi terapeutici | `name`, `line`, `components`, `intent`, `notes` |
| `PatientTherapy` | `clinic/models.py:136` | Timeline terapeutica | `patient` (FK), `regimen` (FK), `start_date`, `end_date`, `outcome`, `adverse_events` |
| `SymptomAssessment` | `clinic/models_symptoms.py:18` | Tracking sintomi strutturato | `patient` (FK), `assessment_date`, `bone_pain_nrs` (NRS 0-10), `fatigue_total` (FACIT-F 0-52), `neuropathy_grade` (CTCAE 0-4), flag CRAB (hypercalcemia, renal, anemia, bone), `infection_details`, `ecog_status` |

#### Controllo accessi

Implementato tramite:
- `Patient.owner` FK → User (`clinic/models.py:19`)
- Logica `can_edit_patient()` in `clinic/views.py:207`: staff vede tutto, non-staff vede solo pazienti posseduti + DEMO
- `simulator/access.py`: funzioni `accessible_assessments(user)` e `get_accessible_assessment_by_id(user, id)` con policy DEMO_MRN_PREFIX

#### API REST

`clinic/api.py` espone `Patient` e `Assessment` via DRF `ModelViewSet` con serializzazione completa.

#### Flussi view

- `dashboard`: statistiche aggregate, distribuzione R-ISS, pazienti che necessitano attenzione (`clinic/views.py:213`)
- `patient_list`: filtri per R-ISS e citogenetica high-risk, ricerca per cognome (`clinic/views.py:265`)
- `patient_detail`: grafici temporali biomarker (Plotly JSON), timeline terapeutica con overlay, effetti terapia (delta M-protein, delta risk_score), simulazione più recente collegata (`clinic/views.py:292`)
- `assessment_new`: creazione assessment (`clinic/views.py:597`)
- `prognosis_timeline`: stima sopravvivenza R-ISS-based (`clinic/views.py:660`)

---

### 2.2 Layer follow-up longitudinale

**Status: Strutturalmente supportato, ma non consumato longitudinalmente dal twin/simulazione.**

Il modello `Assessment` supporta osservazioni ripetute per paziente con `ordering = ["-date"]`. La view `patient_detail` renderizza grafici temporali e calcola delta per-terapia (`therapy_effects`):

```python
# clinic/views.py:398-440
# Per ogni PatientTherapy:
#   baseline = ultimo assessment <= start_date
#   follow = ultimo assessment <= end_date|ora
#   delta_m = follow.m_protein - baseline.m_protein
#   delta_risk = build_patient_twin(follow).risk_score - build_patient_twin(baseline).risk_score
```

**Limiti critici:**
- `build_patient_twin()` accetta esattamente **un** `Assessment` alla volta (`simulator/twin.py:11`)
- Nessun codice aggrega la storia assessment in uno stato twin multi-timepoint
- Nessun meccanismo confronta traiettorie predette vs assessment reali successivi
- La timeline terapeutica è renderizzata nella patient detail ma **mai consumata come input di simulazione** — il simulatore usa i propri `Scenario.drug_doses` e `SimulationAttempt.parameters`

---

### 2.3 Layer twin corrente

**Status: Implementato come funzione stateless di risk-mapping. Non è un digital twin dinamico.**

L'intera logica twin risiede in `simulator/twin.py` — una singola funzione:

```python
def build_patient_twin(assessment: Assessment) -> Dict[str, float]:
```

#### Caratterizzazione precisa

| Domanda | Risposta | Evidenza |
|---------|----------|----------|
| È statico? | **Sì.** Funzione pura senza side-effect. | `simulator/twin.py:11-46` — nessun accesso DB in scrittura |
| Basato su un assessment? | **Sì.** Accetta un singolo oggetto `Assessment`. | Firma funzione, nessuna queryset |
| Produce parametri simulatore? | **Sì.** Restituisce `risk_score`, `tumor_growth_rate`, `healthy_growth_rate`, `carrying_capacity_tumor`, `carrying_capacity_healthy`, `immune_compromise_index` | `simulator/twin.py:36-45` |
| È un modello di stato dinamico? | **No.** Nessuna variabile latente, nessuna ODE, nessuna dinamica temporale. | Assenza completa di queste strutture nel file |
| È persistito? | **No.** Il dict twin è calcolato on-the-fly. Un file JSON è salvato come artifact di una simulation run. | `simulator/models.py:680` — `twin_params_{pk}.json` — nessun modello Django `PatientTwinState` |
| È ricalibrato nel tempo? | **No.** Ogni chiamata mappa indipendentemente i biomarker di un assessment a parametri. | Nessun confronto con parametri twin precedenti |
| Confronta predetto vs osservato? | **No.** | Assenza di qualsiasi logica di residui |

#### Meccanismo di mapping

Configurazione YAML (`simulator/presets/twin_risk.yaml`):

```yaml
weights:
  riss: 0.35       # 35%
  ldh: 0.2         # 20%
  beta2m: 0.25     # 25%
  flc_ratio: 0.2   # 20%

growth_mapping:
  tumor_growth_rate: {min: 0.017, max: 0.04}
  healthy_growth_rate: {min: 0.012, max: 0.02}
carrying_capacity:
  tumor: {min: 5.0e8, max: 2.5e9}
  healthy: {min: 2.5e11, max: 6.0e11}
immune_compromise_index: {min: 0.85, max: 1.3}
```

Flusso: Assessment → 4 component score (R-ISS, LDH, β2M, FLC) → weighted sum → risk_score ∈ [0,1] → interpolazione lineare (lerp) in parametri del solver.

#### Path di utilizzo

`SimulationAttempt.run_model()` (`simulator/models.py:335`) chiama opzionalmente `build_patient_twin(assessment)` quando `use_twin=True` e `twin_assessment_id` è valido, poi fonde i parametri restituiti negli input del solver via `_merge_twin_parameters()`.

L'endpoint API `simulator/api_twin.py:twin_preview` restituisce una risposta JSON con input assessment e parametri twin derivati, ma **non effettua persistenza né aggiornamento**.

---

### 2.4 Layer simulazione matematica

**Status: Implementato — sistema ODE accoppiato PK/PD + dinamica tumorale/sana.**

Esistono due implementazioni parallele:

#### Solver primario: `MathematicalModel` dataclass (`simulator/models_simulation.py`)

**Vettore di stato ODE:** `[tumor_cells, healthy_cells, drug_1_concentration, ..., drug_N_concentration]`

**Equazioni:**

$$\frac{dT}{dt} = r_T \cdot T \cdot \left(1 - \frac{T}{K_T}\right) - \sum_i E_i \cdot T$$

$$\frac{dH}{dt} = r_H \cdot H \cdot \left(1 - \frac{H}{K_H}\right) - \sigma_{immune} \cdot \bar{E} \cdot H$$

$$\frac{dC_i}{dt} = -k_{elim,i} \cdot C_i + f_{dose,i}(t)$$

$$E_i = \frac{E_{max,i} \cdot C_i}{EC_{50,i} + C_i} \quad \text{(Emax model)}$$

Dove:
- $T$ = cellule tumorali, $H$ = cellule sane
- $r_T, r_H$ = growth rate, $K_T, K_H$ = carrying capacity
- $C_i$ = concentrazione farmaco $i$, $E_i$ = effetto PD
- $\sigma_{immune}$ = immune compromise index
- $\bar{E}$ = media degli effetti farmaco
- $f_{dose,i}(t)$ = funzione schedule posologico

**Solver:** `scipy.integrate.solve_ivp` con RK45, `rtol=1e-6`, `atol=1e-8`, `t_eval = 200 punti`

**Output:** `pd.DataFrame` con colonne `time`, `tumor_cells`, `healthy_cells`, `{drug}_concentration`

#### Modello scientifico esteso: `MultipleMyelomaODE` (`simulator/mathematical_models.py`)

**Non collegato al flusso principale.** Usa:
- Crescita tumorale **Gompertziana** ($dT/dt = r_T \cdot T \cdot \ln(K_T/T)$) anziché logistica
- Modello di sorveglianza immunitaria (`ImmuneResponseModel`)
- Dataclass PK/PD strutturate con coefficiente di Hill
- Modello interazione farmaci secondo Greco et al.
- `scipy.integrate.odeint`

`SimulationAttempt.run_model()` usa `MathematicalModel`, **non** `MultipleMyelomaODE`. Il modello esteso è disponibile ma non integrato.

#### Supporto schedule farmaci

Il registro farmaci (`simulator/pharmaco/registry.py`) carica preset YAML per-farmaco da `simulator/presets/drugs/`:
- Validazione profilo: richiede `half_life`, `Vd`, `Emax`, `EC50`, `dose_range`
- Tipi schedule: `continuous`, `cycle` (on/off ciclico), `weekly` (giorni specifici), `interval` (intervallo fisso), `pulsed`
- Costruisce funzione callable `ScheduleFn = Callable[[float], float]` per l'ODE solver

#### Metriche summary

`_summarize_trajectory()` (`simulator/models.py:400`) calcola:
- `tumor_reduction`: $1 - T_{end}/T_{start}$
- `healthy_loss`: $1 - H_{end}/H_{start}$
- `time_to_recurrence`: primo tempo post-nadir con $T > 0.5 \cdot T_{start}$
- `durability_index`: $\int \mathbb{1}(T < T_{start}) dt / T_{horizon}$
- `nadir`: giorno, cellule, riduzione percentuale
- `milestones`: snapshot a giorno 30, 60, 90, fine
- `auc`: area sotto la curva per-farmaco ($\int C_i(t) dt$)

#### Simulazione coorte

Due meccanismi:
1. **Interno a `run_model()`** (`simulator/models.py:510`): quando `cohort_size > 1`, perturbazione log-normale dei parametri biologici (σ=0.25–0.50), produce media/p05/p95 delle distribuzioni
2. **Modulo dedicato** (`simulator/cohort.py`): genera popolazioni sintetiche con `sample_patient_params()`, esegue simulazioni aggregate

#### Ottimizzazione multi-obiettivo

`simulator/optim.py`: usa Optuna `TPESampler` per ottimizzazione Pareto su efficacy, safety, exposure. Lo spazio di ricerca è definito in `simulator/search_space.py`.

#### Modulo design

`simulator/design/` implementa simulazione research-grade per modalità avanzate (small molecule, ADC, CAR-T):
- Logica gate antigenica (`design/targeting.py`)
- Evoluzione tumorale con escape antigenico (`design/evolution.py`)
- Decomposizione tossicità (`design/toxicity.py`)
- Architetturalmente indipendente dal pipeline di simulazione principale

---

### 2.5 Layer scenario/educativo

**Status: Implementato — separato dai dati paziente reale, collegato via bridge twin assessment.**

#### Modello `Scenario` (`simulator/models.py:60`)

Memorizza casi educativi con:
- Stadio clinico, R-ISS, citogenetica (6 flag boolean: `del17p`, `t_4_14`, `t_14_16`, `gain_1q21`, `hyperdiploid`, `t_11_14`)
- Parametri biologia tumorale (`tumor_cell_count`, `tumor_growth_rate`, `carrying_capacity`)
- Caratteristiche paziente (`patient_age`, `ecog_performance_status`, `charlson_comorbidity_index`, `patient_archetype`)
- Valori laboratorio (`creatinine_clearance`, `albumin`, `beta2_microglobulin`, `ldh`, `hemoglobin`, `calcium`)
- Regimi raccomandati (M2M → `clinic.Regimen`)
- Note linee guida, risposta attesa IMWG
- Difficulty scoring calcolato (`difficulty_score`, `difficulty_level`)

#### Modello `SimulationAttempt` (`simulator/models.py:280`)

Collega l'azione utente a uno scenario:
- `scenario` (FK), `user` (FK), `selected_regimen` (FK)
- `parameters` (JSONField), `results` (JSONField), `results_summary` (JSONField)
- `artifacts` (JSONField), `seed`, `cohort_size`, `submitted`, `confidence`

#### Bridge twin

Il parametro `twin_assessment_id` in `SimulationAttempt.parameters` crea il ponte: un assessment di un paziente reale può essere usato per derivare parametri twin dentro una simulazione di scenario. Questo è **l'unico collegamento** tra dati paziente reale e esecuzione simulazione.

#### Virtual patients

`simulator/virtual_patients.py`: generatore sintetico basato su archetipi (`PatientArchetype`) con distribuzioni statistiche per tutti i parametri clinici. Usato per trial in-silico, non per pazienti reali.

---

### 2.6 Layer decision support

**Status: Implementato — euristico, basato su soglie.**

#### Implementazioni

1. **Interpretazione contesto-paziente**: `_interpret_latest_simulation()` (`clinic/views.py:21`)
   - Input: dict summary simulazione
   - Soglie: `tumor_reduction ≥ 0.50` = good, `healthy_loss < 0.20` = acceptable, `TTR ≥ 180` = good
   - Output: raccomandazioni prioritizzate bilingui con azioni, razionali, icone
   - Usato nella view `patient_detail`

2. **Trasparenza algoritmica**: `simulator/decision_algorithm.py`
   - Dict frozen con tutte le soglie, regole, evidenze
   - Versione `1.0.0`, 5 regole decisionali (RULE_001–RULE_005)
   - Esposto via API `simulator/api_algorithm.py` e pagina trasparenza

3. **Suggerimento regimi**: `simulator/regimen_suggester.py`
   - Database curato di regimi frontline e relapsed con trial evidence e controindicazioni
   - 7+ regimi frontline (VRd, Dara-VRd, KRd, DRd, VMP, Dara-VMP, Rd)
   - Regimi relapsed (DPd, DVd, etc.)

4. **Stima prognostica**: `simulator/prognosis.py`
   - Stime sopravvivenza R-ISS-based (Palumbo et al. JCO 2015)
   - Hazard ratio citogenetici (IMWG 2016, mSMART 3.0)
   - Modificatori età, ECOG, risposta terapeutica, MRD
   - Output: `PrognosisEstimate` dataclass con PFS/OS a 12/24/36/60 mesi

**Tutto il decision support è euristico**: soglie hardcoded e logica a regole. Nessuna inversione model-based, nessuna probabilità a posteriori, nessun intervallo di confidenza da calibrazione.

---

### 2.7 Provenance e artifact

**Status: Parziale.**

#### Cosa viene salvato per run

- `SimulationAttempt.parameters` (JSONField): parametri completi della simulazione
- `SimulationAttempt.results` (JSONField): URL file CSV/HTML/JSON generati
- `SimulationAttempt.results_summary` (JSONField): metriche summary
- `SimulationAttempt.artifacts` (JSONField): mappa nome→URL
- `SimulationAttempt.seed`: seme per riproducibilità
- `SimulationAttempt.submitted`: timestamp
- `SimulationAttempt.user`: FK utente

#### File generati

| Tipo | Path | Contenuto |
|------|------|-----------|
| CSV trajectory | `media/sim_data/attempt_{pk}.csv` | Colonne time, tumor_cells, healthy_cells, {drug}_concentration |
| Plot HTML | `media/sim_plots/attempt_{pk}.html` | Plotly full-HTML con JS inline, 2 subplot |
| Twin params JSON | `media/sim_data/attempt_{pk}_twin.json` | Parametri twin + assessment_id |

#### Cosa manca

- **Nessun versioning del modello**: nessun tracking della versione del codice o del modello usato
- **Nessun versioning delle configurazioni**: nessun hash di `twin_risk.yaml` o dei drug preset al momento della run
- **Nessun `ModelVersion` model**: nessun meccanismo per congelare e tracciare versioni
- **Nessuna catena di provenance formale**: nessun FK che colleghi patient → twin derivation → simulation → recommendation
- **Nessun `ProvenanceArtifact` model**: nessun modello dedicato
- **Riproducibilità parziale**: dipende dal seed salvato, ma non dalla registrazione della versione esatta del codebase, della config twin, dei preset farmaci

---

## Feature Verification Matrix

| Feature | Evidenza nel codice | Status | Note |
|---------|-------------------|--------|------|
| **Demografiche paziente** | `clinic.Patient` model | ✅ Implementato | MRN, nome, DOB, sesso, data diagnosi, owner FK |
| **Timeline assessment** | `clinic.Assessment` con `ordering=["-date"]` | ✅ Implementato | Osservazioni ripetute, pannello IMWG completo |
| **Storico biomarker** | `patient_detail` view con grafici temporali | ✅ Implementato | M-protein, FLC, LDH, β2M plottati nel tempo |
| **Timeline terapeutica** | `clinic.PatientTherapy` model + rendering dettaglio | ✅ Implementato | Regimen start/end/outcome/AE; overlay nel grafico |
| **Citogenetica** | `PatientCytogenetics` model | ✅ Implementato | Per-paziente, per-data, metodo registrato |
| **Generazione twin paziente** | `build_patient_twin()` in `simulator/twin.py` | ⚠️ Parziale | Risk-mapping da singolo assessment; non è un modello di stato |
| **Persistenza twin** | Artifact file JSON salvato per run | ⚠️ Parziale | JSON come side-effect della simulazione; nessun modello Django |
| **Aggiornamento twin nel tempo** | Non implementato | ❌ Mancante | Nessun codice aggiorna uno stato twin da nuovi assessment |
| **Modello PK** | 1-compartimento, eliminazione 1° ordine in `MathematicalModel` | ✅ Implementato | Schedule ciclico/settimanale/intervallo via YAML |
| **Modello PD** | Emax in `MathematicalModel`; Hill-Emax in `PharmacodynamicModel` | ✅ Implementato | Hill coefficient solo nel modello esteso (non collegato) |
| **Modello tumorale** | Logistico (principale) / Gompertziano (esteso, non usato) | ✅ Implementato | Principale: logistico. Esteso: Gompertziano ma non connesso |
| **Modello compartimento sano** | Logistico con accoppiamento tossicità | ✅ Implementato | `immune_compromise_index * mean(effects) * H` |
| **Gestione schedule farmaci** | `pharmaco.registry` + preset YAML | ✅ Implementato | Cicli on/off, settimanale, intervallo |
| **Simulazione controfattuale** | Non implementato | ❌ Mancante | Nessun meccanismo per branching da stato reale |
| **Simulazione paziente reale** | Bridge via `twin_assessment_id` | ⚠️ Parziale | Usa parametri twin da un assessment; non consuma storia completa |
| **Simulazione scenario** | `Scenario` + `SimulationAttempt` | ✅ Implementato | Caso d'uso educativo principale |
| **Engine raccomandazioni** | `_interpret_latest_simulation()` + `decision_algorithm.py` | ✅ Implementato | Soglie euristiche, non paziente-specifiche |
| **Generazione artifact** | CSV, HTML plot, JSON twin params per run | ✅ Implementato | `media/sim_data` e `media/sim_plots` |
| **Provenance** | Seed, timestamp, user, parameters per attempt | ⚠️ Parziale | No versione modello, no versione config, no catena formale |
| **Modello di osservazione** | Non implementato | ❌ Mancante | Nessun mapping stato latente → biomarker predetti |
| **Calibrazione** | Non implementato | ❌ Mancante | Nessun fitting parametri su storia paziente |
| **State estimation** | Non implementato | ❌ Mancante | No Kalman, no particle filter, no calibrazione bayesiana |
| **Confronto reale-vs-predetto** | Non implementato | ❌ Mancante | Nessun codice confronta traiettoria simulata con assessment successivi |

---

## Current vs Target Twin Assessment

### Cosa il sistema attuale può già fare

1. Memorizzare demografiche paziente reale, assessment lab ripetuti, citogenetica, timeline terapeutica, sintomi in uno schema relazionale ben strutturato
2. Mappare i biomarker di un singolo assessment in parametri simulatore via risk-weighting configurabile
3. Eseguire simulazione ODE accoppiata PK/PD con scheduling farmaci realistico
4. Generare KPI summary (tumor reduction, healthy loss, TTR, durability, nadir) e salvare dati traiettoria
5. Eseguire replicazione stocastica per coorte con bande di incertezza
6. Presentare raccomandazioni euristiche e stime prognostiche
7. Visualizzare grafici longitudinali biomarker con overlay terapeutici nella patient detail view
8. Iniettare opzionalmente parametri twin-derived in simulazione scenario via `twin_assessment_id`

### Cosa impedisce di essere un vero digital twin paziente-specifico

1. **Nessuna definizione di stato latente.** Il twin è una funzione stateless, non un vettore di stato che evolve nel tempo.
2. **Nessuna persistenza dello stato twin.** Non esiste un modello `PatientTwinState`. L'artifact JSON è un effetto collaterale effimero.
3. **Nessun modello di osservazione.** Non esiste mapping da stato latente del modello (tumor burden, growth rate, concentrazioni farmaco) a biomarker osservati (M-protein, FLC ratio, Hgb, creatinina).
4. **Nessun loop di aggiornamento.** Quando un nuovo assessment viene registrato, nessun codice aggiorna il twin.
5. **Nessuna calibrazione da storia paziente.** Il risk-mapping usa breakpoint globali, non parametri fittati paziente-specifici.
6. **Nessun confronto reale-vs-predetto.** La simulazione produce una traiettoria, ma non viene mai confrontata con assessment reali successivi.
7. **Nessun branching controfattuale.** Non esiste meccanismo per prendere lo stato di un paziente reale al tempo $t$ e simulare terapie alternative.
8. **Input singolo-assessment.** `build_patient_twin` usa un solo assessment; la traiettoria longitudinale completa è ignorata.
9. **Nessun tracking residui.** Non viene calcolato né memorizzato errore tra valori predetti e osservati.
10. **Storia terapeutica disconnessa dalla simulazione.** I record `PatientTherapy` non vengono mai consumati come input di simulazione.

### Layer mancante più critico

**Il modello di osservazione + loop di aggiornamento stato.** Senza modello di osservazione (mapping latente → biomarker osservato), non è possibile confrontare predizioni con realtà. Senza loop di aggiornamento, il twin non apprende mai da nuovi dati. Questi due componenti sono il minimo necessario per la transizione da layer di risk-mapping a vero digital twin.

---

## Architettura target

### Layer 1 — Real Patient

**Scopo:** Ground truth dei dati per follow-up paziente.

**Dati richiesti:** Demografiche paziente, assessment (pannelli lab timestampati), timeline terapeutica (regimen, date, dosi), citogenetica, assessment sintomi.

**Relazione con il codice corrente:** L'app `clinic` fornisce già questo layer con `Patient`, `Assessment`, `PatientTherapy`, `PatientCytogenetics`, `SymptomAssessment`. Nessun cambiamento strutturale necessario.

**Cosa aggiungere:** Un campo `doses` JSONField su `PatientTherapy` (o un modello correlato `PatientTherapyDose`) per registrare le dosi effettive per-farmaco per ciclo, abilitando il consumo diretto da parte del simulatore.

```python
# Esempio estensione PatientTherapy
class PatientTherapy(models.Model):
    # ... campi esistenti ...
    doses = models.JSONField(
        default=dict, blank=True,
        help_text='{"lenalidomide": 25.0, "bortezomib": 1.3, "dexamethasone": 40.0}'
    )
    cycle_length_days = models.IntegerField(null=True, blank=True)
    days_on = models.IntegerField(null=True, blank=True)
```

---

### Layer 2 — Twin State

**Scopo:** Stato biologico paziente-specifico persistente e con versioning che evolve nel tempo.

**Dati richiesti:**
- Vettore di stato latente (e.g., tumor burden stimato, growth rate, competenza immunitaria)
- Stime parametriche (modificatori PK/PD)
- Confidenza/incertezza su ogni parametro
- Provenance della derivazione (quale assessment, quale versione config)

**Relazione con il codice corrente:** Estende l'output di `build_patient_twin()` in un modello persistente. La funzione attuale diventa l'*inizializzatore* per il primo stato twin.

**Cosa aggiungere:** Modello Django `PatientTwinState` e modello `TwinUpdateEvent` per registrare ogni transizione di stato. Un servizio `TwinStateManager` per creare, leggere, aggiornare stati.

---

### Layer 3 — Observation

**Scopo:** Mappare dallo stato latente del modello a valori predetti di biomarker osservabili, e calcolare residui.

**Dati richiesti:**
- Funzione di osservazione $h(\mathbf{x}) \to \hat{\mathbf{y}}$
- Valori osservati reali $\mathbf{y}_{obs}$
- Residuo $\mathbf{r} = \mathbf{y}_{obs} - \hat{\mathbf{y}}$

**Relazione con il codice corrente:** Attualmente assente. Il simulatore produce curve di traiettoria ma non le mappa mai su scale di biomarker osservabili.

**Cosa aggiungere:**
- `observation_model.py`: funzioni di mapping stato modello → biomarker predetti
- Modello Django `ObservationResidual`

Mappature minime iniziali:

| Stato latente | Biomarker predetto | Relazione |
|--------------|-------------------|-----------|
| `tumor_cells` | M-protein (g/dL) | $\hat{M} = \alpha \cdot T + \beta$ (lineare) o $\hat{M} = \alpha \cdot \log_{10}(T) + \beta$ (log) |
| `tumor_cells` | FLC ratio | $\hat{R}_{FLC} = \gamma \cdot (T/T_{ref})^{\delta}$ |
| `healthy_cells` | Hemoglobin (g/dL) | $\hat{Hgb} = Hgb_{baseline} \cdot (H/H_{baseline})^{\eta}$ |
| `healthy_cells` | Albumin (g/dL) | Correlazione indiretta |

---

### Layer 4 — Simulation

**Scopo:** Eseguire simulazioni di traiettoria forward guidate da input terapeutici.

**Dati richiesti:** Stato twin (condizioni iniziali), schedule terapeutico (da `PatientTherapy` reale o regime ipotetico), parametri PK/PD, configurazione modello.

**Relazione con codice corrente:** `MathematicalModel.simulate()` e `SimulationAttempt.run_model()` forniscono già questo. Il cambiamento primario è: accettare condizioni iniziali dallo stato twin + accettare schedule farmaci dalla timeline terapeutica reale.

**Cosa aggiungere:** Una funzione servizio che costruisce un `MathematicalModel` da un `PatientTwinState` + sequenza `PatientTherapy`, senza richiedere un `Scenario`. L'attuale logica `_resolve_solver_inputs` viene estratta in un servizio riutilizzabile.

---

### Layer 5 — Counterfactual

**Scopo:** Branching dallo stato twin calibrato di un paziente e simulazione di regimi terapeutici alternativi.

**Dati richiesti:** Stato twin base al punto di branching, definizione regime alternativo, risultati simulazione per confronto.

**Relazione con codice corrente:** Nessuna implementazione attuale. Costruito sopra Layer 4 eseguendo multiple simulazioni dallo stesso stato iniziale con terapie diverse.

**Cosa aggiungere:** Modello `CounterfactualRun` che collega paziente, stato twin di branching, regime alternativo, e output simulazione risultante. Una view per confronto side-by-side.

---

### Layer 6 — Decision Support

**Scopo:** Interpretare risultati simulazione e stato twin per generare raccomandazioni paziente-specifiche.

**Dati richiesti:** Summary simulazione, residui osservazione (qualità calibrazione), confronti controfattuali, profilo rischio paziente.

**Relazione con codice corrente:** `_interpret_latest_simulation()` e `simulator/decision_algorithm.py` forniscono la base euristica. Vanno estesi con etichette di provenance (e.g., "basato su twin calibrato" vs "basato su mapping singolo assessment").

**Cosa aggiungere:** Labeling confidenza sulle raccomandazioni. Distinzione tra output euristico e model-based. Soglie di qualità-residui che gatano la confidenza delle raccomandazioni.

---

### Layer 7 — Reporting & Provenance

**Scopo:** Tracciabilità completa da dati paziente attraverso derivazione twin, simulazione, fino a raccomandazione.

**Dati richiesti:** Versione modello, versione config twin, versioni preset farmaci, hash commit del codice, metadata run.

**Relazione con codice corrente:** `SimulationAttempt` memorizza parameters, results, artifacts, seed. Necessita formalizzazione con `SimulationRunMetadata`, `ModelVersion`, e `ProvenanceArtifact`.

**Cosa aggiungere:** Snapshot di configurazione versionati. Un modello `ModelVersion`. Una catena di provenance FK da `SimulationAttempt` → `PatientTwinState` → `Assessment`.

---

## Minimal Viable Real-Patient Digital Twin

Il twin minimo viable che può essere costruito sul repo esistente:

### Stati latenti minimi

| Stato | Simbolo | Unità | Origine |
|-------|---------|-------|---------|
| Tumor burden stimato al tempo $t$ | $T_0$ | cellule | Inizializzato da `baseline_tumor_cells` via twin mapping |
| Growth rate tumorale paziente-specifico | $r_T$ | 1/giorno | Inizializzato da `tumor_growth_rate` via twin mapping, poi calibrato |
| Indice compromissione immunitaria | $\sigma_{immune}$ | adimensionale | Inizializzato da `immune_compromise_index` via twin mapping, poi calibrato |

### Biomarker osservati minimi

| Biomarker | Unità | Campo Assessment | Proxy per |
|-----------|-------|-----------------|-----------|
| M-protein | g/dL | `m_protein_g_dl` | Tumor burden (relazione log-lineare approssimata) |
| FLC ratio | κ/λ | `flc_ratio` | Attività tumorale secondaria |
| Hemoglobin | g/dL | `hemoglobin_g_dl` | Funzione midollare sana / tossicità |

### Feature terapeutiche minime

| Feature | Fonte | Note |
|---------|-------|------|
| Identità farmaco | `PatientTherapy.regimen.components` | Parse dei componenti |
| Dose giornaliera approssimata | Nuovo campo `PatientTherapy.doses` o inferenza da `Regimen.components` | JSON struct |
| Date inizio/fine terapia | `PatientTherapy.start_date/end_date` | Già disponibili |

### Loop di aggiornamento minimo

```
1. INIT: Al primo assessment, inizializza stato twin da build_patient_twin()
         → Persisti PatientTwinState(is_current=True)

2. FORWARD: Ad ogni nuovo assessment successivo:
   a) Carica stato twin corrente (is_current=True)
   b) Costruisci schedule farmaci da PatientTherapy tra ultimo e nuovo assessment
   c) Esegui simulazione forward dallo stato twin precedente
      fino alla data del nuovo assessment
   d) Calcola biomarker predetti via observation model:
      - M-protein_pred = h_M(T_simulated)
      - Hgb_pred = h_H(H_simulated)
   e) Calcola residui vs valori reali del nuovo assessment:
      - r_M = M_protein_obs - M_protein_pred
      - r_Hgb = Hgb_obs - Hgb_pred
   f) Calibra: ottimizza r_T e σ_immune per minimizzare |r|
      (scipy.optimize.minimize, gradient-free)
   g) Persisti:
      - Nuovo PatientTwinState(is_current=True)
      - TwinUpdateEvent con residui prima/dopo
      - ObservationResidual con predetti/osservati
   h) Setta stato precedente: is_current=False
```

### Output simulazione minimi

- Traiettoria forward: cellule tumorali, cellule sane, concentrazioni farmaco nel tempo
- Biomarker predetti alle date degli assessment
- Residui (predetto − osservato)
- KPI summary (tumor reduction, healthy loss, TTR)

---

## Modelli Django proposti

### `PatientTwinState`

**Scopo:** Stato twin persistente e versionato per ogni paziente in un momento temporale.

```python
class PatientTwinState(models.Model):
    """Stato twin persistente paziente-specifico."""
    
    patient = models.ForeignKey(
        "clinic.Patient", on_delete=models.CASCADE,
        related_name="twin_states"
    )
    assessment = models.ForeignKey(
        "clinic.Assessment", on_delete=models.SET_NULL,
        null=True, blank=True,
        help_text="Assessment da cui questo stato è stato derivato"
    )
    created_at = models.DateTimeField(auto_now_add=True)
    
    # Vettore di stato
    state_vector = models.JSONField(
        help_text='{"tumor_burden": float, "growth_rate": float, '
                  '"immune_index": float, "carrying_capacity_tumor": float, ...}'
    )
    
    # Incertezza per-parametro
    parameter_confidence = models.JSONField(
        default=dict, blank=True,
        help_text="Stime incertezza per ogni parametro"
    )
    
    risk_score = models.FloatField()
    
    twin_config_version = models.CharField(
        max_length=64,
        help_text="Hash/tag di twin_risk.yaml al momento della derivazione"
    )
    
    is_current = models.BooleanField(
        default=True,
        help_text="Marca lo stato più recente per questo paziente"
    )
    
    method = models.CharField(
        max_length=32, default="initial_mapping",
        help_text="Metodo di derivazione: initial_mapping, "
                  "residual_minimization, manual_override"
    )
    
    class Meta:
        ordering = ["-created_at"]
        indexes = [
            models.Index(fields=["patient", "-created_at"]),
            models.Index(fields=["patient", "is_current"]),
        ]
```

### `TwinUpdateEvent`

**Scopo:** Registrare ogni transizione di stato twin con provenance.

```python
class TwinUpdateEvent(models.Model):
    """Log di aggiornamento stato twin."""
    
    patient = models.ForeignKey(
        "clinic.Patient", on_delete=models.CASCADE,
        related_name="twin_updates"
    )
    previous_state = models.ForeignKey(
        PatientTwinState, on_delete=models.SET_NULL,
        null=True, blank=True, related_name="outgoing_updates"
    )
    new_state = models.ForeignKey(
        PatientTwinState, on_delete=models.CASCADE,
        related_name="incoming_update"
    )
    trigger_assessment = models.ForeignKey(
        "clinic.Assessment", on_delete=models.SET_NULL,
        null=True, blank=True
    )
    
    method = models.CharField(max_length=32)
    residuals_before = models.JSONField(default=dict, blank=True)
    residuals_after = models.JSONField(default=dict, blank=True)
    
    calibration_iterations = models.IntegerField(null=True, blank=True)
    calibration_converged = models.BooleanField(null=True, blank=True)
    
    created_at = models.DateTimeField(auto_now_add=True)
    
    class Meta:
        ordering = ["-created_at"]
```

### `ObservationResidual`

**Scopo:** Registrare biomarker predetti vs osservati per ogni assessment.

```python
class ObservationResidual(models.Model):
    """Confronto predetto vs osservato per assessment."""
    
    twin_state = models.ForeignKey(
        PatientTwinState, on_delete=models.CASCADE,
        related_name="residuals"
    )
    assessment = models.ForeignKey(
        "clinic.Assessment", on_delete=models.CASCADE,
        related_name="twin_residuals"
    )
    
    predicted_values = models.JSONField(
        help_text='{"m_protein": float, "flc_ratio": float, "hemoglobin": float}'
    )
    observed_values = models.JSONField(
        help_text="Stessa struttura, valori reali dall'assessment"
    )
    residuals = models.JSONField(
        help_text="Differenza per-biomarker"
    )
    
    rmse = models.FloatField(help_text="Residuo aggregato")
    
    created_at = models.DateTimeField(auto_now_add=True)
    
    class Meta:
        ordering = ["-created_at"]
        unique_together = ("twin_state", "assessment")
```

### `CounterfactualRun`

**Scopo:** Registrare simulazioni terapeutiche alternative contro uno stato twin calibrato.

```python
class CounterfactualRun(models.Model):
    """Simulazione 'what-if' da stato twin calibrato."""
    
    patient = models.ForeignKey(
        "clinic.Patient", on_delete=models.CASCADE,
        related_name="counterfactual_runs"
    )
    base_twin_state = models.ForeignKey(
        PatientTwinState, on_delete=models.CASCADE,
        related_name="counterfactual_runs"
    )
    actual_therapy = models.ForeignKey(
        "clinic.PatientTherapy", on_delete=models.SET_NULL,
        null=True, blank=True
    )
    alternative_regimen = models.ForeignKey(
        "clinic.Regimen", on_delete=models.CASCADE
    )
    alternative_parameters = models.JSONField(
        default=dict, blank=True,
        help_text="Override dosi per il regime alternativo"
    )
    
    simulation_result = models.JSONField(help_text="KPI summary")
    trajectory_csv = models.FileField(
        upload_to="counterfactual/", null=True, blank=True
    )
    comparison_metrics = models.JSONField(
        default=dict, blank=True,
        help_text="Delta metriche vs terapia attuale"
    )
    
    created_at = models.DateTimeField(auto_now_add=True)
    user = models.ForeignKey(
        "auth.User", on_delete=models.SET_NULL,
        null=True, blank=True
    )
    
    class Meta:
        ordering = ["-created_at"]
```

### `SimulationRunMetadata`

**Scopo:** Provenance formale per ogni esecuzione di simulazione.

```python
class SimulationRunMetadata(models.Model):
    """Metadata di provenance per simulazione."""
    
    attempt = models.OneToOneField(
        "simulator.SimulationAttempt",
        on_delete=models.CASCADE,
        related_name="run_metadata"
    )
    
    model_version = models.CharField(max_length=64)
    twin_config_hash = models.CharField(
        max_length=64,
        help_text="SHA256 di twin_risk.yaml al momento della run"
    )
    drug_preset_hashes = models.JSONField(
        default=dict,
        help_text="SHA256 per-drug preset YAML"
    )
    solver_method = models.CharField(max_length=32, default="RK45")
    solver_rtol = models.FloatField(default=1e-6)
    solver_atol = models.FloatField(default=1e-8)
    
    execution_time_ms = models.IntegerField(null=True, blank=True)
    
    twin_state = models.ForeignKey(
        PatientTwinState, on_delete=models.SET_NULL,
        null=True, blank=True
    )
    
    code_commit_hash = models.CharField(
        max_length=40, blank=True,
        help_text="Git commit hash del codebase"
    )
    
    created_at = models.DateTimeField(auto_now_add=True)
```

### `ModelVersion`

**Scopo:** Tracciare versioni del modello matematico e preset di configurazione.

```python
class ModelVersion(models.Model):
    """Versione congelata del modello e configurazioni."""
    
    version_tag = models.CharField(max_length=32, unique=True)
    description = models.TextField(blank=True)
    
    twin_config_snapshot = models.JSONField(
        help_text="Copia congelata di twin_risk.yaml"
    )
    drug_presets_snapshot = models.JSONField(
        help_text="Copia congelata di tutti i drug YAML"
    )
    model_equations_hash = models.CharField(
        max_length=64,
        help_text="Hash di models_simulation.py"
    )
    
    is_active = models.BooleanField(default=True)
    created_at = models.DateTimeField(auto_now_add=True)
    
    class Meta:
        ordering = ["-created_at"]
    
    def __str__(self):
        return self.version_tag
```

---

## Architettura moduli/servizi proposta

```
twin_engine/
    __init__.py
    state_model.py          # Definizione stato latente + inizializzazione
    observation_model.py    # Mapping stato → biomarker osservabili
    updater.py              # Loop aggiornamento stato
    calibration.py          # Fitting parametri su dati osservati
    therapy_schedule.py     # Conversione PatientTherapy → dose_functions
    counterfactual.py       # Branching e confronto scenari
    provenance.py           # Catena provenance e versionamento
    report_builder.py       # Report paziente-specifico
```

### `state_model.py`

**Responsabilità:**
- Definisce dataclass `TwinState` che rappresenta il vettore di stato latente
- Fornisce `initialize_from_assessment(assessment: Assessment) -> TwinState`: wrappa l'attuale `build_patient_twin()`, aggiunge persistenza via `PatientTwinState`
- Serializzazione/deserializzazione da/verso modello `PatientTwinState`
- Gestione flag `is_current`

```python
@dataclass
class TwinState:
    tumor_burden: float          # cellule stimate
    growth_rate: float           # r_T (1/giorno)
    immune_index: float          # σ_immune
    carrying_capacity_tumor: float
    carrying_capacity_healthy: float
    healthy_growth_rate: float
    risk_score: float
    confidence: dict[str, float]  # incertezza per-parametro

def initialize_from_assessment(assessment: Assessment) -> TwinState:
    """Wrappa build_patient_twin() con stato strutturato."""
    raw = build_patient_twin(assessment)
    return TwinState(
        tumor_burden=raw.get("carrying_capacity_tumor", 1e9) * 0.1,
        growth_rate=raw["tumor_growth_rate"],
        immune_index=raw["immune_compromise_index"],
        carrying_capacity_tumor=raw["carrying_capacity_tumor"],
        carrying_capacity_healthy=raw["carrying_capacity_healthy"],
        healthy_growth_rate=raw["healthy_growth_rate"],
        risk_score=raw["risk_score"],
        confidence={k: 0.5 for k in raw},  # inizialmente bassa
    )
```

### `observation_model.py`

**Responsabilità:**
- Implementa funzioni di osservazione che mappano stato latente → biomarker predetti
- Genera residui
- Configurabile con parametri di calibrazione

```python
def predict_m_protein(
    tumor_burden: float,
    calibration_alpha: float = 1e-10,
    calibration_beta: float = 0.1,
) -> float:
    """M-protein ∝ tumor burden (relazione log-lineare)."""
    if tumor_burden <= 0:
        return 0.0
    return calibration_alpha * tumor_burden + calibration_beta

def predict_hemoglobin(
    healthy_cells: float,
    baseline_healthy: float = 5e11,
    baseline_hgb: float = 13.0,
) -> float:
    """Hemoglobin ∝ (H/H_baseline)^η."""
    ratio = max(healthy_cells, 0) / max(baseline_healthy, 1e-6)
    return baseline_hgb * (ratio ** 0.3)

def compute_residuals(
    predicted: dict[str, float],
    observed: dict[str, float],
) -> dict[str, float]:
    """Residui per-biomarker."""
    residuals = {}
    for key in predicted:
        if key in observed and observed[key] is not None:
            residuals[key] = observed[key] - predicted[key]
    return residuals
```

### `updater.py`

**Responsabilità:**
- Implementa il loop di aggiornamento stato:
  1. Carica stato twin corrente
  2. Esegue simulazione forward fino a nuovo assessment
  3. Chiama observation model per predire biomarker
  4. Chiama calibration per aggiustare parametri
  5. Persiste nuovo `PatientTwinState`, `TwinUpdateEvent`, `ObservationResidual`

```python
def update_twin(
    patient: Patient,
    new_assessment: Assessment,
    *,
    current_state: PatientTwinState | None = None,
) -> tuple[PatientTwinState, ObservationResidual]:
    """Aggiorna stato twin dato un nuovo assessment."""
    if current_state is None:
        current_state = PatientTwinState.objects.filter(
            patient=patient, is_current=True
        ).first()
    
    if current_state is None:
        # Primo stato: inizializza da assessment
        twin_state = initialize_from_assessment(new_assessment)
        # ... persisti e ritorna
    
    # Forward simulation dal tempo dello stato precedente
    # al tempo del nuovo assessment
    # ... costruisci MathematicalModel, simula, predici, calibra, persisti
```

### `calibration.py`

**Responsabilità:**
- Fitting parametri contro dati osservati
- Inizialmente: `scipy.optimize.minimize` gradient-free su $r_T$ e $\sigma_{immune}$
- Futuro: interfaccia per Kalman filter, particle filter, MCMC

```python
def calibrate(
    twin_state: TwinState,
    observation_history: list[dict],
    therapy_history: list[PatientTherapy],
) -> TwinState:
    """Calibra parametri twin per minimizzare residui."""
    from scipy.optimize import minimize
    
    def objective(params):
        r_T, sigma = params
        # Simula forward con r_T e sigma
        # Calcola residui vs observation_history
        # Restituisci RMSE totale
        ...
    
    result = minimize(
        objective,
        x0=[twin_state.growth_rate, twin_state.immune_index],
        method="Nelder-Mead",
        options={"maxiter": 100},
    )
    
    # Aggiorna e restituisci
    ...
```

### `therapy_schedule.py`

**Responsabilità:**
- Converte queryset `PatientTherapy` → `dose_functions: Dict[str, Callable[[float], float]]` per l'ODE solver
- Gestisce gap inter-linea, periodi mantenimento, riduzioni dose
- Interfaccia con il registry farmaci per recuperare schedule type

```python
def build_dose_functions(
    therapies: QuerySet[PatientTherapy],
    reference_date: date,
) -> dict[str, Callable[[float], float]]:
    """Converte timeline terapeutica reale in funzioni dose per ODE."""
    drug_schedule: dict[str, list[tuple[float, float, float]]] = {}
    # Per ogni terapia: calcola intervallo temporale relativo,
    # identifica farmaci dal regimen, assegna dosi
    ...
```

### `counterfactual.py`

**Responsabilità:**
- Branching da stato twin calibrato con regime alternativo
- Confronto metriche (delta tumor reduction, delta healthy loss, delta TTR)

```python
def run_counterfactual(
    twin_state: PatientTwinState,
    alternative_regimen: Regimen,
    time_horizon: float = 180.0,
    alternative_doses: dict | None = None,
) -> CounterfactualResult:
    """Esegui simulazione controfattuale."""
    # Costruisci MathematicalModel dallo stato twin + regime alternativo
    # Simula, genera summary, confronta con terapia attuale
    ...
```

### `provenance.py`

**Responsabilità:**
- Costruzione catena provenance
- Hash di configurazioni
- Snapshot versione modello

```python
def record_run_provenance(
    attempt: SimulationAttempt,
    twin_state: PatientTwinState | None,
) -> SimulationRunMetadata:
    """Registra metadata provenance per una run."""
    import hashlib
    from pathlib import Path
    
    twin_config_path = Path(__file__).parent.parent / "simulator" / "presets" / "twin_risk.yaml"
    config_hash = hashlib.sha256(
        twin_config_path.read_bytes()
    ).hexdigest()[:16]
    
    # ... crea e restituisci SimulationRunMetadata
```

### `report_builder.py`

**Responsabilità:**
- Generazione report twin paziente-specifico
- Aggregazione storia twin, qualità calibrazione, predizioni, confronti controfattuali

```python
def build_twin_report(patient: Patient) -> dict:
    """Aggrega storia twin, calibrazione, predizioni, controfattuali."""
    states = PatientTwinState.objects.filter(patient=patient).order_by("created_at")
    residuals = ObservationResidual.objects.filter(
        twin_state__patient=patient
    ).order_by("assessment__date")
    counterfactuals = CounterfactualRun.objects.filter(
        patient=patient
    ).order_by("-created_at")
    
    return {
        "patient_id": patient.pk,
        "twin_history": [...],
        "calibration_quality": {...},
        "latest_predictions": {...},
        "counterfactual_comparisons": [...],
    }
```

---

## Piano di integrazione con il repo esistente

### Codice da riutilizzare senza modifiche

| Componente | File | Motivazione |
|------------|------|-------------|
| Modelli paziente | `clinic/models.py` (`Patient`, `Assessment`, `PatientTherapy`, `PatientCytogenetics`, `SymptomAssessment`) | Layer dati paziente completo e ben strutturato |
| ODE solver | `simulator/models_simulation.py` (`MathematicalModel.simulate()`) | Engine di simulazione funzionante e testato |
| Registry farmaci | `simulator/pharmaco/registry.py` | Caricamento preset e schedule |
| Preset YAML | `simulator/presets/` | Configurazione farmaci e twin |
| Decision algorithm | `simulator/decision_algorithm.py` | Definizioni trasparenza (augmentare, non sostituire) |
| Prognosis | `simulator/prognosis.py` | Stime sopravvivenza |
| View CRUD clinic | `clinic/views.py` (dashboard, patient_list, patient_new, assessment_new, etc.) | Flussi UI funzionanti |
| Regimen suggester | `simulator/regimen_suggester.py` | Database regimi curato |
| Virtual patients | `simulator/virtual_patients.py` | Generazione trial in-silico |
| Modulo design | `simulator/design/` | Simulazione modalità avanzate |

### Codice da wrappare (adapter di interfaccia)

| Componente | Azione | Nuovo wrapper |
|------------|--------|---------------|
| `build_patient_twin()` | Wrappare in `twin_engine.state_model.initialize_from_assessment()`. Funzione esistente diventa implementazione interna; il wrapper aggiunge persistenza via `PatientTwinState` | `twin_engine/state_model.py` |
| `_resolve_solver_inputs()` | Estrarre in `twin_engine.simulation_bridge.build_solver_inputs(twin_state, therapy_schedule)` riutilizzabile senza `Scenario` | `twin_engine/simulation_bridge.py` |
| `_summarize_trajectory()` | Estrarre in `twin_engine.report_builder.summarize_trajectory()` | `twin_engine/report_builder.py` |

### Codice da refactorare

| Componente | File | Motivazione | Azione |
|------------|------|-------------|--------|
| `SimulationAttempt.run_model()` | `simulator/models.py:335` | >200 righe in un singolo metodo | Decomporre in: risoluzione parametri, derivazione twin, costruzione modello, esecuzione simulazione, calcolo summary, salvataggio artifact |
| `SimulationAttempt.parameters` | `simulator/models.py:310` | JSONField non tipizzato | Aggiungere validazione schema (JSON Schema o Pydantic) per distinguere run scenario-driven da real-patient-driven |
| `_interpret_latest_simulation()` | `clinic/views.py:21` | Non distingue fonte | Refactorare per accettare parametro `source` ("heuristic" vs "calibrated_twin") e aggiustare etichette confidenza |
| `SimulationAttempt.scenario` FK | `simulator/models.py:285` | FK non-nullable → richiede Scenario per ogni run | Rendere nullable per supportare run paziente-driven senza Scenario |

### Naming fuorviante da chiarire

| Termine corrente | Cosa significa realmente | Azione |
|------------------|------------------------|--------|
| "Patient Twin" nel README, help text | Single-assessment risk-mapping | Distinguere "Patient Risk Profile" (attuale) da "Patient Digital Twin" (target) |
| `build_patient_twin()` | Derivazione risk-profile | Opzione: rinominare a `derive_risk_profile()` o mantenere con docstring chiara |
| `twin_assessment_id` parametro | Assessment usato per risk-profile derivation | Documentare che NON è una sorgente di calibrazione twin |
| "Digital Patient Twins — Virtual patient modeling" nel README | Claim aspirazionale | Ammorbidire a "Patient Risk Profiling" fino a operatività twin engine |

---

## Backlog GitHub Issues

### Issue 1 — Separare flusso simulazione paziente reale da scenario educativo

| Campo | Valore |
|-------|--------|
| **Titolo** | `[Architecture] Creare path simulazione paziente-driven indipendente da Scenario` |
| **Goal** | Abilitare esecuzione simulazione dallo stato twin + storia terapeutica di un paziente senza richiedere un Scenario shell |
| **Perché è importante** | Tutte le simulazioni passano attualmente per `SimulationAttempt → Scenario`. Le simulazioni paziente reale sono forzate in questo scaffolding educativo via `twin_assessment_id`, conflando due workflow distinti |
| **File rilevanti** | `simulator/models.py`, `simulator/views_manage.py`, `clinic/views.py` |
| **Scope implementazione** | (1) Estrarre `_resolve_solver_inputs()` e `_summarize_trajectory()` da `run_model()` in servizio standalone. (2) Rendere `SimulationAttempt.scenario` FK nullable per supportare run scenario-free. (3) Aggiungere pulsante "Run Simulation" nella `patient_detail` che costruisce solver inputs da `PatientTwinState` + `PatientTherapy`. (4) Mantenere flusso scenario invariato per uso educativo. |
| **Acceptance Criteria** | Una simulazione può essere lanciata da `patient_detail` senza selezionare uno Scenario. La simulazione usa la timeline terapeutica reale del paziente come input schedule. I risultati sono salvati con provenance che linka al paziente. |
| **Priorità** | **P0** |

---

### Issue 2 — Persistere stato twin paziente

| Campo | Valore |
|-------|--------|
| **Titolo** | `[Twin] Aggiungere modello PatientTwinState e layer persistenza` |
| **Goal** | Sostituire computazione twin effimera con stato persistente e versionato per paziente |
| **Perché è importante** | Il twin attuale è ricalcolato da zero ad ogni run. Nessuno storico stato twin è disponibile per confronto, analisi trend, o riferimento calibrazione |
| **File rilevanti** | `simulator/twin.py`, `simulator/models.py`, `simulator/presets/twin_risk.yaml` |
| **Scope implementazione** | (1) Creare modello Django `PatientTwinState` (vedi spec sopra). (2) Modificare `build_patient_twin()` per persistere anche il risultato. (3) Gestione flag `is_current` (set precedente False, nuovo True). (4) Management command `update_all_twin_states` per popolazione batch da assessment esistenti. (5) Vista storico stato twin nella patient detail. |
| **Acceptance Criteria** | Un record `PatientTwinState` è creato per ogni derivazione twin da assessment. `is_current` marca correttamente uno e un solo stato per paziente. Lo storico è interrogabile via ORM. |
| **Priorità** | **P0** |

---

### Issue 3 — Aggiornare twin da nuovi assessment paziente

| Campo | Valore |
|-------|--------|
| **Titolo** | `[Twin] Implementare loop aggiornamento stato twin su creazione nuovo assessment` |
| **Goal** | Aggiornare automaticamente lo stato twin quando un nuovo assessment viene registrato |
| **Perché è importante** | Fondamento di un digital twin dinamico: il sistema deve reagire a nuove osservazioni cliniche |
| **File rilevanti** | `clinic/views.py` (`assessment_new`), `simulator/twin.py`, nuovo `twin_engine/updater.py` |
| **Scope implementazione** | (1) Creare `twin_engine/updater.py` con `update_twin(patient, new_assessment)`. (2) Nel view `assessment_new` (POST success path), chiamare `update_twin()`. (3) Persistere `TwinUpdateEvent`. (4) Inizialmente: ri-derivare da nuovo assessment via `build_patient_twin()` (update stateless). (5) Successivamente: simulazione forward + residui (dipendenza Issue #5). |
| **Acceptance Criteria** | Creare un nuovo assessment triggera aggiornamento stato twin. Un record `TwinUpdateEvent` è creato con `method="initial_mapping"`. Lo stato precedente `is_current` è settato a False. |
| **Priorità** | **P0** |

---

### Issue 4 — Integrare timeline PatientTherapy come input simulazione

| Campo | Valore |
|-------|--------|
| **Titolo** | `[Simulation] Consumare storia terapeutica reale come input schedule farmaci` |
| **Goal** | Usare record `PatientTherapy` reali (regimen, date, dosi) per costruire funzioni schedule dose per la simulazione |
| **Perché è importante** | Le dosi farmaci attualmente sono inserite manualmente nel form simulazione o derivate da preset. La storia terapeutica reale è memorizzata ma mai usata come input di simulazione |
| **File rilevanti** | `clinic/models.py` (`PatientTherapy`, `Regimen`), `simulator/pharmaco/registry.py`, `simulator/models_simulation.py` |
| **Scope implementazione** | (1) Aggiungere campo `doses` JSONField a `PatientTherapy` per registrare dosi per-farmaco. (2) Creare `twin_engine/therapy_schedule.py` che converte queryset `PatientTherapy` → `dose_functions: Dict[str, Callable]` per ODE solver. (3) Gestire gap inter-linea, mantenimento, riduzioni dose. (4) UI per review/edit dello schedule inferito prima della simulazione. |
| **Acceptance Criteria** | La storia terapeutica di un paziente può essere automaticamente convertita in funzioni dose solver-compatible. Terapia multi-regimen sequenziale è supportata. Lo schedule generato è previsualizzato nella UI. |
| **Priorità** | **P1** |

---

### Issue 5 — Aggiungere confronto traiettoria reale vs simulata

| Campo | Valore |
|-------|--------|
| **Titolo** | `[Twin] Implementare observation model e tracking residui` |
| **Goal** | Predire valori biomarker dall'output simulazione e confrontare con valori assessment reali |
| **Perché è importante** | Loop di feedback critico che abilita calibrazione. Senza confronto predetto-vs-osservato, il twin non può apprendere |
| **File rilevanti** | Nuovo modulo `twin_engine/observation_model.py`, nuovo modello `ObservationResidual` |
| **Scope implementazione** | (1) Implementare funzioni osservazione: tumor_cells → M-protein predetto, healthy_cells → Hgb predetto. (2) Creare modello `ObservationResidual`. (3) Dopo simulazione forward a data assessment, chiamare observation model, calcolare residui, persistere. (4) Display residui nella patient detail (overlay grafico: predetto vs osservato). |
| **Acceptance Criteria** | Per ogni assessment, valori biomarker predetti e osservati sono memorizzati. RMSE per assessment è calcolato e memorizzato. Residui sono visualizzati nella UI paziente. |
| **Priorità** | **P1** |

---

### Issue 6 — Aggiungere branching controfattuale per pazienti reali

| Campo | Valore |
|-------|--------|
| **Titolo** | `[Twin] Implementare simulazione controfattuale da stato twin calibrato` |
| **Goal** | Abilitare simulazioni "what-if" con regimi alternativi dallo stato twin corrente di un paziente reale |
| **Perché è importante** | Il ragionamento controfattuale è la proposizione di valore clinico core di un digital twin — rispondere a "cosa succederebbe se cambiassimo terapia?" |
| **File rilevanti** | Nuovo `twin_engine/counterfactual.py`, nuovo modello `CounterfactualRun` |
| **Scope implementazione** | (1) Creare modello `CounterfactualRun`. (2) Implementare `run_counterfactual(twin_state, alternative_regimen)`. (3) UI nella patient detail: dropdown "Simula Regime Alternativo" → run → display confronto side-by-side. (4) Memorizzare tutte le run con provenance completa. |
| **Acceptance Criteria** | Un utente può selezionare un regime alternativo per un paziente e lanciare simulazione controfattuale. Risultati sono visualizzati accanto alla traiettoria terapia attuale. Tutte le run sono persistite con provenance. |
| **Priorità** | **P1** |

---

### Issue 7 — Aggiungere provenance e tracciabilità per run twin/simulazione

| Campo | Valore |
|-------|--------|
| **Titolo** | `[Provenance] Aggiungere SimulationRunMetadata e ModelVersion` |
| **Goal** | Tracciabilità formale di ogni esecuzione di simulazione e del modello usato |
| **Perché è importante** | Riproducibilità scientifica richiede sapere esattamente quale modello, configurazione, e dati hanno prodotto ogni risultato |
| **File rilevanti** | `simulator/models.py`, `simulator/presets/`, nuovo `twin_engine/provenance.py` |
| **Scope implementazione** | (1) Creare modelli `SimulationRunMetadata` e `ModelVersion`. (2) Hash automatico di `twin_risk.yaml` e drug preset ad ogni run. (3) FK da metadata → `PatientTwinState` quando applicabile. (4) Endpoint API per query provenance. |
| **Acceptance Criteria** | Ogni `SimulationAttempt` ha un `SimulationRunMetadata` associato. Hash config sono registrati. Versione modello è tracciabile. |
| **Priorità** | **P1** |

---

### Issue 8 — Aggiungere dashboard paziente-twin

| Campo | Valore |
|-------|--------|
| **Titolo** | `[UI] Dashboard paziente-twin integrata` |
| **Goal** | Pannello unificato che mostra storia twin, qualità calibrazione, predizioni, controfattuali |
| **Perché è importante** | I clinici hanno bisogno di una vista consolidata dello "stato" del paziente secondo il digital twin |
| **File rilevanti** | `clinic/views.py` (`patient_detail`), template `clinic/templates/clinic/patient_detail.html` |
| **Scope implementazione** | (1) Sezione dedicata "Digital Twin" nella patient detail con: (a) timeline stati twin con risk_score trend, (b) grafico residui (predetto vs osservato), (c) ultima predizione forward, (d) lista run controfattuali con confronto metriche, (e) indicatore qualità calibrazione (RMSE trend). (2) Tab separato o accordion per non sovraccaricare la vista esistente. |
| **Acceptance Criteria** | La patient detail mostra una sezione twin con storia, residui, e predizioni. Per pazienti senza stato twin: messaggio "Twin non ancora inizializzato" con CTA. |
| **Priorità** | **P2** |

---

### Issue 9 — Refactoring engine raccomandazioni

| Campo | Valore |
|-------|--------|
| **Titolo** | `[Decision Support] Distinguere output euristico da paziente-specifico` |
| **Goal** | `_interpret_latest_simulation()` deve distinguere raccomandazioni euristiche da quelle basate su twin calibrato |
| **Perché è importante** | La confidenza clinica in una raccomandazione dipende dalla qualità della fonte: un risk-mapping generico vs un twin calibrato su più assessment hanno affidabilità radicalmente diverse |
| **File rilevanti** | `clinic/views.py` (`_interpret_latest_simulation`), `simulator/decision_algorithm.py` |
| **Scope implementazione** | (1) Aggiungere parametro `source` alla funzione di interpretazione. (2) Quando `source="calibrated_twin"`: aggiungere confidenza basata su RMSE residui. (3) Quando `source="heuristic"`: etichettare esplicitamente come "basato su soglie generali, non calibrato". (4) Badge di confidenza nella UI (alto/medio/basso). |
| **Acceptance Criteria** | Ogni raccomandazione ha un'etichetta fonte e confidenza. La UI mostra badge differenti per raccomandazioni euristiche vs calibrate. |
| **Priorità** | **P2** |

---

### Issue 10 — Preparare architettura per state estimation

| Campo | Valore |
|-------|--------|
| **Titolo** | `[Twin] Preparare interfaccia per metodi di stima stato futuri` |
| **Goal** | Definire interfacce astratte in `twin_engine/calibration.py` che permettano plug-in di EKF, particle filter, MCMC senza refactoring |
| **Perché è importante** | Il metodo iniziale (Nelder-Mead) è quick-start; metodi avanzati richiedono la stessa interfaccia ma implementazione diversa |
| **File rilevanti** | Nuovo `twin_engine/calibration.py` |
| **Scope implementazione** | (1) Definire Protocol `CalibrationMethod` con metodo `calibrate(state, observations, therapy_history) -> CalibratedState`. (2) Implementare `NelderMeadCalibrator` (iniziale). (3) Stub per `EKFCalibrator`, `ParticleFilterCalibrator`. (4) Documentare contratto e requisiti per ognuno. |
| **Acceptance Criteria** | L'interfaccia `CalibrationMethod` è definita e `NelderMeadCalibrator` è funzionante. Uno stub per EKF è presente con docstring che specifica cosa serve per implementarlo. |
| **Priorità** | **P2** |

---

### Issue 11 — Aggiornare documentazione "patient twin"

| Campo | Valore |
|-------|--------|
| **Titolo** | `[Docs] Descrivere accuratamente patient twin: stato attuale vs target` |
| **Goal** | La documentazione deve distinguere chiaramente risk-profile mapping (attuale) da digital twin dinamico (target) |
| **Perché è importante** | Naming fuorviante crea aspettative errate e potenziale rischio comunicativo in contesto clinico/regolatorio |
| **File rilevanti** | `README.md`, `docs/guides/architecture.md`, `simulator/twin.py` docstring, template help text |
| **Scope implementazione** | (1) Aggiornare README: "Patient Risk Profiling" anziché "Digital Patient Twins" finché twin engine non è operativo. (2) Aggiornare docstring di `build_patient_twin()`. (3) Creare guida dedicata `docs/guides/current_twin_architecture.md`. (4) Aggiornare help text UI per distinguere "Twin Assessment" da "Twin Calibration". |
| **Acceptance Criteria** | Nessun punto della documentazione o UI claim capacità di digital twin dinamico non ancora implementata. La guida spiega chiaramente cosa è implementato e il roadmap. |
| **Priorità** | **P1** |

---

### Issue 12 — Aggiungere reporting/export per simulazioni follow-up paziente-specifiche

| Campo | Valore |
|-------|--------|
| **Titolo** | `[Reporting] Report twin paziente-specifico esportabile` |
| **Goal** | Generare report HTML/PDF che mostri storia twin, calibrazione, predizioni, controfattuali per un paziente |
| **Perché è importante** | Documentazione clinica richiede output strutturato e archiviabile delle simulazioni paziente-specifiche |
| **File rilevanti** | Nuovo `twin_engine/report_builder.py`, template dedicati |
| **Scope implementazione** | (1) Implementare `build_twin_report(patient)` che aggrega tutto. (2) Template HTML dedicato per visualizzazione e stampa. (3) Export PDF (via weasyprint o simile). (4) Endpoint API per download. |
| **Acceptance Criteria** | Un report twin completo può essere generato per qualsiasi paziente con storico twin. Il report include: grafico traiettoria predetta, tabella residui, confronti controfattuali, metadati provenance. Il formato è stampabile. |
| **Priorità** | **P2** |

---

## Documentation Roadmap

### `docs/guides/current_twin_architecture.md`

| Campo | Valore |
|-------|--------|
| **Perché serve** | La documentazione attuale non distingue chiaramente il risk-mapping dal concetto di digital twin. Serve un documento di riferimento tecnico che descriva cosa è realmente implementato |
| **Contenuto** | (1) Cos'è `build_patient_twin()` e cosa non è. (2) Flusso dati: Assessment → risk_score → parametri simulatore. (3) Configurazione YAML (peso, breakpoint, mapping). (4) Limitazioni esplicite: stateless, singolo assessment, nessuna calibrazione. (5) Diagramma data flow. |
| **File da referenziare** | `simulator/twin.py`, `simulator/presets/twin_risk.yaml`, `simulator/models.py` (run_model, _merge_twin_parameters) |

### `docs/guides/real_patient_followup_workflow.md`

| Campo | Valore |
|-------|--------|
| **Perché serve** | Non esiste documentazione che descriva il flusso completo di follow-up paziente reale end-to-end |
| **Contenuto** | (1) Creazione paziente e primo assessment. (2) Registrazione terapie e assessment successivi. (3) Come i dati paziente reali si collegano alla simulazione (bridge twin_assessment_id). (4) Visualizzazione timeline nella patient_detail. (5) Interpretazione raccomandazioni. (6) Limiti attuali del collegamento paziente→simulazione. |
| **File da referenziare** | `clinic/views.py`, `clinic/models.py`, `simulator/api_twin.py` |

### `docs/guides/patient_twin_target_architecture.md`

**Questo documento** (il presente file). Quando il twin engine sarà implementato, aggiornare per riflettere lo stato corrente piuttosto che il target.

### `docs/guides/counterfactual_simulation.md`

| Campo | Valore |
|-------|--------|
| **Perché serve** | Il ragionamento controfattuale è la killer feature del digital twin e va documentato per clinici e sviluppatori |
| **Contenuto** | (1) Cos'è una simulazione controfattuale. (2) Prerequisiti (stato twin calibrato). (3) Flusso utente: selezione regime alternativo → simulazione → confronto. (4) Interpretazione risultati. (5) Limitazioni (incertezza parametrica, assunzioni modello). (6) Esempi con dati fittizi. |
| **File da referenziare** | `twin_engine/counterfactual.py` (futuro), `simulator/regimen_suggester.py` |

### `docs/guides/provenance_and_traceability.md`

| Campo | Valore |
|-------|--------|
| **Perché serve** | La riproducibilità scientifica e la compliance regolatoria richiedono documentazione della catena di provenance |
| **Contenuto** | (1) Catena provenance: Patient → Assessment → TwinState → Simulation → Recommendation. (2) Cosa viene tracciato ad ogni step. (3) Versionamento: modello, config, preset. (4) Come riprodurre una run data il suo metadata. (5) Schema FK di provenance. |
| **File da referenziare** | `SimulationRunMetadata`, `ModelVersion`, `twin_engine/provenance.py` (futuri) |

### `docs/research/state_estimation_readiness.md`

| Campo | Valore |
|-------|--------|
| **Perché serve** | Documentare la preparazione architetturale per metodi di stima stato avanzati senza confonderli con feature implementate |
| **Contenuto** | (1) Definizione del problema di state estimation nel contesto MM. (2) Metodi candidati: EKF, particle filter, MCMC, Bayesian calibration. (3) Requisiti per-metodo: modello di transizione stato, modello di osservazione lineare/non-lineare, distribuzione a priori. (4) Mapping ai componenti del repo: observer model → `observation_model.py`, transition → `MathematicalModel`, prior → `twin_risk.yaml`. (5) Roadmap di implementazione incrementale. (6) Riferimenti bibliografici chiave. |
| **File da referenziare** | `simulator/models_simulation.py`, `twin_engine/calibration.py` (futuro), `twin_engine/observation_model.py` (futuro) |

---

## Conclusione tecnica finale

### Stato attuale

bmyCure4MM è una piattaforma Django ben strutturata con:
- Layer dati paziente reale **completo e funzionante** (demografiche, assessment IMWG, terapie, citogenetica, sintomi)
- Simulatore ODE **robusto** (PK/PD accoppiato, schedule configurabili, coorte stocastica, ottimizzazione Pareto)
- Decision support **euristico ma trasparente** (soglie documentate, regole auditabili, regimen evidence-based)
- Layer "twin" che è in realtà un **risk-mapping stateless da singolo assessment a parametri simulatore**

### Gap strutturale

Il gap tra lo stato attuale e un vero digital twin paziente-specifico è **strutturale, non incrementale**. Non si tratta di aggiungere feature a un sistema quasi-completo, ma di costruire un intero layer architetturale mancante:

1. **Stato latente persistente** — il twin deve essere un oggetto con stato che evolve nel tempo, non una funzione pura
2. **Modello di osservazione** — deve esistere un mapping formale $h(\mathbf{x}) \to \hat{\mathbf{y}}$ da stato latente a biomarker predetti
3. **Loop di aggiornamento** — il sistema deve reagire a nuovi dati paziente aggiornando lo stato twin
4. **Calibrazione** — i parametri del modello devono essere fittati sulla storia individuale del paziente, non derivati da breakpoint globali
5. **Confronto predetto-vs-osservato** — i residui devono essere calcolati, memorizzati, e usati per valutare la qualità del twin
6. **Branching controfattuale** — il valore clinico emerge dal "what-if" su terapie alternative

### Strategia raccomandata

1. **P0 (Fondamenta)**: Issues 1-3 — separare flusso paziente da scenario, persistere stato twin, implementare loop di aggiornamento base
2. **P1 (Core twin)**: Issues 4-7 — integrare terapia come input, observation model + residui, controfattuale, provenance
3. **P2 (Maturità)**: Issues 8-12 — dashboard, confidenza raccomandazioni, preparazione state estimation, reporting, documentazione

L'architettura proposta (`twin_engine/`) è progettata per essere **incrementalmente implementabile**: ogni modulo può essere aggiunto e testato indipendentemente, e il sistema esistente continua a funzionare invariato durante la transizione.

### Punto chiave

Il repository ha un'eccellente base per la transizione. Il layer dati paziente è completo, il simulatore ODE è solido, e l'infrastruttura Django (forms, views, API, template) è matura. Il lavoro richiesto è **concentrato nel layer intermedio** tra dati paziente e simulazione — esattamente dove risiede il digital twin.
