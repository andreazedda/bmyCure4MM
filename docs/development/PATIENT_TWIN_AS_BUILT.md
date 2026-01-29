# Patient Twin — Implementazione attuale (as-built)

Data: 2026-01-04  
**Aggiornato:** 2026-01-05 (aggiunta guida pratica UI)

Questo documento descrive **esattamente** come è gestito oggi il "Patient Twin" nel repository, **cos'è**, **come viene calcolato**, **come entra nel simulatore**, **come viene serializzato**, e **se/come è temporale**.

> **Per iniziare subito:** vedi il [Tutorial Pratico (IT)](../it/TUTORIAL_PRIMO_USO.md) o [Tutorial (EN)](../en/TUTORIAL_FIRST_RUN.md)  
> Nota tecnica breve: [docs/en/patient_twin.md](../en/patient_twin.md)

---

## 0) Come si usa (workflow UI pratico)

**TL;DR:** Il Patient Twin NON è un paziente virtuale che crei da zero — è un **motore di calcolo** che trasforma valori di laboratorio *reali* in parametri biologici personalizzati per la simulazione.

### Flusso Pratico (5 passi):

```
1️⃣ Crea Paziente (Clinic → Patients → New Patient)
          ↓
2️⃣ Aggiungi Assessment con Lab (R-ISS, LDH, β2M, FLC)
          ↓        ← Questo è lo "snapshot" per il Twin
3️⃣ Vai al Simulatore (Simulator → Browse Scenarios)
          ↓
4️⃣ Abilita Patient Twin (sezione "Advanced & Twin")
          ↓        ← Seleziona l'Assessment + modalità "Auto"
5️⃣ Lancia & Vedi Risultati (curve tumorali, twin_params.json)
```

**Dove trovi il Tutorial:**
- **Dashboard:** card blu in alto con "🚀 New to the Platform? Start Here!"
- **Menu:** link "🎓 Tutorial" in alto a destra
- **Simulatore:** bottone "📖 How to use Twin?" nella sezione Advanced & Twin

---

## 1) Cos’è il Patient Twin (oggi)

**Definizione operativa:** il Patient Twin è un **trasformatore deterministico** che prende **un singolo** record di laboratorio/valutazione clinica (`clinic.models.Assessment`) e produce un dizionario di **parametri scalari** (float) che possono essere usati dal simulatore PK/PD come “override” dei parametri biologici base.

In pratica:

- Input: un `Assessment` (snapshot di laboratorio/clinica) con campi come R-ISS, LDH, β2M, FLC ratio.
- Output: un payload `dict[str, float]` con:
  - `risk_score` in [0, 1]
  - `tumor_growth_rate`
  - `healthy_growth_rate`
  - `carrying_capacity_tumor`
  - `carrying_capacity_healthy`
  - `immune_compromise_index`

File principali:

- Implementazione core: [simulator/twin.py](https://github.com/andreazedda/bmyCure4MM/blob/main/simulator/twin.py)
- Config dei pesi e range: [simulator/presets/twin_risk.yaml](https://github.com/andreazedda/bmyCure4MM/blob/main/simulator/presets/twin_risk.yaml)
- Applicazione al modello e serializzazione artifact: [simulator/models.py](https://github.com/andreazedda/bmyCure4MM/blob/main/simulator/models.py)

---

## 2) È “temporale”? (risposta precisa)

### 2.1 Stato attuale: **non** è un modello temporale interno

Il Patient Twin **non** è una serie temporale e **non** è un modello dinamico che evolve autonomamente.

- Usa **una sola** istanza `Assessment` (che ha un campo `date`, ma la data **non** entra nel calcolo).
- Produce parametri costanti che vengono usati per **tutta** la durata della simulazione (orizzonte in giorni).

### 2.2 Effetto sul tempo: indiretto, tramite parametri costanti

L’effetto “nel tempo” esiste solo perché quei parametri (es. growth rate) entrano nelle equazioni differenziali/iterative del `MathematicalModel`:

- Se il Twin aumenta `tumor_growth_rate`, la curva tumorale cresce più rapidamente lungo l’asse tempo.
- Se il Twin riduce `carrying_capacity_healthy`, la componente healthy può saturare prima.

Ma i parametri del Twin **non cambiano durante la simulazione** (non ci sono update per step).

---

## 3) Modello dati: quali campi usa

### 3.1 Modello `Assessment`

Il modello è in [clinic/models.py](https://github.com/andreazedda/bmyCure4MM/blob/main/clinic/models.py) (`class Assessment`). Campi rilevanti:

- `r_iss` (string, choices I/II/III)
- `ldH_u_l` (DecimalField)  
  Nota: il nome del campo è letteralmente `ldH_u_l` (H maiuscola) e viene letto così anche nel Twin.
- `beta2m_mg_l` (DecimalField)
- `flc_ratio` (DecimalField)

### 3.2 Casting e valori null

In [simulator/twin.py](https://github.com/andreazedda/bmyCure4MM/blob/main/simulator/twin.py) la funzione `_maybe_float(value)` gestisce:

- `None` o stringa vuota → `None`
- altrimenti prova `float(value)`
- in caso di `TypeError/ValueError` → `None`

Quindi i `DecimalField` di Django vengono convertiti a float quando possibile.

---

## 4) Configurazione (YAML) e parametri

### 4.1 File di config

Il Twin legge **sempre** questo file:

- [simulator/presets/twin_risk.yaml](https://github.com/andreazedda/bmyCure4MM/blob/main/simulator/presets/twin_risk.yaml)

Caricamento:

- `simulator/twin.py::_load_twin_config()` usa `Path(__file__).parent / "presets" / "twin_risk.yaml"`
- Se il file non esiste → `FileNotFoundError` (hard-fail)
- È cache-ato con `@lru_cache` (quindi letto una sola volta per processo Python)

### 4.2 Contenuto YAML (attuale)

Pesi (`weights`):

- `riss: 0.35`
- `ldh: 0.2`
- `beta2m: 0.25`
- `flc_ratio: 0.2`

Map R-ISS (`r_iss_map`):

- `I → 0.1`
- `II → 0.6`
- `III → 1.0`

LDH breakpoints (`ldh_breakpoints`):

- `normal_upper: 250`
- `high: 550`

Beta2M range (`beta2m_range`):

- `min: 2.0`
- `max: 15.0`

FLC ratio range (`flc_ratio_range`):

- normalità: `0.26–1.65`
- estremi: `min: 0.05`, `max: 20.0`

Mapping a parametri simulatore:

- `growth_mapping.tumor_growth_rate: min 0.017, max 0.04`
- `growth_mapping.healthy_growth_rate: min 0.012, max 0.02`
- `carrying_capacity.tumor: min 5.0e8, max 2.5e9`
- `carrying_capacity.healthy: min 2.5e11, max 6.0e11`
- `immune_compromise_index: min 0.85, max 1.3`

---

## 5) Algoritmo: come si calcola `risk_score`

La funzione entry-point è:

- `simulator.twin.build_patient_twin(assessment: Assessment) -> Dict[str, float]`

### 5.1 Componenti normalizzate (0..1)

#### (A) R-ISS

In `simulator/twin.py::_score_riss()`:

- legge `assessment.r_iss`, lo upper-case-a
- mapping di default: `{I:0.2, II:0.6, III:1.0}` ma **nel repo** viene override-ato dal YAML (I=0.1)
- se valore mancante o non riconosciuto → usa il valore di `II` (default 0.6)

Output: float tipicamente in {0.1, 0.6, 1.0}

#### (B) LDH

In `simulator/twin.py::_score_ldh()`:

- legge `assessment.ldH_u_l`
- se `None` → ritorna `0.5` (valore neutro)
- se `LDH <= normal_upper` → `0.1`
- se `LDH >= high` → `1.0`
- altrimenti interpolazione lineare tra `0.1` e `1.0`:

$$
score = 0.1 + 0.9 \cdot \frac{LDH - normal}{high - normal}
$$

Con `normal=250`, `high=550`.

#### (C) β2-microglobulina

In `simulator/twin.py::_score_beta2m()`:

- legge `assessment.beta2m_mg_l`
- se `None` → `0.5`
- se `<= min` → `0.1`
- se `>= max` → `1.0`
- altrimenti interpolazione lineare **0..1** tra `min` e `max`:

$$
score = \frac{\beta2m - min}{max - min}
$$

Nota: qui il range scala fino a 1.0 ma non forza un offset 0.1 come LDH; quindi per valori intermedi può restituire 0.2, 0.3, ecc.

#### (D) FLC ratio

In `simulator/twin.py::_score_flc_ratio()`:

- legge `assessment.flc_ratio`
- se `None` → `0.4`
- se nel range di normalità `0.26–1.65` → `0.1`
- se sopra 1.65 → cresce linearmente fino a 1.0 a `max`:

$$
score = \min(1.0,\ 0.1 + 0.9\cdot \frac{FLC - high}{max - high})
$$

- se sotto 0.26 → simmetrico verso il basso (considera “abnorme” anche ratio molto basso):

$$
score = \min(1.0,\ 0.1 + 0.9\cdot \frac{low - FLC}{low - min})
$$

### 5.2 Combinazione pesata

In `build_patient_twin()`:

- calcola `total_weight = sum(weights)`, fallback 1.0
- calcola:

$$
weighted\_sum = w_{riss} \cdot riss + w_{ldh} \cdot ldh + w_{beta2m}\cdot beta2m + w_{flc}\cdot flc
$$

- normalizza:

$$
risk\_score = clamp(0,1,\ \frac{weighted\_sum}{total\_weight})
$$

Con i pesi attuali, `total_weight = 1.0`.

---

## 6) Mapping: da `risk_score` a parametri del simulatore

Una volta ottenuto `risk_score`, `build_patient_twin()` calcola:

- `tumor_growth_rate = lerp(growth_mapping.tumor_growth_rate, risk_score)`
- `healthy_growth_rate = lerp(growth_mapping.healthy_growth_rate, 1 - risk_score)`
- `carrying_capacity_tumor = lerp(carrying_capacity.tumor, risk_score)`
- `carrying_capacity_healthy = lerp(carrying_capacity.healthy, 1 - risk_score)`
- `immune_compromise_index = lerp(immune_compromise_index, risk_score, default=1.0)`

Dove `lerp` è in `simulator/twin.py::_lerp()`:

- legge `min` e `max` dal config
- clamp di `weight` in [0,1]
- restituisce:

$$
value = min + (max - min)\cdot weight
$$

**Interpretazione importante:**

- rischio alto (risk_score→1) →
  - `tumor_growth_rate` verso il max
  - `healthy_growth_rate` verso il min (perché usa 1-risk)
  - `carrying_capacity_tumor` verso il max
  - `carrying_capacity_healthy` verso il min (perché usa 1-risk)
  - `immune_compromise_index` verso il max

---

## 7) Integrazione nel simulatore: quando viene applicato

### 7.1 Punto di integrazione

L’integrazione avviene in:

- `simulator.models.SimulationAttempt.run_model()` ([simulator/models.py](https://github.com/andreazedda/bmyCure4MM/blob/main/simulator/models.py))

Flusso (semplificato ma fedele):

1. `params = dict(self.parameters or {})`
2. `resolved_params = dict(params)`
3. `use_predlab = settings.PREDLAB_V2` (default `False` se env non la abilita)
4. Se `resolved_params.get("use_twin", True)`:
   - prende `assessment_id = resolved_params.get("twin_assessment_id") or resolved_params.get("assessment_id")`
   - se esiste, prova `Assessment.objects.get(pk=assessment_id)`

   **Nota permessi (stato attuale PR1.1):** l’Assessment viene caricato solo se:

   - `Patient.owner == attempt.user`, oppure
   - l’utente è `is_staff` / `Simulator Editors`.

   Se `Patient.owner` è `NULL`, l’Assessment è accessibile solo a staff/editor (fail-closed).

   - se trovato:
     - `twin_payload = build_patient_twin(twin_assessment)`
     - se `twin_biology_mode == "auto"`: `resolved_params = _merge_twin_parameters(resolved_params, twin_payload)`

### 7.2 Gating: `PREDLAB_V2`

Con PR1, il Patient Twin **non è più gated** da `PREDLAB_V2`:

- `PREDLAB_V2 = os.environ.get("PREDLAB_V2", "0") == "1"` in [mmportal/settings.py](https://github.com/andreazedda/bmyCure4MM/blob/main/mmportal/settings.py)
- Anche con `PREDLAB_V2=0`, il Twin può essere calcolato e serializzato come artifact, e (in modalità `auto`) può applicare override biologici.

**Cosa resta gated da `PREDLAB_V2`:** la risoluzione PK/PD “avanzata” tramite `pharmaco_registry.resolve(...)` (se `PREDLAB_V2` è off, vengono usati i parametri PK/PD di default hard-coded).

### 7.3 Come si passa l’Assessment (stato attuale reale)

Per applicare il Twin serve un `assessment_id` nei parametri.

Con PR1, `SimulationParameterForm` espone:

- `twin_assessment_id` (selector)
- `twin_biology_mode` (`auto|manual`)

La view `simulate_scenario` salva `parameters=form.cleaned_data` (quindi i campi UI arrivano in `SimulationAttempt.parameters`).

Nota permessi (stato attuale PR1): per evitare enumerazione via ID, il selector e la preview sono abilitati solo per utenti `is_staff` o nel gruppo `Simulator Editors`.

Con PR1.1 la policy è:

- utenti non-staff: vedono/possono usare solo gli Assessment dei propri pazienti (`Patient.owner == request.user`)
- staff/editor: vedono anche record con `owner=NULL`

Preview endpoint (PR1):

- `/sim/api/twin/preview/?id=<assessment_id>`
- Alias: `/simulator/api/twin/preview/?id=<assessment_id>`

---

## 8) Merge policy: quando il Twin sovrascrive i parametri

La logica di merge è in:

- `SimulationAttempt._merge_twin_parameters(params, twin_payload)` ([simulator/models.py](https://github.com/andreazedda/bmyCure4MM/blob/main/simulator/models.py))

Regole:

- Solo queste chiavi sono “ammissibili”:
  - `tumor_growth_rate`
  - `healthy_growth_rate`
  - `carrying_capacity_tumor`
  - `carrying_capacity_healthy`
  - `immune_compromise_index`

- Per ogni chiave, il Twin **imposta** `merged[key] = twin_value` solo se:
  - la chiave **non** è presente nei params, **oppure**
  - il valore corrente è uno dei sentinel `{"auto","AUTO","Auto",""}`, **oppure**
  - il valore corrente è `None`

Questo è coerente con l’idea “override solo se marcato auto”.

**Dettaglio importante (implicazione):**

- Nel form UI, campi come `tumor_growth_rate`/`healthy_growth_rate` sono `FloatField` e vengono tipicamente sempre inviati come numeri → quindi **non** sono `auto` e **non** vengono sovrascritti.
- Il test [simulator/tests/test_twin.py](https://github.com/andreazedda/bmyCure4MM/blob/main/simulator/tests/test_twin.py) passa un `parameters` senza `tumor_growth_rate` ecc., quindi lì il Twin viene applicato.

---

## 9) Serializzazione e artifact `twin_params.json`

Se `twin_payload` è non-vuoto, `run_model()` scrive un file JSON su disco:

- directory: `MEDIA_ROOT/sim_data/`
- filename: `attempt_{attempt_id}_twin.json`

Contenuto scritto:

- tutti i campi di `twin_payload` (incluso `risk_score`)
- più `assessment_id` (id dell’Assessment usato)

Link nel risultato:

- `attempt.results["twin_params.json"] = MEDIA_URL + "/sim_data/" + filename`

Inoltre:

- `attempt.results` contiene anche il CSV e il plot HTML
- `attempt.artifacts` copia i link (eccetto `generated_at`) e aggiunge `seed` se presente

---

## 10) Test esistenti e cosa verificano

File test:

- [simulator/tests/test_twin.py](https://github.com/andreazedda/bmyCure4MM/blob/main/simulator/tests/test_twin.py)

Casi principali:

1. `test_risk_score_monotonic_with_stage`
   - costruisce due Assessment (low vs high)
   - verifica monotonia: tumor growth aumenta, healthy carrying diminuisce

2. `test_twin_serialization_written_to_results`
   - `@override_settings(PREDLAB_V2=True)`
   - crea `SimulationAttempt` con `twin_assessment_id` e `use_twin=True`
   - esegue `run_model()`
   - verifica che `attempt.results` contenga `twin_params.json` e che il file contenga `risk_score`

---

## 11) Limitazioni/assunzioni attuali (importanti per “come è gestito adesso”)

- **Gating:** con PR1 il Twin non è più gated da `PREDLAB_V2`; resta gated solo la risoluzione PK/PD avanzata (vedi §7.2).
- **Selezione da UI:** con PR1 esiste `twin_assessment_id` nel form; quindi il flusso UI→Attempt→Twin è percorribile.
- **Access control minimale (PR1.1):** un utente non-staff può usare Twin solo sugli Assessment dei pazienti di cui è owner (`Patient.owner`). Staff/editor possono vedere tutto.
- **Override condizionale:** l’override biologico avviene solo se `twin_biology_mode == "auto"`; in modalità manual i valori restano quelli immessi.
- **Non temporale:** usa un singolo snapshot; non usa `Assessment.date` né la storia.
- **Hard dependency dal YAML:** mancando `twin_risk.yaml`, il Twin lancia eccezione.

---

## 12) Checklist rapida per riprodurre l’uso reale (oggi)

Per vedere il Twin applicato **con l’implementazione attuale**, devi:

1. Avere un `Assessment` nel DB.
2. Creare una `SimulationAttempt` con `parameters` che includano:
   - `twin_assessment_id` (o `assessment_id`)
   - `use_twin: True` (opzionale; default True se chiave assente)
   - **non** includere (o impostare a `"auto"`) i parametri che vuoi far derivare:
     - `tumor_growth_rate`, `healthy_growth_rate`, `carrying_capacity_*`, `immune_compromise_index`
3. Avviare con `PREDLAB_V2=1`.

A quel punto troverai `twin_params.json` tra i risultati e l’override avrà effetto sul modello.
