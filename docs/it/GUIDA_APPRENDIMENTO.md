# Guida Rapida per Imparare a Usare la Piattaforma

**Per chi vuole imparare come simulare trattamenti per Mieloma Multiplo**

---

## 🎯 Obiettivo

Capire **come funziona la simulazione di trattamenti** e **come i pazienti rispondono ai farmaci** usando la piattaforma bmyCure4MM.

---

## ✨ Setup Rapidissimo (2 minuti)

### 1. Installa e carica dati demo

```bash
git clone https://github.com/andreazedda/bmyCure4MM.git
cd bmyCure4MM
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
./quick_start.sh
python manage.py runserver
```

### 2. Apri il browser

Vai su: **http://127.0.0.1:8000**

---

## 📚 Cosa Trovi Già Pronto

### ✅ 3 Pazienti Demo (per vedere diversi scenari di rischio)

1. **Mario Rossi** (MRN: DEMO001) - **Alto Rischio (R-ISS III)**
   - LDH: 600 U/L (molto alto)
   - β2-microglobulina: 12.0 mg/L (molto alto)
   - FLC ratio: 8.5 (molto sbilanciato)
   - → Tumore aggressivo, necessita terapia intensiva

2. **Anna Bianchi** (MRN: DEMO002) - **Rischio Intermedio (R-ISS II)**
   - LDH: 320 U/L (moderatamente elevato)
   - β2-microglobulina: 5.5 mg/L (elevato)
   - FLC ratio: 3.2 (sbilanciato)
   - → Scenario tipico, buon candidato per triplet standard

3. **Giuseppe Verdi** (MRN: DEMO003) - **Basso Rischio (R-ISS I)**
   - LDH: 200 U/L (normale)
   - β2-microglobulina: 2.8 mg/L (quasi normale)
   - FLC ratio: 1.1 (normale)
   - → Prognosi favorevole, risponde bene a terapie standard

---

## 🧪 Il Tuo Primo Test: Confronta Come Rispondono i Pazienti

### Scenario: Voglio vedere come **LO STESSO REGIME** (es. VRd: Bortezomib + Lenalidomide + Dexamethasone) funziona su pazienti con rischio diverso

### 📋 Procedura:

#### PASSO 1: Vai al Simulatore
- Dashboard → Card "🧪 Treatment Simulator" **O**
- Menu in alto → "🧪 Simulator"

#### PASSO 2: Scegli uno Scenario
- Clicca su uno scenario esistente (es. "VRd Triplet")
- Questo apre il form di simulazione

#### PASSO 3: Configura con Patient Twin (primo test: Mario Rossi - alto rischio)

1. **Sezione "1. Baseline"**
   - Preset: Scegli "standard_triplet" (popola dosi automaticamente)

2. **Scorri giù a "4. Advanced & Twin"**
   - ✅ Abilita "Use Patient Twin"
   - Seleziona Assessment: **"Mario Rossi - 2026-01-01"**
   - Twin Biology Mode: **"Auto"** (lascia che il Twin calcoli parametri biologici)

3. **Time Horizon:** 180 giorni (6 mesi)

4. **Cohort Size:** 1 (per ora)

#### PASSO 4: Lancia la Simulazione
- Clicca "Run Simulation" in basso
- Aspetta 10-20 secondi

#### PASSO 5: Analizza i Risultati (Mario Rossi - alto rischio)

Vedrai:
- **Grafico:** Curve rosse (tumore) e verdi (cellule sane) nel tempo
- **KPI:**
  - Tumor Reduction: probabilmente **~85-90%** (discreto, ma il tumore è aggressivo)
  - Healthy Loss: probabilmente **~20-25%** (tossicità moderata)
  - Time to Recurrence: potrebbe essere **breve** (tumore cresce veloce)

- **File `twin_params.json`:**
  - `risk_score`: ~0.80-0.85 (alto rischio)
  - `tumor_growth_rate`: ~0.035-0.038 (vicino al max 0.04 → cresce veloce!)
  - `healthy_growth_rate`: ~0.013 (vicino al min → recupero lento)

**Interpretazione:** Tumore risponde ai farmaci, ma con questo profilo di rischio alto serve monitoraggio stretto e potenzialmente intensificazione terapia.

---

#### PASSO 6: Ripeti con Giuseppe Verdi (basso rischio)

Stessa procedura, ma:
- Twin Assessment: **"Giuseppe Verdi - 2026-01-03"**

**Risultati attesi (Giuseppe - basso rischio):**
- Tumor Reduction: probabilmente **~95-98%** (ottima risposta!)
- Healthy Loss: probabilmente **~10-15%** (tossicità bassa)
- Time to Recurrence: probabilmente **molto lungo** (tumore cresce lento)

- `twin_params.json`:
  - `risk_score`: ~0.15-0.25 (basso rischio)
  - `tumor_growth_rate`: ~0.018-0.020 (vicino al min → cresce lento)
  - `healthy_growth_rate`: ~0.018-0.019 (vicino al max → recupero rapido)

**Interpretazione:** Ottima risposta, bassa tossicità. Paziente ideale per questo regime.

---

## 🔍 Cosa Hai Imparato

### 1. **Lo stesso regime ha effetti diversi su pazienti diversi**
   - Mario (alto rischio): risposta discreta, tossicità moderata
   - Giuseppe (basso rischio): risposta ottima, tossicità bassa

### 2. **Il Patient Twin personalizza la simulazione**
   - Prende valori lab reali (R-ISS, LDH, β2M, FLC)
   - Calcola parametri biologici (growth rate, carrying capacity)
   - Modifica la dinamica tumorale nella simulazione

### 3. **I parametri del Twin spiegano il perché**
   - `tumor_growth_rate` alto → tumore cresce veloce → serve terapia aggressiva
   - `healthy_growth_rate` basso → recupero lento → attenzione a tossicità cumulativa

---

## 🚀 Esperimenti Successivi

### Esperimento 2: Cambia le Dosi
- Riduci Lenalidomide da 25 a 15 mg (per Mario)
- Confronta: efficacia scende? Tossicità migliora?
- **Obiettivo:** Capire il trade-off efficacia/tossicità

### Esperimento 3: Cambia il Regime
- Prova uno scenario diverso (es. Rd doublet - senza Bortezomib)
- Confronta KPI con VRd triplet
- **Obiettivo:** Capire quale regime è meglio per quale paziente

### Esperimento 4: Cohort Size > 1
- Imposta Cohort Size = 50
- Vedrai **bande di variabilità** nel grafico
- **Obiettivo:** Capire l'incertezza stocastica del modello

---

## ❓ FAQ per Principianti

### Q: Ma sto simulando su pazienti virtuali o dati reali?
**A:** I pazienti demo (Mario, Anna, Giuseppe) sono **dati sintetici realistici**. I valori lab sono tipici di pazienti reali con alto/medio/basso rischio. La simulazione è un modello matematico (equazioni differenziali) che predice dinamica tumorale.

### Q: Perché non vedo la curve cambiare se cambio solo il Twin?
**A:** Il Twin cambia **solo i parametri biologici** (growth rate, carrying capacity). Le curve cambiano se quei parametri influenzano significativamente la dinamica. In alcuni regimi, l'effetto farmacologico "domina" e la differenza è sottile — guarda i KPI numerici!

### Q: Cosa significa "Time Horizon"?
**A:** La **durata della simulazione** in giorni. 180 giorni = 6 mesi di terapia simulata.

### Q: Cosa significa "Cohort Size"?
**A:** Numero di **"repliche" della simulazione** per stimare variabilità stocastica. Cohort Size = 50 → il modello gira 50 volte con piccole variazioni casuali e ti mostra la media + banda di confidenza.

### Q: Posso creare i miei pazienti?
**A:** **Sì!** Vai su Clinic → Patients → New Patient, poi aggiungi Assessment con valori lab. Poi usa quello nel simulatore come Twin.

---

## 📖 Link Utili

- **Tutorial Completo:** Dashboard → card blu "🚀 New to Platform?"
- **Documentazione Tecnica Twin:** `docs/development/PATIENT_TWIN_AS_BUILT.md`
- **Help Contestuale:** Clicca i "?" accanto ai campi nel form simulatore

---

**🎉 Buon apprendimento! Se qualcosa non è chiaro, apri un issue su GitHub o chiedi in chat.**
