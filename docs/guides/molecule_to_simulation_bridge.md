# Molecule → Evidence → Simulation (bridge)

=== "IT"
    Questa guida spiega come leggere **insieme**:

    - **ChemTools**: struttura molecolare (3D / binding) + profili evidence-based (`mm_efficacy_profile`, `survival_impact`, `toxicity_profile`).
    - **Simulator**: KPI della simulazione (es. tumor reduction, healthy loss, AUC) per un regime/dosi specifiche.

    ## Cosa esiste oggi (e cosa no)

    - ✅ **Struttura + profilo MM** (ChemTools) sono disponibili nelle pagine Job (binding) e via endpoint asincrono.
    - ✅ **KPI simulati** (Simulator) sono disponibili per scenario/attempt.
    - ❌ **Non esiste ancora un join automatico** “questa molecola = questo attempt” (nessuna schermata unica che faccia match 1:1).

    In pratica: oggi il confronto si fa come **workflow riproducibile** (stesso nome farmaco → due viste diverse).

    ## Workflow consigliato (riproducibile)

    1) **Scegli un farmaco** (es. lenalidomide / bortezomib / daratumumab / carfilzomib).

    2) **ChemTools → struttura e profilo evidence-based**

       - Avvia un job ChemTools “binding” con un PDB ID + ligando (se lo hai) oppure usa SMILES/CID per proprietà/similarità.
       - Nella pagina job, le sezioni arricchite vengono caricate via AJAX e includono:
         - `mm_efficacy_profile`
         - `survival_impact`
         - `toxicity_profile`

       Approfondimenti:
       - Integrazione API e viewer: ../API_INTEGRATION_GUIDE.md
       - Stima efficacia: ../MM_EFFICACY_ESTIMATION.md
       - Stima survival/tossicità: ../SURVIVAL_TOXICITY_ESTIMATION.md

    3) **Simulator → KPIs su scenario e dosi**

       - Apri uno scenario nel Simulator e lancia una simulazione.
       - KPI principali:
         - tumor reduction
         - healthy loss
         - AUC per farmaco (se attivo)

       Approfondimento: Mathematical Models (Guides) e la gallery: simulations_gallery.md

    4) **Confronta (stesso farmaco, due “lenti”)**

       - ChemTools: “quanto è supportato da evidenza/letteratura e che profilo rischio-beneficio emerge”.
       - Simulator: “cosa succede nel modello dinamico (T/H) con quella dose/schedule/biologia”.

    ## Mappatura pratica (nomi comuni)

    Questa tabella è utile per mantenere coerenza tra linguaggio clinico, YAML PK/PD e ligandi tipici in PDB.

    | Farmaco (Simulator / YAML) | Classe | Esempi label/ligando (indicativi) |
    | --- | --- | --- |
    | lenalidomide | IMiD | `LEN`, `CC4` |
    | bortezomib | PI | `BTZ`, `BORTEZOMIB` |
    | daratumumab | anti-CD38 mAb | non tipicamente “ligando piccolo” (strutture sono complessi Ab/Ag) |
    | carfilzomib | PI (2ª gen) | `CFZ`, `CARFILZOMIB` |

    Nota: i codici ligando dipendono dall’entry PDB; ChemTools può anche inferire o arricchire via fonti esterne.

    ## Limiti e interpretazione

    - ChemTools produce stime *evidence-based* (non un trial simulator).
    - Simulator produce un modello dinamico *toy/educational* con guardrail numerici; non sostituisce decisioni cliniche.

=== "EN"
    This guide explains how to read **together**:

    - **ChemTools**: molecular structure (3D / binding) + evidence-based profiles (`mm_efficacy_profile`, `survival_impact`, `toxicity_profile`).
    - **Simulator**: simulation KPIs (e.g., tumor reduction, healthy loss, AUC) for a given regimen/doses.

    ## What exists today (and what doesn’t)

    - ✅ **Structure + MM evidence profile** (ChemTools) are available in binding Job pages and via async enriched sections.
    - ✅ **Simulated KPIs** (Simulator) are available per scenario/attempt.
    - ❌ There is **no automatic 1:1 join** yet (“this molecule = this attempt”).

    Today, comparison is done as a **reproducible workflow** (same drug name → two different views).

    ## Recommended workflow (reproducible)

    1) **Pick a drug** (lenalidomide / bortezomib / daratumumab / carfilzomib).

    2) **ChemTools → structure & evidence-based profile**

       - Start a ChemTools “binding” job using a PDB ID + ligand (if available), or use SMILES/CID for properties/similarity.
       - The enriched sections load via AJAX and include:
         - `mm_efficacy_profile`
         - `survival_impact`
         - `toxicity_profile`

       Deep dives:
       - API integration + viewer: ../API_INTEGRATION_GUIDE.md
       - MM efficacy estimation: ../MM_EFFICACY_ESTIMATION.md
       - Survival/toxicity estimation: ../SURVIVAL_TOXICITY_ESTIMATION.md

    3) **Simulator → KPIs on scenario & doses**

       - Open a scenario and run a simulation.
       - Core KPIs:
         - tumor reduction
         - healthy loss
         - AUC per active drug

       See the modeling guides and the gallery: simulations_gallery.md

    4) **Compare (same drug, two “lenses”)**

       - ChemTools: evidence/clinical literature support + risk/benefit summary.
       - Simulator: dynamic behavior (T/H) given dose/schedule/biology.

    ## Practical mapping (common names)

    | Drug (Simulator / YAML) | Class | Ligand labels (indicative) |
    | --- | --- | --- |
    | lenalidomide | IMiD | `LEN`, `CC4` |
    | bortezomib | PI | `BTZ`, `BORTEZOMIB` |
    | daratumumab | anti-CD38 mAb | not typically a small-molecule “ligand” |
    | carfilzomib | PI (2nd gen) | `CFZ`, `CARFILZOMIB` |

    Note: ligand codes depend on the specific PDB entry; ChemTools can infer/enrich via external sources.

    ## Limits & interpretation

    - ChemTools outputs *evidence-based* estimates (not a clinical trial simulator).
    - Simulator outputs an educational dynamic model with numerical guardrails; it’s not a clinical decision engine.
