"""Clinical regimen presets for simulator workflows."""

from __future__ import annotations

PRESETS: dict[str, dict[str, object]] = {
    "VRd": {
        "label": "VRd (Bortezomib + Lenalidomide + Dexamethasone)",
        "description_en": "Triplet induction with proteasome inhibitor + IMiD. Standard for transplant-eligible patients.",
        "description_it": "Tripletta di induzione con inibitore del proteasoma + IMiD. Standard per pazienti trapiantabili.",
        
        # 🎮 Educational Story
        "story_en": {
            "title": "Virtual Patient: Alex - The Young Warrior 🎯",
            "background": "Meet Alex, a 58-year-old teacher newly diagnosed with multiple myeloma. Alex is fit, active, and a good candidate for stem cell transplant later. We need to reduce the tumor burden quickly and safely before the transplant procedure.",
            "challenge": "Your goal: Achieve >80% tumor reduction while keeping healthy cells above 60%. This 'triplet therapy' combines three weapons: a proteasome inhibitor (the shredder ✂️), an immune modulator (the booster 🚀), and steroids (the firefighters 🧯).",
            "why_these_drugs": "Bortezomib blocks cancer's waste disposal (proteasomes), Lenalidomide wakes up the immune system, and Dexamethasone reduces inflammation. Together, they're like a coordinated SWAT team!",
            "expected_outcome": "In 6 months (168 days), we expect 70-90% tumor reduction. Watch the neuropathy meter - Bortezomib can damage nerves, so we give it only twice weekly.",
            "learning_points": [
                "This is the 'standard of care' - proven effective in thousands of patients",
                "Intermittent dosing (days 1,4,8,11) prevents toxicity accumulation",
                "Lenalidomide given daily for 3 weeks, then 1 week rest (21/7 schedule)",
                "High interaction strength (0.08) means drugs synergize well",
            ],
        },
        "story_it": {
            "title": "Paziente Virtuale: Alex - Il Giovane Guerriero 🎯",
            "background": "Conosci Alex, un insegnante di 58 anni con diagnosi recente di mieloma multiplo. Alex è in forma, attivo e un buon candidato per trapianto di cellule staminali. Dobbiamo ridurre il carico tumorale rapidamente e in sicurezza prima del trapianto.",
            "challenge": "Tuo obiettivo: Raggiungere >80% di riduzione tumorale mantenendo le cellule sane sopra il 60%. Questa 'tripletta terapeutica' combina tre armi: un inibitore del proteasoma (il tritatore ✂️), un modulatore immunitario (il potenziatore 🚀) e steroidi (i pompieri 🧯).",
            "why_these_drugs": "Bortezomib blocca lo smaltimento rifiuti del cancro (proteasomi), Lenalidomide sveglia il sistema immunitario, e Desametasone riduce l'infiammazione. Insieme, sono come una squadra SWAT coordinata!",
            "expected_outcome": "In 6 mesi (168 giorni), ci aspettiamo 70-90% di riduzione tumorale. Guarda il metro della neuropatia - Bortezomib può danneggiare i nervi, quindi lo diamo solo due volte a settimana.",
            "learning_points": [
                "Questo è lo 'standard of care' - provato efficace in migliaia di pazienti",
                "Dosaggio intermittente (giorni 1,4,8,11) previene accumulo di tossicità",
                "Lenalidomide dato giornalmente per 3 settimane, poi 1 settimana di riposo (schema 21/7)",
                "Alta forza di interazione (0.08) significa che i farmaci si sinergizzano bene",
            ],
        },
        
        "bounds_pct": 20,
        "default_params": {
            "baseline_tumor_cells": 1.2e9,
            "baseline_healthy_cells": 4.5e11,
            "lenalidomide_dose": 25.0,
            "bortezomib_dose": 1.3,
            "daratumumab_dose": 0.0,
            "time_horizon": 168.0,
            "tumor_growth_rate": 0.023,
            "healthy_growth_rate": 0.015,
            "interaction_strength": 0.08,
        },
        "schedule": {
            "lenalidomide_dose": {"type": "days_on_off", "on": 21, "off": 7},
            "bortezomib_dose": {"type": "days_list", "days": [1, 4, 8, 11]},
            "daratumumab_dose": {"type": "days_list", "days": []},
        },
    },
    "Dara-Rd": {
        "label": "Daratumumab + Lenalidomide",
        "description_en": "Antibody + IMiD backbone for high-depth responses with lower neuropathy risk.",
        "description_it": "Combinazione anticorpo + IMiD per risposte profonde con minore rischio di neuropatia.",
        
        # 🎮 Educational Story
        "story_en": {
            "title": "Virtual Patient: Maria - The Precision Strike 🎯",
            "background": "Maria is a 65-year-old librarian with myeloma. She has mild neuropathy from diabetes, so we want to avoid Bortezomib (which can worsen nerve damage). Good news: Daratumumab is like a 'smart missile' that targets cancer cells precisely!",
            "challenge": "Your mission: Achieve deep remission (>85% tumor reduction) without triggering neuropathy. We're using a 'gentler but precise' approach - no nerve-toxic drugs, just targeted antibody + immune modulation.",
            "why_these_drugs": "Daratumumab is an antibody that finds CD38 proteins on myeloma cells and marks them for destruction. Lenalidomide boosts the immune system to finish the job. Think of it as 'guided missiles + army reinforcements'.",
            "expected_outcome": "In 7 months (196 days), expect 75-95% tumor reduction with minimal neuropathy. Daratumumab infusions are weekly at first (to saturate targets), then monthly for maintenance.",
            "learning_points": [
                "Daratumumab targets CD38 - found on 95% of myeloma cells",
                "Weekly infusions for month 1-2 (loading phase), then less frequent",
                "Lower neuropathy risk makes this ideal for patients with nerve issues",
                "Longer treatment duration (196 vs 168 days) for sustained responses",
            ],
        },
        "story_it": {
            "title": "Paziente Virtuale: Maria - Il Colpo di Precisione 🎯",
            "background": "Maria è una bibliotecaria di 65 anni con mieloma. Ha una lieve neuropatia dal diabete, quindi vogliamo evitare Bortezomib (che può peggiorare il danno nervoso). Buone notizie: Daratumumab è come un 'missile intelligente' che colpisce le cellule cancerose con precisione!",
            "challenge": "La tua missione: Raggiungere remissione profonda (>85% riduzione tumorale) senza scatenare neuropatia. Stiamo usando un approccio 'più gentile ma preciso' - nessun farmaco neurotossico, solo anticorpo mirato + modulazione immunitaria.",
            "why_these_drugs": "Daratumumab è un anticorpo che trova proteine CD38 sulle cellule mielomatose e le marca per la distruzione. Lenalidomide potenzia il sistema immunitario per finire il lavoro. Pensalo come 'missili guidati + rinforzi dell'esercito'.",
            "expected_outcome": "In 7 mesi (196 giorni), aspettati 75-95% di riduzione tumorale con neuropatia minima. Le infusioni di Daratumumab sono settimanali all'inizio (per saturare i bersagli), poi mensili per mantenimento.",
            "learning_points": [
                "Daratumumab colpisce CD38 - presente sul 95% delle cellule mielomatose",
                "Infusioni settimanali per mese 1-2 (fase di carico), poi meno frequenti",
                "Rischio di neuropatia più basso lo rende ideale per pazienti con problemi nervosi",
                "Durata trattamento più lunga (196 vs 168 giorni) per risposte sostenute",
            ],
        },
        
        "bounds_pct": 20,
        "default_params": {
            "baseline_tumor_cells": 9.5e8,
            "baseline_healthy_cells": 5.0e11,
            "lenalidomide_dose": 25.0,
            "bortezomib_dose": 0.0,
            "daratumumab_dose": 16.0,
            "time_horizon": 196.0,
            "tumor_growth_rate": 0.02,
            "healthy_growth_rate": 0.016,
            "interaction_strength": 0.12,
        },
        "schedule": {
            "lenalidomide_dose": {"type": "days_on_off", "on": 21, "off": 7},
            "bortezomib_dose": {"type": "days_list", "days": []},
            "daratumumab_dose": {"type": "qweekly", "weeks": [1, 2, 3, 4], "day": 1},
        },
    },
    "KRd": {
        "label": "Carfilzomib + Lenalidomide",
        "description_en": "Second-generation PI + IMiD for aggressive or high-risk disease.",
        "description_it": "Inibitore del proteasoma di seconda generazione + IMiD per malattia aggressiva o ad alto rischio.",
        
        # 🎮 Educational Story
        "story_en": {
            "title": "Virtual Patient: Carlos - The Hard Mode Challenge 💪",
            "background": "Carlos is a 52-year-old engineer with aggressive, high-risk myeloma (unfavorable genetics). Standard therapy might not be enough - we need the 'heavy artillery'. Carfilzomib is Bortezomib's stronger cousin, with more firepower but requires careful monitoring.",
            "challenge": "HARD MODE activated! 🔥 Achieve >90% tumor reduction in aggressive disease while managing toxicity. Carfilzomib hits harder but has cardiac and renal risks. This is like playing on 'Expert difficulty' - precision timing is everything!",
            "why_these_drugs": "Carfilzomib is a 'second-generation proteasome inhibitor' - more selective, more potent, and given twice weekly for maximum effect. Combined with Lenalidomide, it's our best weapon against high-risk disease. Think 'precision bombing campaign'.",
            "expected_outcome": "In 6 months (168 days), aim for 85-95% tumor reduction. Twice-weekly dosing (days 1,2,8,9,15,16) = intensive schedule. Monitor blood pressure and kidney function closely!",
            "learning_points": [
                "Carfilzomib binds irreversibly to proteasomes (permanent blockade)",
                "Twice-weekly dosing for first 3 weeks of each 28-day cycle",
                "Lower neuropathy than Bortezomib, but higher cardiovascular risk",
                "Reserved for high-risk or relapsed disease - it's your 'ace card'",
            ],
        },
        "story_it": {
            "title": "Paziente Virtuale: Carlos - La Sfida Modalità Difficile 💪",
            "background": "Carlos è un ingegnere di 52 anni con mieloma aggressivo e ad alto rischio (genetica sfavorevole). La terapia standard potrebbe non bastare - serve l''artiglieria pesante'. Carfilzomib è il cugino più forte di Bortezomib, con più potenza di fuoco ma richiede monitoraggio attento.",
            "challenge": "MODALITÀ DIFFICILE attivata! 🔥 Raggiungi >90% di riduzione tumorale in malattia aggressiva gestendo la tossicità. Carfilzomib colpisce più forte ma ha rischi cardiaci e renali. È come giocare in 'Difficoltà Esperto' - il tempismo preciso è tutto!",
            "why_these_drugs": "Carfilzomib è un 'inibitore del proteasoma di seconda generazione' - più selettivo, più potente, e dato due volte a settimana per massimo effetto. Combinato con Lenalidomide, è la nostra arma migliore contro malattia ad alto rischio. Pensa 'campagna di bombardamento di precisione'.",
            "expected_outcome": "In 6 mesi (168 giorni), punta a 85-95% di riduzione tumorale. Dosaggio bisettimanale (giorni 1,2,8,9,15,16) = programma intensivo. Monitora pressione sanguigna e funzione renale attentamente!",
            "learning_points": [
                "Carfilzomib si lega irreversibilmente ai proteasomi (blocco permanente)",
                "Dosaggio bisettimanale per prime 3 settimane di ogni ciclo di 28 giorni",
                "Neuropatia più bassa di Bortezomib, ma rischio cardiovascolare più alto",
                "Riservato per malattia ad alto rischio o recidiva - è la tua 'carta asso'",
            ],
        },
        
        "bounds_pct": 20,
        "default_params": {
            "baseline_tumor_cells": 1.0e9,
            "baseline_healthy_cells": 4.8e11,
            "lenalidomide_dose": 25.0,
            "bortezomib_dose": 1.1,
            "daratumumab_dose": 0.0,
            "time_horizon": 168.0,
            "tumor_growth_rate": 0.024,
            "healthy_growth_rate": 0.014,
            "interaction_strength": 0.1,
        },
        "schedule": {
            "lenalidomide_dose": {"type": "days_on_off", "on": 21, "off": 7},
            "bortezomib_dose": {"type": "days_list", "days": [1, 2, 8, 9, 15, 16]},
            "daratumumab_dose": {"type": "days_list", "days": []},
        },
    },
}
