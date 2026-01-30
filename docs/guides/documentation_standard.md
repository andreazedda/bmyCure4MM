# Documentation Standard (editorial + technical)

=== "IT"
    Questo documento definisce uno **standard professionale** per mantenere la documentazione consistente, navigabile e “a prova di GitHub Pages”.

    ## Se cerchi…

    | Voglio… | Vai a… |
    | --- | --- |
    | roadmap di miglioramento | `Guides → Documentation Roadmap` |
    | struttura consigliata pagine | sezione “Template pagina” qui sotto |
    | checklist di qualità | sezione “Checklist” |

    ## Terminologia

    - **Core page**: pagina guida/reference che deve essere IT/EN.
    - **Strict build**: `mkdocs build --strict` senza warning.
    - **As-built**: documentazione che descrive lo stato reale del codice.

    ## Obiettivi di qualità

    - **Bilingue (IT/EN)**: ogni contenuto “core” deve esistere in entrambe le lingue.
    - **Scansionabile**: tabelle, bullet, callout, “quick links”.
    - **Verificabile**: riferimenti a file/funzioni/oggetti DB; build `mkdocs --strict` sempre ok.
    - **Visuale**: diagrammi Mermaid e figure SVG dove aiutano (con zoom/pan).
    - **Onboarding**: percorsi guidati per persona (clinico/dev/research).

    ## Template pagina (uguale per tutte le guide/reference)

    Ogni pagina guida/reference dovrebbe seguire questa struttura (quando applicabile):

    1. **Scope** (1–3 frasi)
    2. **Se cerchi…** (tabella)
    3. **Terminologia** (glossario locale, 5–15 voci max)
    4. **Teoria** (concetti + formule, se utile)
    5. **Implementazione** (mapping al codice)
    6. **Esempi** (worked example / snippet / output atteso)
    7. **Troubleshooting** (solo se ci sono failure mode comuni)

    ### Blocchi standard (copiabili)

    ```markdown
    === "IT"
        Descrizione breve…

        ## Se cerchi…
        | Voglio… | Vai a… |
        | --- | --- |
        | … | … |

        ## Terminologia
        - …

    === "EN"
        Short description…

        ## Quick find
        | I want… | Go to… |
        | --- | --- |
        | … | … |

        ## Terminology
        - …
    ```

    ## Bilinguismo

    Preferenza: usare i tab di Material:

    - `=== "IT"` per italiano
    - `=== "EN"` per inglese

    Regola: il contenuto deve essere **equivalente** (stessa informazione), anche se adattato per leggibilità.

    ## Tabelle (quando usarle)

    Usare tabelle per:

    - parametri (nome, tipo, default, unità, effetto)
    - mapping (Model ↔ Table ↔ File)
    - decision thresholds (range → badge → meaning)
    - endpoint (path → view → template → payload)

    ## Diagrammi (Mermaid)

    Linee guida:

    - label brevi; evitare caratteri ambigui dentro `[]` (`_`, `/`, backtick, ecc.) se causano parse error
    - preferire `["label"]` con testo semplice o `<br/>` per righe
    - diagrammi grandi: contare sul **click-to-zoom** già attivo nel sito

    ## Figure (SVG)

    Per grafici/plot:

    - preferire SVG generati e committati in `docs/assets/images/...`
    - mantenere script di generazione in `docs/assets/scripts/`
    - non dipendere da runtime Python in GitHub Pages

    ## Checklist (pre-merge)

    - `./venv/bin/mkdocs build --strict` passa
    - link interni ok (no warning)
    - nuove pagine in `mkdocs.yml nav`
    - IT/EN presenti per pagine core
    - diagrammi Mermaid renderizzano (niente “Syntax error in text”)

=== "EN"
    This document defines a **professional standard** to keep the documentation consistent, navigable, and “GitHub Pages friendly”.

    ## Quick find

    | I want… | Go to… |
    | --- | --- |
    | improvement roadmap | `Guides → Documentation Roadmap` |
    | recommended page structure | see “Page template” below |
    | quality checklist | see “Checklist” |

    ## Terminology

    - **Core page**: guide/reference page that must exist in IT/EN.
    - **Strict build**: `mkdocs build --strict` with no warnings.
    - **As-built**: documentation describing the actual code state.

    ## Quality goals

    - **Bilingual (IT/EN)**: every core topic exists in both languages.
    - **Scannable**: tables, bullets, callouts, quick links.
    - **Verifiable**: references to files/functions/DB objects; `mkdocs --strict` stays green.
    - **Visual**: Mermaid diagrams and committed SVG figures (with zoom/pan).
    - **Onboarding**: guided paths per persona (clinician/dev/research).

    ## Page template (same for all guides/reference)

    Each guide/reference page should follow this structure (when applicable):

    1. **Scope** (1–3 sentences)
    2. **Quick find** (table)
    3. **Terminology** (local glossary, 5–15 items)
    4. **Theory** (concepts + equations when useful)
    5. **Implementation** (code mapping)
    6. **Examples** (worked example / snippet / expected output)
    7. **Troubleshooting** (only for common failure modes)

    ## Tables (when to use them)

    Use tables for:

    - parameters (name, type, default, units, effect)
    - mappings (Model ↔ Table ↔ File)
    - thresholds (range → badge → meaning)
    - endpoints (path → view → template → payload)

    ## Mermaid diagrams

    Guidelines:

    - keep labels short; avoid ambiguous characters inside `[]` when they trigger parse errors
    - prefer `["label"]` with plain text or `<br/>`
    - rely on the site’s click-to-zoom for large diagrams

    ## SVG figures

    For plots:

    - prefer committed SVGs under `docs/assets/images/...`
    - keep generators under `docs/assets/scripts/`
    - do not depend on runtime Python in GitHub Pages

    ## Checklist (pre-merge)

    - `./venv/bin/mkdocs build --strict` passes
    - internal links ok (no warnings)
    - new pages included in `mkdocs.yml nav`
    - IT/EN present for core pages
    - Mermaid renders (no “Syntax error in text”)
