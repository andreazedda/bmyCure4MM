# Documentation Standard (editoriale + tecnico)

Questo documento definisce uno **standard professionale** per mantenere la documentazione consistente, navigabile e “a prova di GitHub Pages”.

## Obiettivi di qualità

- **Bilingue (IT/EN)**: ogni contenuto “core” deve esistere in entrambe le lingue.
- **Scansionabile**: tabelle, bullets, callout, “quick links”.
- **Verificabile**: riferimenti a file/funzioni/oggetti DB; build `mkdocs --strict` sempre ok.
- **Visuale**: diagrammi Mermaid e figure SVG dove aiutano (con zoom/pan).
- **Onboarding**: percorsi guidati per persona (clinico/dev/research).

## Struttura pagina (template)

Ogni pagina guida/reference dovrebbe seguire questa struttura (quando applicabile):

1. **Cos’è / Scope** (1–3 frasi)
2. **Quick links / “Se cerchi…”** (tabella)
3. **Terminologia** (glossario locale)
4. **Concetto → Esempio → Implementazione** (sempre)
5. **Diagrammi / Figure**
6. **Riferimenti al codice** (file + simboli)
7. **Edge cases / Troubleshooting** (se utile)

## Bilinguismo

Preferenza: usare i tab di Material.

```markdown
=== "IT"
    Testo in italiano…

=== "EN"
    English content…
```

Regola: il contenuto deve essere **equivalente** (stessa informazione), anche se adattato per leggibilità.

## Tabelle (quando usarle)

Usare tabelle per:

- parametri (nome, tipo, default, unità, effetto)
- mapping (Model ↔ Table ↔ File)
- decision thresholds (range → badge → meaning)
- endpoint (path → view → template)

## Diagrammi (Mermaid)

Linee guida:

- label brevi; evitare caratteri ambigui dentro `[]` (`_`, `/`, backtick, ecc.) se causano parse error
- preferire `["label"]` con testo semplice o `<br/>` per righe
- diagrammi grandi: contare sul **click-to-zoom** già attivo

## Figure (SVG)

Per grafici/plot:

- preferire SVG generati e committati in `docs/assets/images/...`
- mantenere script di generazione in `docs/assets/scripts/`
- non dipendere da CDN o runtime Python in GitHub Pages

## Checklist pre-merge

- `./venv/bin/mkdocs build --strict` passa
- link interni ok (no warning)
- nuove pagine in `mkdocs.yml nav`
- IT/EN presenti per pagine core
- diagrammi Mermaid renderizzano (niente “Syntax error in text”)

