# Documentation Viewer

Sistema completo per visualizzare la documentazione markdown direttamente nel browser.

## 🎯 Caratteristiche

### Sicurezza
- ✅ **Whitelist validation**: Solo file autorizzati sono accessibili
- ✅ **Path traversal protection**: Protezione contro attacchi `../` 
- ✅ **Input sanitization**: Rimozione caratteri pericolosi (null bytes, etc.)
- ✅ **File type validation**: Solo file markdown dalle directory autorizzate

### Funzionalità
- 📖 **Rendering markdown completo**: Tabelle, codice, TOC automatico
- 🎨 **Syntax highlighting**: Codice colorato con Highlight.js
- 🔍 **Ricerca full-text**: Cerca in tutti i documenti
- 📊 **Analytics**: Tracciamento delle visualizzazioni
- ⭐ **Feedback utenti**: Sistema di valutazione con stelle
- 📥 **Download**: Scarica il markdown originale
- 🍞 **Breadcrumbs**: Navigazione gerarchica
- 📱 **Responsive**: Design ottimizzato per mobile

## 🚀 Utilizzo

### Accesso Web
Vai a: `http://localhost:8000/docs/`

### Whitelist File
I seguenti percorsi sono accessibili:
```python
ALLOWED_DOC_PATHS = [
    'docs/',                     # Tutta la directory docs/
    'tests/',                    # Tutti i test
    'README.md',                 # README principale
    'IMPLEMENTATION_LOG.md',     # Log implementazione
    'IMPLEMENTATION_SUMMARY.md', # Sommario implementazione
]
```

Per aggiungere nuovi percorsi, modifica `docs_viewer/utils.py`.

### Validazione Documentazione
Comando management per controlli qualità:

```bash
# Validazione base
python manage.py validate_docs

# Con dettagli completi
python manage.py validate_docs --verbose

# Con suggerimenti per link rotti
python manage.py validate_docs --fix-links
```

Il comando verifica:
- ✅ Rendering corretto di tutti i file markdown
- ✅ Link interni rotti
- ✅ File senza titolo (H1)
- ✅ Encoding UTF-8 valido
- ✅ Sintassi markdown corretta

## 🧪 Test

### Esecuzione Test
```bash
# Tutti i test
pytest docs_viewer/tests/ -v

# Solo test unitari
pytest docs_viewer/tests/test_utils.py -v

# Solo test di integrazione
pytest docs_viewer/tests/test_views.py -v

# Test E2E con Playwright (se disponibili)
pytest tests/e2e/test_docs_viewer.py --headed
```

### Coverage Test
- **40 test totali**: 100% pass rate
- **Test unitari**: 20 test (security, rendering, utilities)
- **Test integrazione**: 20 test (views, caching, search)
- **Test E2E**: Pianificati (Playwright)

## 📁 Struttura

```
docs_viewer/
├── models.py              # DocumentView, DocumentFeedback
├── views.py               # 4 views (index, view, raw, search)
├── utils.py               # Security + markdown rendering
├── urls.py                # URL patterns
├── admin.py               # Admin interface
├── templates/
│   └── docs_viewer/
│       ├── index.html     # Lista documenti
│       ├── view.html      # Viewer markdown
│       └── search.html    # Risultati ricerca
├── tests/
│   ├── test_utils.py      # Test security/utils
│   └── test_views.py      # Test views/integration
└── management/
    └── commands/
        └── validate_docs.py  # Comando validazione
```

## 🔒 Sicurezza

### Layers di Protezione
1. **Input Sanitization**: `sanitize_path()` rimuove caratteri pericolosi
2. **Whitelist Check**: `is_safe_path()` verifica path nella whitelist
3. **Path Validation**: Controlla pattern `../` e path assoluti
4. **File Existence**: Verifica esistenza file prima dell'accesso

### Esempio Attack Prevention
```python
# ❌ Bloccato - Path traversal
/docs/view/../../manage.py

# ❌ Bloccato - Absolute path
/docs/view//etc/passwd

# ❌ Bloccato - File non in whitelist
/docs/view/mmportal/settings.py

# ✅ Permesso - File in whitelist
/docs/view/README.md
/docs/view/docs/en/quickstart.md
```

## 🎨 Personalizzazione

### Aggiungere File alla Whitelist
Modifica `docs_viewer/utils.py`:

```python
ALLOWED_DOC_PATHS = [
    'docs/',
    'tests/',
    'README.md',
    'your_new_directory/',  # ← Aggiungi qui
]
```

### Cambiare Stile Syntax Highlighting
Modifica `docs_viewer/templates/docs_viewer/view.html`:

```html
<!-- Cambia tema highlight.js -->
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/highlight.js/11.9.0/styles/github-dark.min.css">
<!-- Opzioni: github, monokai, atom-one-dark, vs2015, etc. -->
```

### Configurare Cache
Modifica `docs_viewer/views.py`:

```python
@cache_page(60 * 15)  # ← Cambia durata cache (secondi)
def docs_index(request):
    ...
```

## 📊 Analytics

### Modelli Database
```python
# Tracciamento visualizzazioni
DocumentView(path, user, viewed_at)

# Feedback utenti
DocumentFeedback(document_path, user, rating, comment)
```

### Query Admin
Accedi all'interfaccia admin per vedere:
- Documenti più visualizzati
- Rating medio per documento
- Commenti utenti
- Pattern di navigazione

## 🐛 Troubleshooting

### Test Falliti
```bash
# Pulisci cache e riprova
python manage.py clear_cache
pytest docs_viewer/tests/ -v --tb=short
```

### Link Rotti
```bash
# Verifica e mostra link rotti
python manage.py validate_docs --fix-links
```

### File Non Visibili
1. Verifica che il file sia nella whitelist
2. Controlla permessi file (lettura)
3. Verifica encoding UTF-8
4. Esegui: `python manage.py validate_docs --verbose`

## 🔄 CI/CD Integration

### GitHub Actions Example
```yaml
- name: Validate Documentation
  run: |
    python manage.py validate_docs
    
- name: Run Docs Tests
  run: |
    pytest docs_viewer/tests/ -v --cov=docs_viewer
```

## 📝 Changelog

### v1.0.0 (2025-01-XX)
- ✅ Sistema completo documentazione viewer
- ✅ Security whitelist + path traversal protection
- ✅ Markdown rendering con TOC e syntax highlighting
- ✅ Ricerca full-text
- ✅ Analytics e feedback utenti
- ✅ 40 test con 100% pass rate
- ✅ Comando validazione documentazione
- ✅ Admin interface
- ✅ Responsive design bilingue (EN/IT)

## 🤝 Contribuire

Per aggiungere nuove funzionalità:
1. Aggiungi test in `docs_viewer/tests/`
2. Implementa feature
3. Esegui test: `pytest docs_viewer/tests/ -v`
4. Valida docs: `python manage.py validate_docs`
5. Commit e PR

## 📄 License

Parte del progetto bmyCure4MM.
