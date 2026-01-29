# Endpoints (quick reference)

Questa pagina elenca i principali endpoint web “umani” (pagine) e dove cercarli nel codice.

!!! tip "Approccio"
    Se vuoi trovare velocemente “dove vive” una pagina:
    1) apri `mmportal/urls.py` e gli `urls.py` delle app  
    2) cerca la view corrispondente  
    3) segui il template in `templates/` o `app/templates/...`

## Core

- Home: `/` (vedi `mmportal/urls.py`)
- Admin: `/admin/`

## Clinic (pazienti)

Tipicamente:

- Lista pazienti: `/patients/`
- Nuovo paziente: `/patients/new/` (o simile)
- Dettaglio paziente: `/patients/<id>/`

Codice:
- URL: `clinic/urls.py`
- View: `clinic/views.py`
- Template: `clinic/templates/clinic/`

## Simulator (scenari e simulazioni)

Tipicamente:

- Lista scenari: `/simulator/`
- Dettaglio scenario: `/simulator/scenario/<id>/`
- Run simulation (POST/HTMX): endpoint interno nel dettaglio scenario

Codice:
- URL: `simulator/urls.py`
- View: `simulator/views.py`
- Template: `simulator/templates/simulator/`

## ChemTools

Tipicamente:

- Home tools: `/chem/` (o simile)
- Job detail: `/chem/jobs/<id>/`

Codice:
- URL: `chemtools/urls.py`
- View: `chemtools/views.py`
- Template: `chemtools/templates/chemtools/`

## Docs Viewer (interno Django)

Tipicamente:

- Index: `/docs/`
- View: `/docs/view/<path>/`
- Search: `/docs/search/`

Codice:
- URL: `docs_viewer/urls.py`
- View: `docs_viewer/views.py`

