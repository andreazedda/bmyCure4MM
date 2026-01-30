# Endpoints (quick reference)

=== "IT"
    Questa pagina elenca i principali endpoint web “umani” (pagine) e dove cercarli nel codice.

    ## Terminologia

    - **URLConf**: file `urls.py` che mappa path → view.
    - **View**: funzione/classe Django che gestisce la request.
    - **Template**: file HTML in `templates/` o `app/templates/...`.
    - **HTMX**: request parziali (snippet) per UI reattiva.

    ## Se cerchi…

    | Voglio… | Vai a… |
    | --- | --- |
    | routing principale | `mmportal/urls.py` + `*/urls.py` |
    | mapping DB/Models | `Reference → Database Objects` |
    | come partire | `Guides → Learning Path` |

    !!! tip "Approccio"
        Se vuoi trovare velocemente “dove vive” una pagina:
        1) apri `mmportal/urls.py` e gli `urls.py` delle app  
        2) cerca la view corrispondente  
        3) segui il template in `templates/` o `app/templates/...`

    ## Core

    | Area | Path | Codice |
    | --- | --- | --- |
    | Home | `/` | `mmportal/urls.py` |
    | Admin | `/admin/` | Django admin |

    ## Clinic (pazienti)

    | Pagina | Path (tipico) | URLConf | Views | Templates |
    | --- | --- | --- | --- | --- |
    | Lista pazienti | `/patients/` | `clinic/urls.py` | `clinic/views.py` | `clinic/templates/clinic/` |
    | Nuovo paziente | `/patients/new/` | `clinic/urls.py` | `clinic/views.py` | `clinic/templates/clinic/` |
    | Dettaglio | `/patients/<id>/` | `clinic/urls.py` | `clinic/views.py` | `clinic/templates/clinic/` |

    ## Simulator (scenari e simulazioni)

    | Pagina | Path (tipico) | URLConf | Views | Templates |
    | --- | --- | --- | --- | --- |
    | Lista scenari | `/simulator/` | `simulator/urls.py` | `simulator/views.py` | `simulator/templates/simulator/` |
    | Dettaglio scenario | `/simulator/scenario/<id>/` | `simulator/urls.py` | `simulator/views.py` | `simulator/templates/simulator/` |
    | Run simulation | (HTMX/POST) | `simulator/urls.py` | `simulator/views.py` | partials HTMX |

    ## ChemTools

    | Pagina | Path (tipico) | URLConf | Views | Templates |
    | --- | --- | --- | --- | --- |
    | Tools home | `/chem/` | `chemtools/urls.py` | `chemtools/views.py` | `chemtools/templates/chemtools/` |
    | Job detail | `/chem/jobs/<id>/` | `chemtools/urls.py` | `chemtools/views.py` | `chemtools/templates/chemtools/` |

    ## Docs Viewer (interno Django)

    | Pagina | Path | URLConf | Views |
    | --- | --- | --- | --- |
    | Index | `/docs/` | `docs_viewer/urls.py` | `docs_viewer/views.py` |
    | View | `/docs/view/<path>/` | `docs_viewer/urls.py` | `docs_viewer/views.py` |
    | Search | `/docs/search/` | `docs_viewer/urls.py` | `docs_viewer/views.py` |

=== "EN"
    This page lists the main human-facing endpoints and where they live in code.

    ## Terminology

    - **URLConf**: `urls.py` mapping path → view.
    - **View**: Django function/class handling the request.
    - **Template**: HTML file under `templates/` or `app/templates/...`.
    - **HTMX**: partial requests (snippets) for reactive UI.

    ## Quick find

    | I want… | Go to… |
    | --- | --- |
    | main routing | `mmportal/urls.py` + `*/urls.py` |
    | DB/Model mapping | `Reference → Database Objects` |
    | where to start | `Guides → Learning Path` |

    ## Core

    | Area | Path | Code |
    | --- | --- | --- |
    | Home | `/` | `mmportal/urls.py` |
    | Admin | `/admin/` | Django admin |

    ## Clinic

    | Page | Typical path | URLConf | Views | Templates |
    | --- | --- | --- | --- | --- |
    | Patients list | `/patients/` | `clinic/urls.py` | `clinic/views.py` | `clinic/templates/clinic/` |
    | New patient | `/patients/new/` | `clinic/urls.py` | `clinic/views.py` | `clinic/templates/clinic/` |
    | Detail | `/patients/<id>/` | `clinic/urls.py` | `clinic/views.py` | `clinic/templates/clinic/` |

    ## Simulator

    | Page | Typical path | URLConf | Views | Templates |
    | --- | --- | --- | --- | --- |
    | Scenarios list | `/simulator/` | `simulator/urls.py` | `simulator/views.py` | `simulator/templates/simulator/` |
    | Scenario detail | `/simulator/scenario/<id>/` | `simulator/urls.py` | `simulator/views.py` | `simulator/templates/simulator/` |
    | Run simulation | (HTMX/POST) | `simulator/urls.py` | `simulator/views.py` | HTMX partials |

    ## ChemTools

    | Page | Typical path | URLConf | Views | Templates |
    | --- | --- | --- | --- | --- |
    | Tools home | `/chem/` | `chemtools/urls.py` | `chemtools/views.py` | `chemtools/templates/chemtools/` |
    | Job detail | `/chem/jobs/<id>/` | `chemtools/urls.py` | `chemtools/views.py` | `chemtools/templates/chemtools/` |

    ## Docs Viewer

    | Page | Path | URLConf | Views |
    | --- | --- | --- | --- |
    | Index | `/docs/` | `docs_viewer/urls.py` | `docs_viewer/views.py` |
    | View | `/docs/view/<path>/` | `docs_viewer/urls.py` | `docs_viewer/views.py` |
    | Search | `/docs/search/` | `docs_viewer/urls.py` | `docs_viewer/views.py` |
