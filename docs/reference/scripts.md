# Scripts (utility)

=== "IT"
    ## `update_template.py`

    ## Se cerchi…

    | Voglio… | Vai a… |
    | --- | --- |
    | template coinvolti | `clinic/templates/clinic/` |
    | vedere dove appare in UI | `Reference → Endpoints` |

    ## Terminologia

    - **Marker**: stringa/anchor nel template usata per delimitare la sezione da sostituire.
    - **Nesting**: profondità `if/endif` da preservare durante il replace.

    Scopo: sostituire la sezione “simulation attempt” dentro `clinic/templates/clinic/patient_detail.html` usando il contenuto di `clinic/templates/clinic/patient_detail_decision_support.html`.

    Workflow (alto livello):

    ```mermaid
    flowchart TD
      A[Legge patient_detail_decision_support.html] --> B[Legge patient_detail.html]
      B --> C["Cerca marker: if latest_simulation_attempt"]
      C --> D["Conta nesting: if / endif"]
      D --> E[Sostituisce la sezione]
      E --> F[Scrive patient_detail.html]
    ```

    Esecuzione:

    ```bash
    python3 update_template.py
    ```

    !!! warning "Nota"
        Lo script stampa “Backup saved at … .backup” ma **non crea** effettivamente un backup. Se vuoi, posso correggerlo aggiungendo copia preventiva del file.

=== "EN"
    ## `update_template.py`

    ## Quick find

    | I want… | Go to… |
    | --- | --- |
    | involved templates | `clinic/templates/clinic/` |
    | where it appears in UI | `Reference → Endpoints` |

    ## Terminology

    - **Marker**: template anchor string delimiting the section to replace.
    - **Nesting**: `if/endif` depth that must be preserved during replacement.

    Purpose: replace the “simulation attempt” section inside `clinic/templates/clinic/patient_detail.html` using the contents of `clinic/templates/clinic/patient_detail_decision_support.html`.

    High-level workflow:

    ```mermaid
    flowchart TD
      A[Read patient_detail_decision_support.html] --> B[Read patient_detail.html]
      B --> C["Find marker: if latest_simulation_attempt"]
      C --> D["Count nesting: if / endif"]
      D --> E[Replace the section]
      E --> F[Write patient_detail.html]
    ```

    Run:

    ```bash
    python3 update_template.py
    ```

    !!! warning "Note"
        The script prints “Backup saved at … .backup” but does **not** actually create a backup. If you want, we can fix it by copying the file before writing.
