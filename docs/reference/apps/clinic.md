# App `clinic`

## Scopo

Gestisce:

- anagrafica paziente
- assessment (lab snapshot, IMWG)
- terapie (timeline, regimen)
- citogenetica (catalogo + storico)

## Modelli principali

Vedi `clinic/models.py`:

- `Patient`
- `Assessment`
- `Regimen`
- `PatientTherapy`
- `CytogeneticAbnormality`
- `PatientCytogenetics`

!!! note "Accesso minimo"
    `Patient.owner` è una FK a `AUTH_USER_MODEL` e implementa un “ownership” basilare per filtrare/limitare l’accesso (a seconda delle view).

## Template principali

Directory: `clinic/templates/clinic/`

Pagine chiave:
- `patient_list.html`, `patient_detail.html`
- `assessment_form.html`, `patient_form.html`
- `regimen_list.html`, `regimen_form.html`

