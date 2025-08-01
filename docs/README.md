# Documentation

This folder contains guides and references for the **bmyCure4MM** project.

## Contents

- [Project Overview](#project-overview)
- [Getting Started](#getting-started)
- [Pipeline Scripts](#pipeline-scripts)
- [Lab Notebooks](#lab-notebooks)

### Project Overview

`bmyCure4MM` provides utilities for drug discovery and structural biology
experiments related to multiple myeloma.  The repository includes Python
scripts for visualising protein structures, evaluating drug parameters and
managing file paths, together with several Jupyter notebooks used for
research experiments.

### Getting Started

1. **Clone the repository**
   ```bash
   git clone <repo-url>
   cd bmyCure4MM
   ```
2. **Install dependencies**
   Python 3.10+ is required.  All Python dependencies used in the notebooks
   and scripts can be installed with
   ```bash
   pip install -r lab/requirements.txt
   ```
3. **Generate local settings**
   The pipeline scripts expect a `pipelines/configs/general_settings.json`
   file describing output directories.  You can generate it automatically
   by running
   ```bash
   python pipelines/processes/settings_generator.py
   ```
   This will create the `data` folders and a JSON file with absolute paths
   inside `pipelines/configs`.

### Pipeline Scripts

The `pipelines/processes` directory contains small utilities:

- **binding_visualizer.py** – Download a PDB structure and produce an HTML
  visualisation with *py3Dmol*.
- **drug_parameter_evaluator.py** – Fetch a drug from PubChem, calculate
  common molecular descriptors using *RDKit* and render the molecule in 3‑D.
- **paths_generator.py** – Build the folder structure under `pipelines/data`
  and save the resulting paths in `pipelines/data/paths.json`.
- **settings_generator.py** – Helper used during setup to create
  `general_settings.json`.
- **pdb_ligand_cid_finder.py** – Simple utility to list ligands within a PDB file.
- **LigandSimilaritySearcher** – Search PubChem for similar ligands and compute descriptors.

Each script reads its configuration from YAML files in `pipelines/configs`.
Run them with `python <script_name.py>` and check the `pipelines/data/` folder
for outputs and logs.

### Lab Notebooks

The `lab` directory collects exploratory notebooks and prototypes.
They require the packages listed in `lab/requirements.txt` and can be opened
with Jupyter Lab or Jupyter Notebook after installing the environment.  A
short description of each file is available in
[lab_file_overview.md](lab_file_overview.md).
