# Lipinski Analyzer Module

## Overview

This module retrieves a compound from PubChem using its CID and
calculates basic Lipinski descriptors with **RDKit**. Parameters are
printed to the console and saved to a log file.

## Usage

1. Edit `lipinski_analyzer.yaml` to set the PubChem `drug_id` or provide
   a `smiles` string.
2. Run `python lipinski_analyzer.py`.

A log file `lipinski_analyzer.log` will be created in the same folder.
