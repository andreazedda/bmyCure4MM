# LigandSimilaritySearcher

This module searches PubChem for molecules similar to an input SMILES string.
It ranks candidates by Tanimoto similarity and calculates basic descriptors
for the top results.

## Usage

1. Adjust the configuration in `configs/ligand_similarity_searcher.yaml`.
2. Run the searcher from the sources directory:
   ```bash
   python ligand_similarity_searcher.py
   ```
3. Results are written to the CSV specified in the configuration.

## Outputs

The script creates a CSV containing:
- PubChem CID and SMILES
- Tanimoto similarity to the query
- Molecular weight, LogP, polar surface area, H‑bond acceptors and donors
