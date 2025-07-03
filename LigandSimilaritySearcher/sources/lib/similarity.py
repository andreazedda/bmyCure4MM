"""Utility functions for ligand similarity searching."""

from typing import List, Dict
from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys, Descriptors, rdMolDescriptors
from rdkit.DataStructs import FingerprintSimilarity
import pubchempy as pcp
import logging


def get_fingerprint(smiles: str, fingerprint_type: str = "morgan", radius: int = 2):
    """Return a fingerprint for the given SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    if fingerprint_type == "morgan":
        return AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits=2048)
    if fingerprint_type == "maccs":
        return MACCSkeys.GenMACCSKeys(mol)
    raise ValueError(f"Unsupported fingerprint type: {fingerprint_type}")


def tanimoto_similarity(fp1, fp2) -> float:
    """Compute Tanimoto similarity between two fingerprints."""
    return FingerprintSimilarity(fp1, fp2)


def search_similar_compounds(
    smiles: str,
    fingerprint_type: str = "morgan",
    radius: int = 2,
    n_results: int = 10,
) -> List[Dict[str, object]]:
    """Search PubChem for compounds similar to the given SMILES."""
    target_fp = get_fingerprint(smiles, fingerprint_type, radius)
    try:
        raw = pcp.get_compounds(
            smiles,
            namespace="smiles",
            searchtype="similarity",
            listkey_count=n_results * 5,
        )
    except Exception as exc:
        logging.error("PubChem search failed: %s", exc)
        raise

    hits = []
    for comp in raw:
        comp_smiles = comp.isomeric_smiles
        if not comp_smiles:
            continue
        try:
            fp = get_fingerprint(comp_smiles, fingerprint_type, radius)
            sim = tanimoto_similarity(target_fp, fp)
            hits.append({"cid": comp.cid, "smiles": comp_smiles, "similarity": sim})
        except Exception as exc:  # skip bad molecules
            logging.warning("Failed to process CID %s: %s", comp.cid, exc)
    hits.sort(key=lambda x: x["similarity"], reverse=True)
    return hits[:n_results]


def compute_descriptors(smiles: str) -> Dict[str, float]:
    """Calculate common molecular descriptors for a SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES for descriptor calculation: {smiles}")
    return {
        "MW": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "TPSA": rdMolDescriptors.CalcTPSA(mol),
        "HBA": Descriptors.NumHAcceptors(mol),
        "HBD": Descriptors.NumHDonors(mol),
    }
