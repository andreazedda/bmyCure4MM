"""Utility functions for ligand similarity searching."""

import time
from typing import List, Dict
import requests
from urllib.parse import quote
import logging
import contextlib
import io

try:
    import pubchempy as pcp  # type: ignore
except Exception:  # pragma: no cover
    pcp = None


def get_fingerprint(smiles: str, fingerprint_type: str = "morgan", radius: int = 2):
    """Return a fingerprint for the given SMILES string."""
    try:
        with contextlib.redirect_stderr(io.StringIO()), contextlib.redirect_stdout(io.StringIO()):
            from rdkit import Chem  # type: ignore
            from rdkit.Chem import AllChem, MACCSkeys  # type: ignore
    except Exception as exc:  # pragma: no cover
        raise RuntimeError("RDKit is required for fingerprint computation but is not available") from exc

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
    try:
        with contextlib.redirect_stderr(io.StringIO()), contextlib.redirect_stdout(io.StringIO()):
            from rdkit.DataStructs import FingerprintSimilarity  # type: ignore
    except Exception as exc:  # pragma: no cover
        raise RuntimeError("RDKit is required for similarity computation but is not available") from exc
    return FingerprintSimilarity(fp1, fp2)


def search_similar_compounds(
    smiles: str,
    fingerprint_type: str = "morgan",
    radius: int = 2,
    n_results: int = 10,
    threshold: int = 90,
    max_records: int | None = None,
    max_retries: int = 3,
    include_properties: bool = False,
) -> List[Dict[str, object]]:
    """Search PubChem for compounds similar to the given SMILES.

    By default this returns only CIDs (no extra API calls).
    Set `include_properties=True` to fetch IsomericSMILES and compute RDKit similarities.
    """

    max_records = int(max_records) if max_records is not None else int(n_results) * 5
    threshold = int(threshold)

    # NOTE: tests assert query string is embedded in the URL argument (not passed via params)
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/"
        f"fastsimilarity_2d/smiles/{quote(smiles, safe='')}/cids/JSON"
        f"?Threshold={threshold}&MaxRecords={max_records}"
    )
    
    # Retry logic for PubChem API
    cids = []
    last_exception = None
    for attempt in range(max_retries):
        try:
            logging.info(f"Attempting PubChem similarity search (attempt {attempt + 1}/{max_retries})...")
            logging.info(f"URL: {url}")
            response = requests.get(url, timeout=30)
            response.raise_for_status()  # may raise in tests as generic Exception
            
            data = response.json()
            cids = data.get('IdentifierList', {}).get('CID', [])
            logging.info(f"Found {len(cids)} similar compounds from PubChem")
            break  # Success! Exit retry loop
            
        except Exception as exc:
            last_exception = exc
            status_code = getattr(getattr(exc, "response", None), "status_code", "unknown")
            logging.warning(
                f"PubChem API error {status_code} (attempt {attempt + 1}/{max_retries}): {exc}"
            )
            if attempt < max_retries - 1:
                wait_time = (attempt + 1) * 2  # Exponential backoff: 2s, 4s, 6s
                logging.info(f"Waiting {wait_time}s before retry...")
                time.sleep(wait_time)
            else:
                logging.error(f"PubChem search failed after {max_retries} attempts")
                raise RuntimeError(
                    f"PubChem API error after {max_retries} attempts (HTTP {status_code}). "
                    f"This may be a temporary issue or invalid input. "
                    f"Last error: {exc}"
                ) from exc

    if not cids:
        logging.warning("No similar compounds found in PubChem")
        return []

    if not include_properties:
        return [
            {"cid": cid, "smiles": None, "similarity": None}
            for cid in cids[:n_results]
        ]

    target_fp = None
    try:
        target_fp = get_fingerprint(smiles, fingerprint_type, radius)
    except Exception as exc:  # pragma: no cover
        logging.warning(
            "RDKit unavailable; returning properties without similarity scores (%s)",
            exc,
        )
    
    # Fetch compound details (SMILES) using PubChem properties API
    # pubchempy.get_compounds() doesn't populate SMILES by default, so we use the properties endpoint
    logging.info(f"Fetching SMILES for {len(cids)} compounds...")
    hits = []
    
    try:
        # Use PubChem properties API to get SMILES efficiently
        # https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/1,2,3/property/IsomericSMILES/JSON
        cids_to_fetch = cids[: max(n_results * 2, 1)]
        cid_str = ','.join(map(str, cids_to_fetch))
        props_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid_str}/property/IsomericSMILES/JSON"
        
        logging.info(f"Fetching properties from: {props_url}")
        response = requests.get(props_url, timeout=30)
        response.raise_for_status()
        
        props_data = response.json()
        properties = props_data.get('PropertyTable', {}).get('Properties', [])
        
        logging.info(f"Retrieved properties for {len(properties)} compounds")
        
        for prop in properties:
            cid = prop.get('CID')
            # Try IsomericSMILES first, fall back to SMILES
            comp_smiles = prop.get('IsomericSMILES') or prop.get('SMILES')
            
            if not comp_smiles:
                logging.warning(f"No SMILES for CID {cid}")
                continue
                
            sim = None
            if target_fp is not None:
                try:
                    fp = get_fingerprint(comp_smiles, fingerprint_type, radius)
                    sim = tanimoto_similarity(target_fp, fp)
                except Exception as exc:  # pragma: no cover
                    logging.warning("Failed to compute similarity for CID %s: %s", cid, exc)
                    sim = None
            hits.append({"cid": cid, "smiles": comp_smiles, "similarity": sim})
                
    except Exception as exc:
        logging.error(f"Error fetching compound properties: {exc}")
        # If we can't get details, return CIDs without similarity scores
        for cid in cids[:n_results]:
            hits.append({"cid": cid, "smiles": None, "similarity": None})
    
    hits.sort(key=lambda x: x["similarity"] if x["similarity"] is not None else 0, reverse=True)
    return hits[:n_results]


def compute_descriptors(smiles: str) -> Dict[str, float]:
    """Calculate common molecular descriptors for a SMILES string."""
    try:
        with contextlib.redirect_stderr(io.StringIO()), contextlib.redirect_stdout(io.StringIO()):
            from rdkit import Chem  # type: ignore
            from rdkit.Chem import Descriptors, rdMolDescriptors  # type: ignore
    except Exception as exc:  # pragma: no cover
        raise RuntimeError("RDKit is required for descriptor calculation but is not available") from exc

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
