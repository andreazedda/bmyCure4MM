"""Utility functions for ligand similarity searching."""

import time
from typing import List, Dict
from rdkit import Chem
from rdkit.Chem import AllChem, MACCSkeys, Descriptors, rdMolDescriptors
from rdkit.DataStructs import FingerprintSimilarity
import pubchempy as pcp
import requests
from urllib.parse import quote
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
    max_retries: int = 3,
) -> List[Dict[str, object]]:
    """Search PubChem for compounds similar to the given SMILES using fastsimilarity_2d API."""
    target_fp = get_fingerprint(smiles, fingerprint_type, radius)
    
    # Use the correct PubChem REST API endpoint for similarity search
    # According to docs: https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/smiles/{SMILES}/cids/JSON
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/fastsimilarity_2d/smiles/{quote(smiles, safe='')}/cids/JSON"
    params = {
        'Threshold': 90,  # 90% similarity threshold
        'MaxRecords': n_results * 5  # Get extra to filter by Tanimoto
    }
    
    # Retry logic for PubChem API
    cids = []
    last_exception = None
    for attempt in range(max_retries):
        try:
            logging.info(f"Attempting PubChem similarity search (attempt {attempt + 1}/{max_retries})...")
            logging.info(f"URL: {url}")
            logging.info(f"Params: {params}")
            
            response = requests.get(url, params=params, timeout=30)
            response.raise_for_status()  # Raise exception for 4xx/5xx status codes
            
            data = response.json()
            cids = data.get('IdentifierList', {}).get('CID', [])
            logging.info(f"Found {len(cids)} similar compounds from PubChem")
            break  # Success! Exit retry loop
            
        except requests.exceptions.HTTPError as exc:
            last_exception = exc
            status_code = exc.response.status_code if exc.response else 'unknown'
            logging.warning(f"PubChem API HTTP error {status_code} (attempt {attempt + 1}/{max_retries}): {exc}")
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
                
        except (requests.exceptions.Timeout, requests.exceptions.ConnectionError) as exc:
            last_exception = exc
            logging.warning(f"PubChem API network error (attempt {attempt + 1}/{max_retries}): {exc}")
            if attempt < max_retries - 1:
                wait_time = (attempt + 1) * 2
                logging.info(f"Waiting {wait_time}s before retry...")
                time.sleep(wait_time)
            else:
                logging.error(f"PubChem search failed after {max_retries} attempts")
                raise RuntimeError(
                    f"PubChem API unavailable after {max_retries} attempts. "
                    f"Network error: {exc}"
                ) from exc
                
        except Exception as exc:
            logging.error(f"Unexpected error in PubChem search: {exc}")
            raise

    # Now fetch compound details for the CIDs we found
    if not cids:
        logging.warning("No similar compounds found in PubChem")
        return []
    
    # Fetch compound details (SMILES) using PubChem properties API
    # pubchempy.get_compounds() doesn't populate SMILES by default, so we use the properties endpoint
    logging.info(f"Fetching SMILES for {len(cids)} compounds...")
    hits = []
    
    try:
        # Use PubChem properties API to get SMILES efficiently
        # https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/1,2,3/property/IsomericSMILES/JSON
        cids_to_fetch = cids[:n_results * 2]
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
                
            try:
                fp = get_fingerprint(comp_smiles, fingerprint_type, radius)
                sim = tanimoto_similarity(target_fp, fp)
                hits.append({"cid": cid, "smiles": comp_smiles, "similarity": sim})
            except Exception as exc:  # skip bad molecules
                logging.warning(f"Failed to process CID {cid}: {exc}")
                
    except Exception as exc:
        logging.error(f"Error fetching compound properties: {exc}")
        # If we can't get details, return CIDs without similarity scores
        for cid in cids[:n_results]:
            hits.append({"cid": cid, "smiles": None, "similarity": None})
    
    hits.sort(key=lambda x: x["similarity"] if x["similarity"] is not None else 0, reverse=True)
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
