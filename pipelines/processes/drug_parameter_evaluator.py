"""
    Description: This script fetches drug data from PubChem, calculates drug-related parameters, and visualizes the molecular structure of the drug.

    Functions:
    - fetch_drug_data: Fetches SMILES data for the given drug ID from PubChem or another database.
    - calculate_drug_parameters: Calculates drug-related parameters using the SMILES string.
    - visualize_drug_structure: Visualizes the molecular structure of a drug from SMILES.
    - main: Main function to execute drug data fetching, parameter calculation, and visualization.

    Returns:
    - Drug data fetched from PubChem.
    - Calculated drug parameters.
    - Visualization of the drug structure saved as an HTML file.
    
    Example:
    - Fetch drug data for a drug ID.
    - Calculate parameters for the drug.
    - Visualize the molecular structure of the drug.
    
    Usage:
    - Run the script.
    - Check the log file for detailed information.
    - Open the generated HTML file to view the drug structure.
    
    Configurations:
    - drug_id: The PubChem CID of the drug to evaluate.
    - viewer:
        - width: Width of the visualization window.
        - height: Height of the visualization window.
        
    Input data:
    - PubChem CID of the drug.
    
    Output:
    - Drug data fetched from PubChem.
    - Calculated drug parameters.
    - Visualization of the drug structure saved as an HTML file.
    
    Author: Andrea Zedda
    Last modified: 15/10/2024
"""

import py3Dmol
import requests
import logging
from colorama import Fore, Style
import traceback
import yaml
import json
import os
import time
from rdkit import Chem
from rdkit.Chem import Descriptors
import processes_utils as pu

# Load general settings
general_settings = pu.load_general_settings()

# Load configuration settings
module_name = os.path.splitext(os.path.basename(__file__))[0]
with open(os.path.join(general_settings['configs_path'], module_name + ".yaml"), "r") as config_file:
    config = yaml.safe_load(config_file)

PUBCHEM_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
HDRS = {
    "User-Agent": "MM-Portal/1.0 (+https://example.org) Python-requests",
    "Accept": "application/json",
    "Accept-Encoding": "gzip, deflate",
}

def _get_json(url: str, *, tries: int = 3, timeout: float = 10.0) -> dict:
    last = None
    for i in range(tries):
        try:
            r = requests.get(url, headers=HDRS, timeout=timeout)
            r.raise_for_status()
            return r.json()
        except requests.RequestException as e:
            last = e
            time.sleep(0.6 * (i + 1))
    raise last  # type: ignore[misc]


def fetch_drug_data(drug_id: str) -> str:
    """
    Return CanonicalSMILES for a PubChem CID.
    Fallbacks: (1) property Canonical+Isomeric, (2) record JSON parsing, (3) TXT endpoint.
    """
    cid = str(drug_id).strip()
    # 1) canonical + isomeric in one call
    url = f"{PUBCHEM_BASE}/compound/cid/{cid}/property/CanonicalSMILES,IsomericSMILES/JSON"
    try:
        payload = _get_json(url)
        props = payload.get("PropertyTable", {}).get("Properties", [])
        if props:
            rec = props[0]
            smiles = rec.get("CanonicalSMILES") or rec.get("IsomericSMILES")
            if smiles:
                return smiles
    except Exception:
        pass

    # 2) full record JSON (PC_Compounds)
    try:
        rec = _get_json(f"{PUBCHEM_BASE}/compound/cid/{cid}/JSON")
        comps = rec.get("PC_Compounds") or []
        # walk for a SMILES string in computed props
        for c in comps:
            for prop in (c.get("props") or []):
                urn = prop.get("urn", {})
                if urn.get("label") == "SMILES" and urn.get("name") in ("Canonical", "Isomeric"):
                    val = prop.get("value", {})
                    if "sval" in val and val["sval"]:
                        return val["sval"]
    except Exception:
        pass

    # 3) plain text as last resort
    try:
        txt = requests.get(
            f"{PUBCHEM_BASE}/compound/cid/{cid}/property/CanonicalSMILES/TXT",
            headers=HDRS, timeout=10.0
        )
        txt.raise_for_status()
        line = (txt.text or "").strip()
        if line:
            return line
    except Exception:
        pass

    raise ValueError(
        f"Could not retrieve SMILES for PubChem CID {cid}. "
        "Try entering a SMILES manually."
    )


def calculate_drug_parameters(smiles):
    """
    Calculates drug-related parameters using the SMILES string.
    """
    try:
        molecule = Chem.MolFromSmiles(smiles)

        # Calculate parameters
        parameters = {
            'Molecular_Weight': Descriptors.MolWt(molecule),
            'LogP': Descriptors.MolLogP(molecule),
            'Num_H_Donors': Descriptors.NumHDonors(molecule),
            'Num_H_Acceptors': Descriptors.NumHAcceptors(molecule),
            'TPSA': Descriptors.TPSA(molecule),
            'Num_Rotatable_Bonds': Descriptors.NumRotatableBonds(molecule),
            'LogS': -Descriptors.MolLogP(molecule) + 0.5  # Approximate LogS using empirical formula
        }

        # Print parameters with details
        print(Fore.GREEN + "Calculated Parameters:" + Style.RESET_ALL)
        details = {
            'Molecular_Weight': "Sum of atomic masses of all atoms in the molecule. Optimal range: ≤ 500 u.",
            'LogP': "Logarithm of the partition coefficient (octanol/water). Optimal range: 0–3.",
            'Num_H_Donors': "Number of hydrogen bond donors (e.g., –OH or –NH groups). Optimal range: ≤ 5.",
            'Num_H_Acceptors': "Number of hydrogen bond acceptors (e.g., oxygen, nitrogen). Optimal range: ≤ 10.",
            'TPSA': "Topological Polar Surface Area, related to drug permeability. Optimal range: ≤ 140 Å².",
            'Num_Rotatable_Bonds': "Number of rotatable bonds, related to molecular flexibility. Optimal range: ≤ 10.",
            'LogS': "Approximate aqueous solubility. Optimal range: ≥ -5."
        }

        for param, value in parameters.items():
            print(f"{Fore.BLUE}{param}: {value}{Style.RESET_ALL}")
            print(f"  {Fore.YELLOW}{details[param]}{Style.RESET_ALL}\n")

        logging.info("Parameters calculated: %s", parameters)
        return parameters
    except Exception as error:
        logging.error("Error calculating parameters: %s", error)
        raise



def visualize_drug_structure(smiles, width, height, output_html, parameters=None):
    """
    Visualizes the molecular structure of a drug from SMILES and displays calculated parameters.
    """
    logging.info("Initializing 3Dmol viewer for drug visualization.")
    mol_block = Chem.MolToMolBlock(Chem.MolFromSmiles(smiles))

    viewer = py3Dmol.view(width=width, height=height)
    viewer.addModel(mol_block, 'mol')
    viewer.setStyle({'stick': {}})
    viewer.zoomTo()
    viewer.setBackgroundColor('white')
    viewer.render()

    # Create HTML with both 3D viewer and parameters table
    viewer_html = viewer._make_html()
    
    # Build parameters table HTML if parameters provided
    params_html = ""
    if parameters:
        params_html = """
<style>
    body { font-family: Arial, sans-serif; margin: 20px; background: #f5f5f5; }
    .container { max-width: 1200px; margin: 0 auto; background: white; padding: 20px; border-radius: 8px; box-shadow: 0 2px 4px rgba(0,0,0,0.1); }
    h1 { color: #2c3e50; border-bottom: 3px solid #3498db; padding-bottom: 10px; }
    h2 { color: #34495e; margin-top: 30px; }
    .params-table { width: 100%; border-collapse: collapse; margin: 20px 0; }
    .params-table th { background: #3498db; color: white; padding: 12px; text-align: left; font-weight: bold; }
    .params-table td { padding: 10px; border-bottom: 1px solid #ddd; }
    .params-table tr:hover { background: #f8f9fa; }
    .param-name { font-weight: bold; color: #2c3e50; width: 30%; }
    .param-value { color: #27ae60; font-weight: bold; width: 20%; }
    .param-desc { color: #7f8c8d; font-size: 0.9em; }
    .optimal { color: #27ae60; }
    .warning { color: #e67e22; }
    .danger { color: #e74c3c; }
    .viewer-container { margin: 20px 0; text-align: center; }
</style>
<div class="container">
    <h1>💊 Drug Parameter Evaluation</h1>
    
    <h2>Molecular Structure</h2>
    <div class="viewer-container">
        """ + viewer_html + """
    </div>
    
    <h2>Calculated Parameters</h2>
    <table class="params-table">
        <thead>
            <tr>
                <th>Parameter</th>
                <th>Value</th>
                <th>Description & Optimal Range</th>
            </tr>
        </thead>
        <tbody>
"""
        
        # Parameter details and optimal ranges
        param_info = {
            'Molecular_Weight': {
                'desc': "Sum of atomic masses of all atoms in the molecule.",
                'optimal': "≤ 500 Da (Lipinski's Rule)",
                'check': lambda v: v <= 500
            },
            'LogP': {
                'desc': "Logarithm of the partition coefficient (octanol/water), indicates lipophilicity.",
                'optimal': "0–5 (ideally 0–3 for good absorption)",
                'check': lambda v: 0 <= v <= 5
            },
            'Num_H_Donors': {
                'desc': "Number of hydrogen bond donors (e.g., –OH or –NH groups).",
                'optimal': "≤ 5 (Lipinski's Rule)",
                'check': lambda v: v <= 5
            },
            'Num_H_Acceptors': {
                'desc': "Number of hydrogen bond acceptors (e.g., oxygen, nitrogen).",
                'optimal': "≤ 10 (Lipinski's Rule)",
                'check': lambda v: v <= 10
            },
            'TPSA': {
                'desc': "Topological Polar Surface Area, related to drug permeability and blood-brain barrier penetration.",
                'optimal': "≤ 140 Å² for good oral bioavailability",
                'check': lambda v: v <= 140
            },
            'Num_Rotatable_Bonds': {
                'desc': "Number of rotatable bonds, related to molecular flexibility and oral bioavailability.",
                'optimal': "≤ 10 for good oral bioavailability",
                'check': lambda v: v <= 10
            },
            'LogS': {
                'desc': "Approximate aqueous solubility (logarithmic scale).",
                'optimal': "≥ -5 for good solubility",
                'check': lambda v: v >= -5
            }
        }
        
        for param, value in parameters.items():
            if param in param_info:
                info = param_info[param]
                is_optimal = info['check'](value)
                status_class = 'optimal' if is_optimal else 'warning'
                status_icon = '✓' if is_optimal else '⚠'
                
                # Format value
                if isinstance(value, float):
                    value_str = f"{value:.2f}"
                else:
                    value_str = str(value)
                
                params_html += f"""
            <tr>
                <td class="param-name">{param.replace('_', ' ')}</td>
                <td class="param-value {status_class}">{status_icon} {value_str}</td>
                <td class="param-desc">{info['desc']}<br><strong>Optimal:</strong> {info['optimal']}</td>
            </tr>
"""
        
        params_html += """
        </tbody>
    </table>
    
    <h2>Drug-likeness Assessment (Lipinski's Rule of Five)</h2>
    <p style="padding: 15px; background: #ecf0f1; border-left: 4px solid #3498db; margin: 20px 0;">
        <strong>Rule of Five:</strong> A drug-like molecule should have:<br>
        • Molecular Weight ≤ 500 Da<br>
        • LogP ≤ 5<br>
        • H-bond Donors ≤ 5<br>
        • H-bond Acceptors ≤ 10<br>
        <br>
        Violations of these rules may indicate poor oral bioavailability.
    </p>
</div>
"""
    else:
        # Fallback to simple viewer if no parameters
        params_html = viewer_html

    with open(output_html, 'w') as html_file:
        html_file.write(params_html)
    
    logging.info("Visualization saved to %s", output_html)
    print(Fore.GREEN + f"Visualization saved to {output_html}. Open this file in a browser to view the structure." + Style.RESET_ALL)


def main():
    """
    Main function to execute drug data fetching, parameter calculation, and visualization.
    """
    logging.basicConfig(
        filename=os.path.join(general_settings['logs_path'], "drug_evaluator.log"),
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )

    try:
        # Fetch SMILES data for drug
        smiles = fetch_drug_data(config['drug_id'])

        # Calculate drug parameters
        parameters = calculate_drug_parameters(smiles)
        print(Fore.GREEN + "Calculated Parameters:" + Style.RESET_ALL, parameters)

        # Visualize the structure with parameters
        output_html = os.path.join(general_settings['outputs_path'], f"{config['drug_id']}_structure.html")
        visualize_drug_structure(smiles, config['viewer']['width'], config['viewer']['height'], output_html, parameters)

    except Exception as error:
        logging.error("An error occurred in the main function: %s", error)
        print(Fore.RED + "An error occurred. Traceback is shown below:" + Style.RESET_ALL)
        print(Fore.YELLOW + traceback.format_exc() + Style.RESET_ALL)


if __name__ == "__main__":
    main()
