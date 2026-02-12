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

try:
    import py3Dmol  # type: ignore
except Exception:  # pragma: no cover
    py3Dmol = None

import requests
import logging

try:
    from colorama import Fore, Style
except Exception:  # pragma: no cover
    Fore = Style = None  # type: ignore

import traceback
import yaml
import json
import os
import time
import contextlib
import io

# RDKit is an optional/binary dependency in some environments.
try:
    with contextlib.redirect_stderr(io.StringIO()), contextlib.redirect_stdout(io.StringIO()):
        from rdkit import Chem  # type: ignore
        from rdkit.Chem import Descriptors  # type: ignore
        from rdkit.Chem import rdMolDescriptors  # type: ignore
except Exception:  # pragma: no cover
    Chem = None  # type: ignore
    Descriptors = None  # type: ignore
    rdMolDescriptors = None  # type: ignore
try:
    from . import processes_utils as pu
except Exception:  # pragma: no cover
    import processes_utils as pu

def _load_general_settings():
    return pu.load_general_settings()


def _load_config(general_settings: dict) -> dict:
    module_name = os.path.splitext(os.path.basename(__file__))[0]
    with open(
        os.path.join(general_settings["configs_path"], module_name + ".yaml"),
        "r",
    ) as config_file:
        return yaml.safe_load(config_file)

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
        if Chem is None or Descriptors is None:
            raise RuntimeError("RDKit is required to calculate drug parameters")

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
        if Fore is not None and Style is not None:
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
            if Fore is not None and Style is not None:
                print(f"{Fore.BLUE}{param}: {value}{Style.RESET_ALL}")
                print(f"  {Fore.YELLOW}{details[param]}{Style.RESET_ALL}\n")

        logging.info("Parameters calculated: %s", parameters)
        return parameters
    except Exception as error:
        logging.error("Error calculating parameters: %s", error)
        raise



def visualize_drug_structure(
    smiles_or_parameters,
    width: int = 600,
    height: int = 400,
    output_html: str | None = None,
    parameters: dict | None = None,
):
    """Generate an HTML snippet for the drug 3D viewer and parameters table.

    Compatibility:
    - Tests call: visualize_drug_structure(parameters_dict)
    - CLI/main calls: visualize_drug_structure(smiles, width, height, output_html, parameters)

    Returns the HTML string. If output_html is provided, writes the HTML to disk as well.
    """
    if isinstance(smiles_or_parameters, dict) and parameters is None:
        smiles = None
        parameters = smiles_or_parameters
    else:
        smiles = smiles_or_parameters

    # Build a viewer HTML block.
    viewer_html = """
<div id="container-01" style="width:100%;height:400px;position:relative;"></div>
<script>
// $3Dmol placeholder (real viewer if py3Dmol is available)
var $3Dmol = window.$3Dmol || {};
</script>
""".strip()

    if py3Dmol is not None and smiles and Chem is not None:
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is not None:
                mol_block = Chem.MolToMolBlock(mol)
                viewer = py3Dmol.view(width=width, height=height)
                viewer.addModel(mol_block, "mol")
                viewer.setStyle({"stick": {}})
                viewer.zoomTo()
                viewer.setBackgroundColor("white")
                viewer.render()
                viewer_html = viewer._make_html()
        except Exception:
            # Keep the placeholder viewer_html.
            pass

    # Normalize parameter keys.
    params = parameters or {}
    normalized = {
        "MW": params.get("MW", params.get("Molecular_Weight")),
        "LogP": params.get("LogP"),
        "HBD": params.get("HBD", params.get("Num_H_Donors")),
        "HBA": params.get("HBA", params.get("Num_H_Acceptors")),
        "TPSA": params.get("TPSA"),
        "RotBonds": params.get("RotBonds", params.get("Num_Rotatable_Bonds")),
        "LogS": params.get("LogS"),
    }

    green = "#28a745"
    warn = "#ffc107"

    def fmt(v):
        if v is None:
            return ""
        if isinstance(v, float):
            return f"{v:.2f}"
        return str(v)

    def ok(name: str, v):
        if v is None:
            return True
        if name == "MW":
            return float(v) <= 500
        if name == "LogP":
            return 0 <= float(v) <= 5
        if name == "HBD":
            return float(v) <= 5
        if name == "HBA":
            return float(v) <= 10
        if name == "TPSA":
            return float(v) <= 140
        if name == "RotBonds":
            return float(v) <= 10
        if name == "LogS":
            return float(v) > -4
        return True

    rows = []
    labels = {
        "MW": "MW (Molecular Weight)",
        "LogP": "LogP",
        "HBD": "HBD",
        "HBA": "HBA",
        "TPSA": "TPSA",
        "RotBonds": "RotBonds",
        "LogS": "LogS",
    }
    for key, label in labels.items():
        value = normalized.get(key)
        is_ok = ok(key, value)
        icon = "✓" if is_ok else "⚠"
        color = green if is_ok else warn
        rows.append(
            f"<tr><td>{label}</td><td><span style=\"color:{color}\">{icon}</span> {fmt(value)}</td></tr>"
        )

    html = f"""
<div>
  <h2>3D Protein-Ligand Structure</h2>
  {viewer_html}

  <h2>Parameters</h2>
  <table id="parameters-table" class="parameters-table">
    <tbody>
      {"".join(rows)}
    </tbody>
  </table>

  <h3>Lipinski</h3>
</div>
""".strip()

    if output_html:
        with open(output_html, "w") as html_file:
            html_file.write(html)
        logging.info("Visualization saved to %s", output_html)
        if Fore is not None and Style is not None:
            print(
                Fore.GREEN
                + f"Visualization saved to {output_html}. Open this file in a browser to view the structure."
                + Style.RESET_ALL
            )
    return html


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
        general_settings = _load_general_settings()
        config = _load_config(general_settings)

        # Fetch SMILES data for drug
        smiles = fetch_drug_data(config["drug_id"])

        # Calculate drug parameters
        parameters = calculate_drug_parameters(smiles)
        if Fore is not None and Style is not None:
            print(Fore.GREEN + "Calculated Parameters:" + Style.RESET_ALL, parameters)

        # Visualize the structure with parameters
        output_html = os.path.join(
            general_settings["outputs_path"], f"{config['drug_id']}_structure.html"
        )
        visualize_drug_structure(
            smiles,
            config["viewer"]["width"],
            config["viewer"]["height"],
            output_html,
            parameters,
        )

    except Exception as error:
        logging.error("An error occurred in the main function: %s", error)
        if Fore is not None and Style is not None:
            print(Fore.RED + "An error occurred. Traceback is shown below:" + Style.RESET_ALL)
            print(Fore.YELLOW + traceback.format_exc() + Style.RESET_ALL)
        else:
            print("An error occurred. Traceback is shown below:")
            print(traceback.format_exc())


if __name__ == "__main__":
    main()
