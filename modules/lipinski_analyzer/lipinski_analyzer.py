import os
import logging
import traceback
from colorama import Fore, Style
import requests
from rdkit import Chem
from rdkit.Chem import Descriptors
import yaml


def fetch_smiles(drug_id: str) -> str:
    """Fetch SMILES string for a PubChem CID."""
    url = (
        f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{drug_id}/"
        "property/CanonicalSMILES/JSON"
    )
    r = requests.get(url, timeout=10)
    r.raise_for_status()
    return r.json()["PropertyTable"]["Properties"][0]["CanonicalSMILES"]


def calculate_lipinski(smiles: str) -> dict:
    """Compute basic Lipinski descriptors from a SMILES string."""
    mol = Chem.MolFromSmiles(smiles)
    return {
        "MolecularWeight": Descriptors.MolWt(mol),
        "LogP": Descriptors.MolLogP(mol),
        "NumHDonors": Descriptors.NumHDonors(mol),
        "NumHAcceptors": Descriptors.NumHAcceptors(mol),
    }


def main() -> None:
    module_name = os.path.splitext(os.path.basename(__file__))[0]
    config_path = os.path.join(os.path.dirname(__file__), f"{module_name}.yaml")
    with open(config_path, "r") as fh:
        config = yaml.safe_load(fh)

    logging.basicConfig(
        filename=os.path.join(os.path.dirname(__file__), "lipinski_analyzer.log"),
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    try:
        smiles = config.get("smiles")
        if not smiles:
            smiles = fetch_smiles(config["drug_id"])
        params = calculate_lipinski(smiles)
        for k, v in params.items():
            print(f"{Fore.BLUE}{k}: {v}{Style.RESET_ALL}")
        logging.info("Parameters: %s", params)
    except Exception as exc:
        logging.error("Error: %s", exc)
        print(Fore.RED + "An error occurred" + Style.RESET_ALL)
        print(Fore.YELLOW + traceback.format_exc() + Style.RESET_ALL)


if __name__ == "__main__":
    main()
