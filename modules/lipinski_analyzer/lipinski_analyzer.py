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
        "property/SMILES/JSON"
    )
    print(f"{Fore.CYAN}🔎🔍 [DEBUG] Fetching SMILES for CID: {drug_id} {Style.RESET_ALL}")
    print(f"{Fore.CYAN}🌐🌍 [DEBUG] Request URL: {url}{Style.RESET_ALL}")
    r = requests.get(url, timeout=10)
    print(f"{Fore.CYAN}📡 [DEBUG] HTTP status code: {r.status_code}{Style.RESET_ALL}")
    try:
        r.raise_for_status()
    except Exception as e:
        print(f"{Fore.RED}❌🚨 [DEBUG] PubChem response: {r.text}{Style.RESET_ALL}")
        print(f"{Fore.RED}🛑 [DEBUG] Exception: {e}{Style.RESET_ALL}")
        raise
    print(f"{Fore.GREEN}✅🟢 [DEBUG] PubChem response OK!{Style.RESET_ALL}")
    props = r.json()["PropertyTable"]["Properties"][0]
    smiles = props.get("SMILES")
    print(f"{Fore.CYAN}🧬 [DEBUG] SMILES found: {smiles}{Style.RESET_ALL}")
    return smiles


def calculate_lipinski(smiles: str) -> dict:
    """Compute basic Lipinski descriptors from a SMILES string."""
    print(f"{Fore.MAGENTA}🧪🔬 [DEBUG] Calculating Lipinski descriptors for SMILES: {smiles}{Style.RESET_ALL}")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"{Fore.RED}❌💥 [DEBUG] Invalid SMILES string!{Style.RESET_ALL}")
        return {}
    else:
        print(f"{Fore.GREEN}✅🧪 [DEBUG] RDKit molecule created!{Style.RESET_ALL}")
    mw = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hdon = Descriptors.NumHDonors(mol)
    haccept = Descriptors.NumHAcceptors(mol)
    print(f"{Fore.BLUE}📏 [DEBUG] Molecular Weight: {mw}{Style.RESET_ALL}")
    print(f"{Fore.BLUE}💧 [DEBUG] LogP: {logp}{Style.RESET_ALL}")
    print(f"{Fore.BLUE}🧲 [DEBUG] NumHDonors: {hdon}{Style.RESET_ALL}")
    print(f"{Fore.BLUE}🧲 [DEBUG] NumHAcceptors: {haccept}{Style.RESET_ALL}")
    return {
        "MolecularWeight": mw,
        "LogP": logp,
        "NumHDonors": hdon,
        "NumHAcceptors": haccept,
    }


def main() -> None:
    module_name = os.path.splitext(os.path.basename(__file__))[0]
    config_path = os.path.join(os.path.dirname(__file__), f"{module_name}.yaml")
    print(f"{Fore.YELLOW}📂🗂️ [DEBUG] Loading config from: {config_path}{Style.RESET_ALL}")
    with open(config_path, "r") as fh:
        config = yaml.safe_load(fh)
    print(f"{Fore.YELLOW}📝🖊️ [DEBUG] Config loaded: {config}{Style.RESET_ALL}")

    logging.basicConfig(
        filename=os.path.join(os.path.dirname(__file__), "lipinski_analyzer.log"),
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )
    print(f"{Fore.YELLOW}🗒️ [DEBUG] Logging to: lipinski_analyzer.log{Style.RESET_ALL}")

    try:
        print(f"{Fore.YELLOW}🚦 [DEBUG] Starting Lipinski analysis...{Style.RESET_ALL}")
        smiles = config.get("smiles")
        if not smiles:
            print(f"{Fore.CYAN}🔄🔁 [DEBUG] No SMILES in config, fetching from PubChem...{Style.RESET_ALL}")
            smiles = fetch_smiles(config["drug_id"])
        else:
            print(f"{Fore.CYAN}🧬🧫 [DEBUG] Using SMILES from config: {smiles}{Style.RESET_ALL}")
        params = calculate_lipinski(smiles)
        print(f"{Fore.GREEN}🎯🏆 [DEBUG] Lipinski parameters calculated: {params}{Style.RESET_ALL}")
        for k, v in params.items():
            print(f"{Fore.BLUE}🔹 {k}: {v}{Style.RESET_ALL}")
        logging.info("Parameters: %s", params)
        print(f"{Fore.GREEN}✅🎉 [DEBUG] Analysis complete!{Style.RESET_ALL}")
    except Exception as exc:
        logging.error("Error: %s", exc)
        print(Fore.RED + "❗️🚨 An error occurred" + Style.RESET_ALL)
        print(Fore.YELLOW + "🟡 Traceback:" + Style.RESET_ALL)
        print(Fore.YELLOW + traceback.format_exc() + Style.RESET_ALL)


if __name__ == "__main__":
    main()
