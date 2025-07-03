import py3Dmol
import requests
import logging
from colorama import Fore, Style
import traceback
yaml = __import__('yaml')
import json
import os
import datetime
import platform
import getpass
import pkg_resources
import hashlib

print(Fore.CYAN + '[binding_visualizer] Starting script...' + Style.RESET_ALL)

# Load general settings
# general_settings = pu.load_general_settings()

# Load configuration settings
# get name of the current module
module_name = os.path.splitext(os.path.basename(__file__))[0]
this_script_folder_path = os.path.dirname(os.path.realpath(__file__))
config_path = os.path.join(this_script_folder_path, module_name + ".yaml")
print(Fore.CYAN + f"[INFO] Loading configuration from {config_path}" + Style.RESET_ALL)
if not os.path.exists(config_path):
    print(Fore.RED + f"Configuration file {config_path} not found." + Style.RESET_ALL)
    exit(1)
with open(config_path, "r") as config_file:
    config = yaml.safe_load(config_file)
print(Fore.GREEN + "[SUCCESS] Configuration loaded successfully." + Style.RESET_ALL)

# ...existing code...

def fetch_pdb_data(pdb_id):
    """
    Fetches PDB data for the given ID from the RCSB PDB database.
    # ...existing code...
    """
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        print(Fore.CYAN + f"[INFO] Fetching PDB data for ID: {pdb_id}" + Style.RESET_ALL)
        # Log the request attempt
        logging.info("Fetching PDB data for ID: %s", pdb_id)
        # Make a request to the PDB URL
        response = requests.get(pdb_url, timeout=10)
        print(Fore.CYAN + f"[INFO] HTTP GET {pdb_url} status: {response.status_code}" + Style.RESET_ALL)
        # Raise an exception if the request was unsuccessful
        response.raise_for_status()
        print(Fore.GREEN + f"[SUCCESS] PDB data fetched for {pdb_id}." + Style.RESET_ALL)
        # Log successful data fetching
        logging.info("PDB data fetched successfully.")
        return response.text
    except requests.exceptions.RequestException as error:
        # Log the error if fetching fails
        logging.error("Error fetching PDB data: %s", error)
        print(Fore.RED + "Error fetching PDB data. Check the log file for details." + Style.RESET_ALL)
        print(Fore.YELLOW + f"[ERROR] {error}" + Style.RESET_ALL)
        raise

# ...existing code...

def file_hash(path):
    """Return SHA256 hash of a file."""
    print(Fore.CYAN + f"[INFO] Calculating SHA256 for {path}" + Style.RESET_ALL)
    h = hashlib.sha256()
    with open(path, 'rb') as f:
        while True:
            chunk = f.read(8192)
            if not chunk:
                break
            h.update(chunk)
    hash_val = h.hexdigest()
    print(Fore.GREEN + f"[SUCCESS] SHA256 for {path}: {hash_val}" + Style.RESET_ALL)
    return hash_val

# ...existing code...

def parse_pdb_header(pdb_data):
    """Extract method, resolution, ligands, chains from PDB header."""
    print(Fore.CYAN + "[INFO] Parsing PDB header for metadata..." + Style.RESET_ALL)
    method = None
    resolution = None
    ligands = set()
    chains = set()
    for line in pdb_data.splitlines():
        if line.startswith('EXPDTA'):
            method = line[10:].strip()
        elif line.startswith('REMARK   2') and 'RESOLUTION.' in line:
            parts = line.split()
            for i, p in enumerate(parts):
                if p == 'RESOLUTION.':
                    try:
                        resolution = parts[i+1] + ' ' + parts[i+2]
                    except Exception:
                        pass
        elif line.startswith('HET   '):
            het_code = line[7:10].strip()
            if het_code and het_code != 'HOH':
                ligands.add(het_code)
        elif line.startswith('COMPND') and 'CHAIN:' in line:
            chain_part = line.split('CHAIN:')[1].split(';')[0]
            for c in chain_part.split(','):
                chains.add(c.strip())
    print(Fore.GREEN + f"[SUCCESS] Parsed header: method={method}, resolution={resolution}, ligands={ligands}, chains={chains}" + Style.RESET_ALL)
    return {
        'method': method,
        'resolution': resolution,
        'ligands': sorted(ligands),
        'chains': sorted(chains)
    }

# ...existing code...
# (The rest of the script is unchanged and continues as in the original file)
