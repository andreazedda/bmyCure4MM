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
from rdkit import Chem
from rdkit.Chem import Descriptors
import processes_utils as pu

# Load general settings
general_settings = pu.load_general_settings()

# Load configuration settings
module_name = os.path.splitext(os.path.basename(__file__))[0]
with open(os.path.join(general_settings['configs_path'], module_name + ".yaml"), "r") as config_file:
    config = yaml.safe_load(config_file)


def fetch_drug_data(drug_id):
    """
    Fetches SMILES data for the given drug ID from PubChem or another database.
    """
    drug_url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{drug_id}/property/CanonicalSMILES/JSON"
    try:
        logging.info("Fetching drug data for ID: %s", drug_id)
        response = requests.get(drug_url, timeout=10)
        response.raise_for_status()
        smiles_data = response.json()['PropertyTable']['Properties'][0]['CanonicalSMILES']
        logging.info("Drug data fetched successfully.")
        return smiles_data
    except requests.exceptions.RequestException as error:
        logging.error("Error fetching drug data: %s", error)
        print(Fore.RED + "Error fetching drug data. Check the log file for details." + Style.RESET_ALL)
        raise


def calculate_drug_parameters(smiles):
    """
    Calculates drug-related parameters using the SMILES string.
    """
    try:
        molecule = Chem.MolFromSmiles(smiles)
        parameters = {
            'Molecular_Weight': Descriptors.MolWt(molecule),
            'LogP': Descriptors.MolLogP(molecule),
            'Num_H_Donors': Descriptors.NumHDonors(molecule),
            'Num_H_Acceptors': Descriptors.NumHAcceptors(molecule)
        }
        logging.info("Parameters calculated: %s", parameters)
        return parameters
    except Exception as error:
        logging.error("Error calculating parameters: %s", error)
        raise


def visualize_drug_structure(smiles, width, height, output_html):
    """
    Visualizes the molecular structure of a drug from SMILES.
    """
    logging.info("Initializing 3Dmol viewer for drug visualization.")
    mol_block = Chem.MolToMolBlock(Chem.MolFromSmiles(smiles))

    viewer = py3Dmol.view(width=width, height=height)
    viewer.addModel(mol_block, 'mol')
    viewer.setStyle({'stick': {}})
    viewer.zoomTo()
    viewer.setBackgroundColor('white')
    viewer.render()

    with open(output_html, 'w') as html_file:
        html_file.write(viewer._make_html())
    
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

        # Visualize the structure
        output_html = os.path.join(general_settings['outputs_path'], f"{config['drug_id']}_structure.html")
        visualize_drug_structure(smiles, config['viewer']['width'], config['viewer']['height'], output_html)

    except Exception as error:
        logging.error("An error occurred in the main function: %s", error)
        print(Fore.RED + "An error occurred. Traceback is shown below:" + Style.RESET_ALL)
        print(Fore.YELLOW + traceback.format_exc() + Style.RESET_ALL)


if __name__ == "__main__":
    main()
