import py3Dmol
import requests
import logging
from colorama import Fore, Style
import traceback
import yaml
import json
import os
import processes_utils as pu

# Load general settings
general_settings = pu.load_general_settings()

# Load configuration settings
# get name of the current module
module_name = os.path.splitext(os.path.basename(__file__))[0]
with open(os.path.join(general_settings['configs_path'], module_name + ".yaml"), "r") as config_file:
    config = yaml.safe_load(config_file)


def fetch_pdb_data(pdb_id):
    """
    Fetches PDB data for the given ID from the RCSB PDB database.

    Args:
        pdb_id (str): The PDB ID of the structure to be fetched.

    Returns:
        str: PDB file data as a string if the request is successful.

    Raises:
        requests.exceptions.RequestException: If an error occurs during the request.
    """
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        # Log the request attempt
        logging.info("Fetching PDB data for ID: %s", pdb_id)

        # Make a request to the PDB URL
        response = requests.get(pdb_url, timeout=10)

        # Raise an exception if the request was unsuccessful
        response.raise_for_status()

        # Log successful data fetching
        logging.info("PDB data fetched successfully.")
        return response.text
    except requests.exceptions.RequestException as error:
        # Log the error if fetching fails
        logging.error("Error fetching PDB data: %s", error)
        print(Fore.RED + "Error fetching PDB data. Check the log file for details." + Style.RESET_ALL)
        raise

def visualize_structure(pdb_data, width, height, chain_style, residue_style, output_html):
    """
    Visualizes the PDB structure using py3Dmol.

    Args:
        pdb_data (str): The PDB data string that contains the molecular structure.
        width (int): Width of the viewer window.
        height (int): Height of the viewer window.
        chain_style (dict): Visualization style for the chain.
        residue_style (dict): Visualization style for the residue.
        output_html (str): Path to the output HTML file to save the visualization.

    Functionality:
        - Creates a viewer window using py3Dmol.
        - Loads the PDB data into the viewer.
        - Sets visualization styles for different parts of the molecule.
        - Saves the viewer as an HTML file.
    """
    # Log that the viewer is being initialized
    logging.info("Initializing 3Dmol viewer.")

    # Create a viewer with specified dimensions
    viewer = py3Dmol.view(width=width, height=height)

    # Load the PDB data into the viewer
    viewer.addModel(pdb_data, 'pdb')

    # Set visualization styles
    viewer.setStyle({'chain': 'A'}, chain_style)
    viewer.setStyle({'resn': 'BOR'}, residue_style)

    # Adjust the zoom to focus on the loaded structure
    viewer.zoomTo()

    # Set hoverable with a callback to display labels for atoms, and clear labels on unhover
    viewer.setHoverable({}, True, "function(atom, viewer) { \
        if(atom) { \
            viewer.addLabel(atom.chain + ' - ' + atom.resn, { \
                position: { x: atom.x, y: atom.y, z: atom.z }, \
                backgroundColor: 'black', \
                fontColor: 'white', \
                fontSize: 10, \
                showBackground: true \
            }); \
        } \
    }", "function(atom, viewer) { viewer.removeAllLabels(); }")

    # Enable user interaction such as rotation and zoom
    viewer.setBackgroundColor('white')
    viewer.zoom(1.2)
    viewer.render()

    # Save the visualization to an HTML file
    with open(output_html, 'w') as html_file:
        html_file.write(viewer._make_html())

    # Log the completion of visualization
    logging.info("Visualization saved to %s", output_html)
    print(Fore.GREEN + f"Visualization saved to {output_html}. Open this file in a browser to view the structure." + Style.RESET_ALL)

def main():
    """
    Main function to execute the PDB data fetching and visualization workflow.

    Workflow:
        - Uses settings from configuration files for PDB ID, viewer dimensions, and visualization styles.
        - Fetches the PDB data using `fetch_pdb_data`.
        - Visualizes the data using `visualize_structure`.
        - Handles any exceptions, logs errors, and prints a traceback for debugging.
    """
    # Initialize logging configuration
    logging.basicConfig(
        filename=os.path.join(general_settings['logs_path'], "binding_visualizer.log"),
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )

    try:
        # Fetch the PDB data using the helper function
        pdb_data = fetch_pdb_data(config['pdb_id'])

        # Visualize the structure using the fetched data
        output_html = os.path.join(general_settings['outputs_path'], f"{config['pdb_id']}_structure_viewer.html")
        visualize_structure(
            pdb_data, 
            config['viewer']['width'], 
            config['viewer']['height'], 
            config['visualization']['chain_style'], 
            config['visualization']['residue_style'], 
            output_html
        )
    except Exception as error:
        # Log if an error occurs in the main function
        logging.error("An error occurred in the main function: %s", error)

        # Print a red error message with the traceback details
        print(Fore.RED + "An error occurred. Traceback is shown below:" + Style.RESET_ALL)
        print(Fore.YELLOW + traceback.format_exc() + Style.RESET_ALL)

# Entry point for script execution
if __name__ == "__main__":
    main()
