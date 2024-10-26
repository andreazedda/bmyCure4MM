"""
    This module contains the PathsGenerator class, which is responsible for generating the paths for the data files.
    All paths are stored in the settings dictionary, which is then saved to a JSON file.
"""

import processes_utils as pu
import os
import json
import logging

class PathsGenerator:
    """
    This class generates the paths for the data files and saves them to a JSON file.
    """

    def __init__(self):
        """
        Initialize the PathsGenerator class.
        """
        # Load general settings
        self.general_settings = pu.load_general_settings()

        # Initialize the settings dictionary
        self.settings = {}

    def generate_paths(self):
        """
        Generate the paths for the data files.
        """
        # Set the data path
        self.settings["data_path"] = os.path.join(self.general_settings["root_path"], "data")

        # Set the logs path
        self.settings["logs_path"] = os.path.join(self.settings["data_path"], "logs")

        # Set the outputs path
        self.settings["outputs_path"] = os.path.join(self.settings["data_path"], "outputs")

        # Set the pickles path
        self.settings["pickles_path"] = os.path.join(self.settings["outputs_path"], "pickles")

    def save_paths(self):
        """
        Save the paths to a JSON file.
        """
        # Create the data directory if it does not exist
        if not os.path.exists(self.settings["data_path"]):
            os.makedirs(self.settings["data_path"])

        # Create the logs directory if it does not exist
        if not os.path.exists(self.settings["logs_path"]):
            os.makedirs(self.settings["logs_path"])

        # Create the outputs directory if it does not exist
        if not os.path.exists(self.settings["outputs_path"]):
            os.makedirs(self.settings["outputs_path"])

        # Create the pickles directory if it does not exist
        if not os.path.exists(self.settings["pickles_path"]):
            os.makedirs(self.settings["pickles_path"])

        # Save the settings to a JSON file
        with open(os.path.join(self.settings["data_path"], "paths.json"), "w") as settings_file:
            json.dump(self.settings, settings_file, indent=4)
            
        # Log the successful saving of the paths 
        logging.info("Paths saved successfully.")
        
                
def main():
    """
    Main function to execute the paths generation workflow.
    """
    # Initialize the PathsGenerator
    paths_generator = PathsGenerator()

    # Generate the paths
    paths_generator.generate_paths()

    # Save the paths
    paths_generator.save_paths()
    
# Entry point for script execution

if __name__ == "__main__":
    main()