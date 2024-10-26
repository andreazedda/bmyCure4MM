import json
import os

def load_general_settings():
    """
    Load the general settings for the pipeline.

    Returns:
        dict: The general settings for the pipeline.
    """
    
    # get the folder of this script
    script_folder = os.path.dirname(os.path.abspath(__file__))
    # continue to go up one folder until the "pipelines" directory is found
    while os.path.basename(script_folder) != "pipelines":
        script_folder = os.path.dirname(script_folder)
    # find configs directory
    configs_folder = os.path.join(script_folder, "configs")
    # load the general settings from the JSON file
    with open(os.path.join(configs_folder, "general_settings.json"), "r") as settings_file:
        general_settings = json.load(settings_file)
    return general_settings