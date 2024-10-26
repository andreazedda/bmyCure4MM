"""
This module generates the settings for the pipeline, logging configurations, and error handling.
"""

import os
import json
import logging
import traceback
from colorama import init, Fore, Style

# Initialize colorama for Windows compatibility with ANSI colors
init(autoreset=True)

# Get the path from where this script is running
ABSOLUTE_PATH = os.path.dirname(os.path.abspath(__file__))

# Look for parent directory until the "pipelines" directory is found
while os.path.basename(ABSOLUTE_PATH) != "pipelines":
    ABSOLUTE_PATH = os.path.dirname(ABSOLUTE_PATH)

# Find configs directory
CONFIGS_PATH = os.path.join(ABSOLUTE_PATH, "configs")
# Check if the directory exists
if not os.path.exists(CONFIGS_PATH):
    os.makedirs(CONFIGS_PATH)

# Define paths for data, logs, outputs, and pickles directories
DATA_PATH = os.path.join(ABSOLUTE_PATH, "data")
LOGS_PATH = os.path.join(DATA_PATH, "logs")
OUTPUTS_PATH = os.path.join(DATA_PATH, "outputs")
PICKLES_PATH = os.path.join(OUTPUTS_PATH, "pickles")

# Set up logging
if not os.path.exists(LOGS_PATH):
    os.makedirs(LOGS_PATH)

logging.basicConfig(
    filename=os.path.join(LOGS_PATH, 'pipeline_log.log'),
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s'
)

def log_colored_message(level, message):
    color_map = {
        "INFO": Fore.GREEN,
        "WARNING": Fore.YELLOW,
        "ERROR": Fore.RED
    }
    color = color_map.get(level, Fore.WHITE)
    print(f"{color}{level}: {message}{Style.RESET_ALL}")

def save_json_settings():
    """
    Save the pipeline settings to a JSON file.
    """
    try:
        settings = {
            "root_path": ABSOLUTE_PATH,
            "data_path": DATA_PATH,
            "logs_path": LOGS_PATH,
            "outputs_path": OUTPUTS_PATH,
            "pickles_path": PICKLES_PATH,
            "configs_path": CONFIGS_PATH
        }

        with open(os.path.join(CONFIGS_PATH, "general_settings.json"), "w") as settings_file:
            json.dump(settings, settings_file, indent=4)
        logging.info("Pipeline settings saved successfully.")
        log_colored_message("INFO", "Pipeline settings saved successfully.")

    except Exception as e:
        error_message = f"An error occurred: {str(e)}"
        logging.error(error_message)
        log_colored_message("ERROR", error_message)
        log_colored_message("ERROR", traceback.format_exc())

# Execute function
if __name__ == "__main__":
    save_json_settings()
