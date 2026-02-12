import json
import os

def load_general_settings():
    """
    Load the general settings for the pipeline.

    Returns:
        dict: The general settings for the pipeline.
    """
    
    # Discover pipelines root folder.
    pipelines_root = os.path.dirname(os.path.abspath(__file__))
    while os.path.basename(pipelines_root) != "pipelines":
        pipelines_root = os.path.dirname(pipelines_root)

    configs_folder = os.path.join(pipelines_root, "configs")
    settings_path = os.path.join(configs_folder, "general_settings.json")
    with open(settings_path, "r") as settings_file:
        general_settings = json.load(settings_file)

    # Default layout (portable across machines).
    defaults = {
        "root_path": pipelines_root,
        "configs_path": os.path.join(pipelines_root, "configs"),
        "data_path": os.path.join(pipelines_root, "data"),
        "logs_path": os.path.join(pipelines_root, "data", "logs"),
        "outputs_path": os.path.join(pipelines_root, "data", "outputs"),
        "pickles_path": os.path.join(pipelines_root, "data", "outputs", "pickles"),
    }

    # Fill missing keys.
    for key, value in defaults.items():
        general_settings.setdefault(key, value)

    # Expand relative paths (relative to pipelines_root).
    for key in defaults.keys():
        val = general_settings.get(key)
        if isinstance(val, str) and val:
            if not os.path.isabs(val):
                general_settings[key] = os.path.normpath(os.path.join(pipelines_root, val))
            elif key in ("root_path", "configs_path", "data_path", "logs_path", "outputs_path", "pickles_path"):
                # If a machine-specific absolute path doesn't exist, fall back to defaults.
                if not os.path.exists(val):
                    general_settings[key] = defaults[key]

    return general_settings