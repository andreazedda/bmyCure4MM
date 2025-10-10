def ensure_local_3dmoljs(local_js_path):
    if not os.path.exists(local_js_path):
        print(Fore.YELLOW + f"[WARNING] 3Dmol-min.js not found at {local_js_path}. Downloading..." + Style.RESET_ALL)
        url = "https://3Dmol.csb.pitt.edu/build/3Dmol-min.js"
        try:
            r = requests.get(url, timeout=30)
            r.raise_for_status()
            with open(local_js_path, "wb") as f:
                f.write(r.content)
            print(Fore.GREEN + "[SUCCESS] 3Dmol-min.js downloaded and saved." + Style.RESET_ALL)
            # Print first 100 bytes for verification
            with open(local_js_path, "rb") as f:
                first_bytes = f.read(100)
            print(Fore.CYAN + "[DEBUG] First 100 bytes of 3Dmol-min.js: " + repr(first_bytes) + Style.RESET_ALL)
        except Exception as e:
            print(Fore.RED + f"[ERROR] Failed to download 3Dmol-min.js: {e}" + Style.RESET_ALL)
            raise RuntimeError("Failed to download 3Dmol-min.js for offline use.") from e
import py3Dmol
import requests
import logging
from colorama import Fore, Style
import traceback
import yaml
import json
import os
import sys
import datetime
import platform
import getpass
import pkg_resources
import hashlib
import subprocess
import shutil
from jinja2 import Template
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from io import BytesIO
import base64

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

def parse_pdb_header(pdb_data):
    """
    Extract comprehensive metadata from PDB structure data.
    
    Args:
        pdb_data (str): Raw PDB file content
        
    Returns:
        dict: Extracted metadata including structural, experimental, and quality metrics
    """
    print(Fore.CYAN + "[INFO] Parsing PDB header for comprehensive metadata..." + Style.RESET_ALL)
    
    metadata = {
        'method': 'N/A',
        'resolution': 'N/A',
        'space_group': 'N/A',
        'unit_cell': 'N/A',
        'chains': [],
        'ligands': [],
        'water_count': 0,
        'total_atoms': 0,
        'protein_chains': [],
        'structure_quality': 'N/A',
        'deposition_date': 'N/A',
        'authors': 'N/A',
        'organism': 'N/A',
        'expression_system': 'N/A'
    }
    
    ligand_codes = set()
    chain_ids = set()
    protein_chain_ids = set()
    atom_count = 0
    water_count = 0
    
    lines = pdb_data.split('\n')
    
    for line in lines:
        # Basic structure information
        if line.startswith('HEADER'):
            metadata['deposition_date'] = line[50:59].strip() if len(line) > 58 else 'N/A'
            
        elif line.startswith('EXPDTA'):
            metadata['method'] = line[10:].strip()
            
        elif line.startswith('REMARK   2 RESOLUTION'):
            res_text = line[23:].strip()
            if 'ANGSTROM' in res_text:
                try:
                    metadata['resolution'] = f"{float(res_text.split()[0]):.2f} Å"
                except:
                    metadata['resolution'] = res_text
                    
        elif line.startswith('CRYST1'):
            # Unit cell parameters
            try:
                a = float(line[6:15].strip())
                b = float(line[15:24].strip())
                c = float(line[24:33].strip())
                alpha = float(line[33:40].strip())
                beta = float(line[40:47].strip())
                gamma = float(line[47:54].strip())
                space_group = line[55:66].strip()
                
                metadata['unit_cell'] = f"a={a:.1f} b={b:.1f} c={c:.1f} α={alpha:.1f}° β={beta:.1f}° γ={gamma:.1f}°"
                metadata['space_group'] = space_group
            except:
                pass
                
        elif line.startswith('AUTHOR'):
            if metadata['authors'] == 'N/A':
                metadata['authors'] = line[10:].strip()
            else:
                metadata['authors'] += ', ' + line[10:].strip()
                
        elif line.startswith('SOURCE'):
            if 'ORGANISM_SCIENTIFIC:' in line:
                metadata['organism'] = line.split('ORGANISM_SCIENTIFIC:')[1].split(';')[0].strip()
            elif 'EXPRESSION_SYSTEM:' in line:
                metadata['expression_system'] = line.split('EXPRESSION_SYSTEM:')[1].split(';')[0].strip()
                
        elif line.startswith('ATOM') or line.startswith('HETATM'):
            atom_count += 1
            chain_id = line[21].strip()
            chain_ids.add(chain_id)
            
            if line.startswith('ATOM'):
                protein_chain_ids.add(chain_id)
            elif line[17:20].strip() == 'HOH':
                water_count += 1
                
        elif line.startswith('HETATM'):
            ligand_code = line[17:20].strip()
            if ligand_code not in ['HOH', 'WAT']:  # Exclude water
                ligand_codes.add(ligand_code)
        
        elif line.startswith('HET   '):
            het_code = line[7:10].strip()
            if het_code and het_code != 'HOH':
                ligand_codes.add(het_code)
    
    metadata['chains'] = sorted(list(chain_ids))
    metadata['protein_chains'] = sorted(list(protein_chain_ids))
    metadata['ligands'] = sorted(list(ligand_codes))
    metadata['total_atoms'] = atom_count
    metadata['water_count'] = water_count
    
    # Calculate basic quality metrics
    if metadata['resolution'] != 'N/A' and 'Å' in metadata['resolution']:
        try:
            res_value = float(metadata['resolution'].split()[0])
            if res_value <= 1.5:
                metadata['structure_quality'] = 'High (≤1.5Å)'
            elif res_value <= 2.0:
                metadata['structure_quality'] = 'Good (1.5-2.0Å)'
            elif res_value <= 2.5:
                metadata['structure_quality'] = 'Moderate (2.0-2.5Å)'
            else:
                metadata['structure_quality'] = 'Low (>2.5Å)'
        except:
            pass
    
    print(Fore.GREEN + f"[SUCCESS] Parsed comprehensive metadata: method={metadata['method']}, resolution={metadata['resolution']}, ligands={metadata['ligands']}, chains={metadata['chains']}" + Style.RESET_ALL)
    
    return metadata

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
    local_js_path = os.path.join(this_script_folder_path, "3Dmol-min.js")
    ensure_local_3dmoljs(local_js_path)
    print(Fore.CYAN + f"[INFO] Initializing 3Dmol viewer with width={width}, height={height}" + Style.RESET_ALL)
    # Log that the viewer is being initialized
    logging.info("Initializing 3Dmol viewer.")

    # Create a viewer with specified dimensions
    viewer = py3Dmol.view(width=width, height=height)
    print(Fore.CYAN + "[INFO] Adding model to viewer..." + Style.RESET_ALL)
    viewer.addModel(pdb_data, 'pdb')

    print(Fore.CYAN + "[INFO] Setting visualization styles..." + Style.RESET_ALL)
    # Set visualization styles
    viewer.setStyle({'chain': 'A'}, chain_style)
    viewer.setStyle({'resn': 'BOR'}, residue_style)

    print(Fore.CYAN + "[INFO] Zooming and rendering viewer..." + Style.RESET_ALL)
    # Adjust the zoom to focus on the loaded structure
    viewer.zoomTo()
    # Set hoverable with a callback to display labels for atoms, and clear labels on unhover
    print(Fore.CYAN + "[INFO] Setting hoverable labels..." + Style.RESET_ALL)
    viewer.setHoverable({}, True,
        'function(atom, viewer) {\n'
        '  if(atom) {\n'
        '    viewer.addLabel(atom.chain + " - " + atom.resn, {\n'
        '      position: { x: atom.x, y: atom.y, z: atom.z },\n'
        '      backgroundColor: "black",\n'
        '      fontColor: "white",\n'
        '      fontSize: 10,\n'
        '      showBackground: true\n'
        '    });\n'
        '  }\n'
        '}\n',
        'function(atom, viewer) { viewer.removeAllLabels(); }'
    )

    # Enable user interaction such as rotation and zoom
    viewer.setBackgroundColor('white')
    viewer.zoom(1.2)
    viewer.render()
    print(Fore.GREEN + "[SUCCESS] 3Dmol viewer rendered." + Style.RESET_ALL)

    # Prepare additional HTML info with more color and information
    # Fetch PDB link
    pdb_id = config['pdb_id']
    pdb_link = f"https://www.rcsb.org/structure/{pdb_id}"
    config_file_info = config_path
    user = getpass.getuser()
    host = platform.node()
    print(Fore.CYAN + f"[INFO] Collecting dependency versions..." + Style.RESET_ALL)
    deps = {
        'py3Dmol': pkg_resources.get_distribution('py3Dmol').version if pkg_resources.working_set.by_key.get('py3dmol') else 'N/A',
        'requests': pkg_resources.get_distribution('requests').version if pkg_resources.working_set.by_key.get('requests') else 'N/A',
        'pyyaml': pkg_resources.get_distribution('pyyaml').version if pkg_resources.working_set.by_key.get('pyyaml') else 'N/A',
    }
    try:
        import subprocess
        print(Fore.CYAN + "[INFO] Getting git commit hash..." + Style.RESET_ALL)
        git_hash = subprocess.check_output(['git', 'rev-parse', '--short', 'HEAD'], cwd=this_script_folder_path).decode().strip()
    except Exception:
        git_hash = 'N/A'
    pdb_title = ''
    for line in pdb_data.splitlines():
        if line.startswith('TITLE '):
            pdb_title += line[10:].strip() + ' '
    pdb_title = pdb_title.strip()
    print(Fore.CYAN + "[INFO] Parsing PDB header for method, resolution, ligands, chains..." + Style.RESET_ALL)
    # --- New: Parse header ---
    pdb_meta = parse_pdb_header(pdb_data)
    print(Fore.CYAN + "[INFO] Calculating script and config SHA256 hashes..." + Style.RESET_ALL)
    # --- New: Hashes ---
    script_hash = file_hash(__file__)
    config_hash = file_hash(config_path)
    print(Fore.CYAN + "[INFO] Reading inline YAML config..." + Style.RESET_ALL)
    # --- New: Inline YAML ---
    with open(config_path, 'r') as f:
        config_yaml = f.read()
    # --- New: Download timestamp ---
    pdb_download_time = datetime.datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S UTC')
    print(Fore.CYAN + "[INFO] Checking for Bio.PDB (biopython) availability..." + Style.RESET_ALL)
    try:
        from Bio.PDB import PDBParser, NeighborSearch
        biopython_available = True
        print(Fore.GREEN + "[SUCCESS] Bio.PDB available. Binding site mapping enabled." + Style.RESET_ALL)
    except ImportError:
        biopython_available = False
        print(Fore.YELLOW + "[WARNING] Bio.PDB not available. Binding site mapping disabled." + Style.RESET_ALL)
    def get_binding_site_residues(pdb_path, ligand_code, cutoff=5.0):
        if not biopython_available:
            print(Fore.YELLOW + "[binding_visualizer] Skipping binding site mapping (Bio.PDB not available)." + Style.RESET_ALL)
            return []
        print(Fore.CYAN + f"[INFO] Mapping binding site for ligand {ligand_code} with cutoff {cutoff} Å..." + Style.RESET_ALL)
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('PDB', pdb_path)
        atoms = [atom for atom in structure.get_atoms()]
        ns = NeighborSearch(atoms)
        ligand_atoms = [atom for atom in structure.get_atoms() if atom.parent.resname == ligand_code]
        binding_residues = set()
        for atom in ligand_atoms:
            for neighbor in ns.search(atom.coord, cutoff):
                res = neighbor.parent
                if res.get_id()[0] == " " and res.resname != ligand_code:
                    binding_residues.add((res.parent.id, res.id[1], res.resname))
        print(Fore.GREEN + f"[SUCCESS] Found {len(binding_residues)} binding site residues." + Style.RESET_ALL)
        return sorted(binding_residues)
    def fetch_rcsb_validation(pdb_id):
        print(Fore.CYAN + f"[INFO] Fetching RCSB validation for {pdb_id}..." + Style.RESET_ALL)
        try:
            url = f"https://validate-rcsb-1.wwpdb.org/api/validation/entry/{pdb_id.lower()}"
            resp = requests.get(url, timeout=10)
            if resp.status_code == 200:
                print(Fore.GREEN + f"[SUCCESS] Validation data fetched." + Style.RESET_ALL)
                return resp.json()
        except Exception:
            print(Fore.YELLOW + f"[WARNING] Validation fetch failed." + Style.RESET_ALL)
            pass
        return None
    def fetch_rcsb_citation(pdb_id):
        print(Fore.CYAN + f"[INFO] Fetching RCSB citation for {pdb_id}..." + Style.RESET_ALL)
        try:
            url = f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.lower()}"
            resp = requests.get(url, timeout=10)
            if resp.status_code == 200:
                print(Fore.GREEN + f"[SUCCESS] Citation data fetched." + Style.RESET_ALL)
                return resp.json()
        except Exception:
            print(Fore.YELLOW + f"[WARNING] Citation fetch failed." + Style.RESET_ALL)
            pass
        return None
    ligand_code = config.get('ligand', 'BOR')
    binding_site_cfg = config.get('binding_site_detection', {'enabled': False})
    mutation_list = config.get('mutations', [])
    therapies = config.get('therapies', [])
    pathways = config.get('pathways', [])
    literature = config.get('literature', [])
    print(Fore.CYAN + f"[INFO] Saving PDB to temp file for biopython: {pdb_id}.pdb" + Style.RESET_ALL)
    # --- Save PDB to temp file for biopython ---
    pdb_temp_path = os.path.join(this_script_folder_path, f"{pdb_id}.pdb")
    with open(pdb_temp_path, 'w') as f:
        f.write(pdb_data)
    binding_residues = []
    if binding_site_cfg.get('enabled', False) and biopython_available:
        cutoff = float(binding_site_cfg.get('cutoff_angstrom', 5.0))
        print(Fore.CYAN + f"[INFO] Running binding site detection (cutoff={cutoff})..." + Style.RESET_ALL)
        # --- Binding site detection ---
        binding_residues = get_binding_site_residues(pdb_temp_path, ligand_code, cutoff)
    else:
        print(Fore.YELLOW + "[WARNING] Binding site detection not enabled or Bio.PDB unavailable." + Style.RESET_ALL)
    print(Fore.CYAN + f"[INFO] Fetching RCSB validation/citation..." + Style.RESET_ALL)
    # --- RCSB validation/citation ---
    validation = fetch_rcsb_validation(pdb_id)
    citation = fetch_rcsb_citation(pdb_id)
    mutation_table_html = ""
    if mutation_list:
        print(Fore.CYAN + f"[INFO] Overlaying {len(mutation_list)} mutations..." + Style.RESET_ALL)
        mutation_table_html = "<table style='border-collapse:collapse;width:100%;margin-bottom:8px;'><tr><th>Chain</th><th>ResNum</th><th>Mutation</th><th>Effect</th></tr>"
        for m in mutation_list:
            mutation_table_html += f"<tr><td>{m.get('chain')}</td><td>{m.get('resnum')}</td><td>{m.get('mutation')}</td><td>{m.get('effect','')}</td></tr>"
        mutation_table_html += "</table>"
    binding_table_html = ""
    if binding_residues:
        print(Fore.CYAN + f"[INFO] Generating binding site table..." + Style.RESET_ALL)
        binding_table_html = "<table style='border-collapse:collapse;width:100%;margin-bottom:8px;'><tr><th>Chain</th><th>ResNum</th><th>ResName</th></tr>"
        for chain, resnum, resname in binding_residues:
            binding_table_html += f"<tr><td>{chain}</td><td>{resnum}</td><td>{resname}</td></tr>"
        binding_table_html += "</table>"
    therapies_html = ""
    if therapies:
        print(Fore.CYAN + f"[INFO] Adding {len(therapies)} therapies to table..." + Style.RESET_ALL)
        therapies_html = "<table style='border-collapse:collapse;width:100%;margin-bottom:8px;'><tr><th>Name</th><th>Ligand</th><th>Phase</th><th>MM Relevance</th><th>Resistance Mutations</th></tr>"
        for t in therapies:
            therapies_html += f"<tr><td>{t.get('name')}</td><td>{t.get('pdb_ligand')}</td><td>{t.get('clinical_phase','')}</td><td>{t.get('mm_relevance','')}</td><td>{', '.join(t.get('resistance_mutations',[]))}</td></tr>"
        therapies_html += "</table>"
    pathways_html = ""
    if pathways:
        print(Fore.CYAN + f"[INFO] Adding {len(pathways)} pathways to table..." + Style.RESET_ALL)
        pathways_html = "<ul>"
        for p in pathways:
            kegg = f"<a href='https://www.kegg.jp/dbget-bin/www_bget?{p.get('kegg_id')}' target='_blank'>{p.get('kegg_id')}</a>" if p.get('kegg_id') else ''
            reactome = f"<a href='https://reactome.org/content/detail/{p.get('reactome_id')}' target='_blank'>{p.get('reactome_id')}</a>" if p.get('reactome_id') else ''
            pathways_html += f"<li>{p.get('name','')} {kegg} {reactome}</li>"
        pathways_html += "</ul>"
    literature_html = ""
    if literature:
        print(Fore.CYAN + f"[INFO] Adding {len(literature)} literature references..." + Style.RESET_ALL)
        literature_html = "<ul>"
        for l in literature:
            doi = f"<a href='https://doi.org/{l.get('doi')}' target='_blank'>{l.get('doi')}</a>" if l.get('doi') else ''
            pmid = f"<a href='https://pubmed.ncbi.nlm.nih.gov/{l.get('pmid')}' target='_blank'>{l.get('pmid')}</a>" if l.get('pmid') else ''
            literature_html += f"<li>{l.get('title','')} {doi} {pmid}</li>"
        literature_html += "</ul>"
    validation_html = ""
    if validation:
        outliers = validation.get('geometry_quality', {}).get('ramachandran_outliers', {}).get('percent', None)
        validation_html = f"<div style='color:#b03a2e;'><b>Validation:</b> Ramachandran outliers: {outliers}%</div>" if outliers is not None else ''
        print(Fore.CYAN + f"[INFO] Validation summary: Ramachandran outliers: {outliers}%" + Style.RESET_ALL)
    citation_html = ""
    if citation:
        try:
            pub = citation['rcsb_primary_citation']
            title = pub.get('title','')
            doi = pub.get('doi','')
            citation_html = f"<div><b>Citation:</b> <a href='https://doi.org/{doi}' target='_blank'>{title}</a></div>"
            print(Fore.CYAN + f"[INFO] Citation: {title} (DOI: {doi})" + Style.RESET_ALL)
        except Exception:
            print(Fore.YELLOW + f"[WARNING] Citation parsing failed." + Style.RESET_ALL)
            pass
    print(Fore.CYAN + f"[INFO] Preparing HTML info panel..." + Style.RESET_ALL)
    info_html = f"""
    <div style='font-family:Arial,sans-serif;font-size:14px;color:#333;margin-bottom:10px;'>
        <h2 style='color:#2a5298;margin-top:0;'>PDB Structure Viewer</h2>
        <div style='margin-bottom:8px;'><span style='color:#1e8449;font-weight:bold;'>PDB ID:</span> <span style='color:#154360;'><a href='https://www.rcsb.org/structure/{pdb_id}' target='_blank'>{pdb_id}</a></span></div>
        <div style='margin-bottom:8px;'><span style='color:#76448a;font-weight:bold;'>Title:</span> <span style='color:#154360;'>{pdb_title}</span></div>
        <div style='margin-bottom:8px;'><span style='color:#76448a;font-weight:bold;'>Experimental Method:</span> <span style='color:#154360;'>{pdb_meta['method']}</span></div>
        <div style='margin-bottom:8px;'><span style='color:#76448a;font-weight:bold;'>Resolution:</span> <span style='color:#154360;'>{pdb_meta['resolution']}</span></div>
        <div style='margin-bottom:8px;'><span style='color:#76448a;font-weight:bold;'>Ligands:</span> <span style='color:#154360;'>{', '.join(pdb_meta['ligands'])}</span></div>
        <div style='margin-bottom:8px;'><span style='color:#76448a;font-weight:bold;'>Chains:</span> <span style='color:#154360;'>{', '.join(pdb_meta['chains'])}</span></div>
        <div style='margin-bottom:8px;'><span style='color:#b9770e;font-weight:bold;'>Viewer size:</span> <span style='color:#154360;'>{width} x {height}</span></div>
        <div style='margin-bottom:8px;'><span style='color:#2874a6;font-weight:bold;'>Chain style:</span> <span style='color:#154360;'>{chain_style}</span></div>
        <div style='margin-bottom:8px;'><span style='color:#a93226;font-weight:bold;'>Residue style:</span> <span style='color:#154360;'>{residue_style}</span></div>
        <div style='margin-bottom:8px;'><span style='color:#884ea0;font-weight:bold;'>Generated:</span> <span style='color:#154360;'>{datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</span></div>
        <div style='margin-bottom:8px;'><span style='color:#117864;font-weight:bold;'>Script:</span> <span style='color:#154360;'>{os.path.basename(__file__)}</span></div>
        <div style='margin-bottom:8px;'><span style='color:#2874a6;font-weight:bold;'>Config file:</span> <span style='color:#154360;'>{os.path.basename(config_path)}</span></div>
        <div style='margin-bottom:8px;'><span style='color:#b03a2e;font-weight:bold;'>Python version:</span> <span style='color:#154360;'>{os.sys.version.split()[0]}</span></div>
        <div style='margin-bottom:8px;'><span style='color:#229954;font-weight:bold;'>Platform:</span> <span style='color:#154360;'>{platform.system()} {platform.release()}</span></div>
        <div style='margin-bottom:8px;'><span style='color:#d35400;font-weight:bold;'>User/Host:</span> <span style='color:#154360;'>{user}@{host}</span></div>
        <div style='margin-bottom:8px;'><span style='color:#2874a6;font-weight:bold;'>PDB Downloaded:</span> <span style='color:#154360;'>{pdb_download_time}</span></div>
        <div style='margin-bottom:8px;'><span style='color:#2874a6;font-weight:bold;'>Script SHA256:</span> <span style='color:#154360;'>{script_hash}</span></div>
        <div style='margin-bottom:8px;'><span style='color:#2874a6;font-weight:bold;'>Config SHA256:</span> <span style='color:#154360;'>{config_hash}</span></div>
        <div style='margin-bottom:8px;'><span style='color:#2874a6;font-weight:bold;'>Instructions:</span> <span style='color:#154360;'>Drag to rotate, scroll to zoom, double-click to center. Hover over atoms for details.</span></div>
        <div style='margin-bottom:8px;'><span style='color:#2874a6;font-weight:bold;'>Dependencies:</span> <span style='color:#154360;'>{', '.join([f'{k} {v}' for k,v in deps.items()])}</span></div>
        <div style='margin-bottom:8px;'><span style='color:#2874a6;font-weight:bold;'>Git commit:</span> <span style='color:#154360;'>{git_hash}</span></div>
        <div style='margin-bottom:8px;'><span style='color:#2874a6;font-weight:bold;'>Downloads:</span> <span style='color:#154360;'><a href='https://files.rcsb.org/download/{pdb_id}.pdb'>PDB file</a> | <a href='binding_visualizer.yaml'>YAML config</a></span></div>
        <details style='margin-bottom:8px;'><summary style='font-weight:bold;color:#2874a6;'>Show YAML Config</summary><pre style='background:#fff;border:1px solid #ccc;border-radius:4px;padding:8px;'>{config_yaml}</pre></details>
        {validation_html}
        {citation_html}
        {binding_table_html}
        {mutation_table_html}
        {therapies_html}
        {pathways_html}
        {literature_html}
        <div style='font-size:12px;color:#888;'>Visualization generated by <b>binding_visualizer.py</b> using <b>py3Dmol</b>.<br>Contact: <a href='mailto:your.email@example.com'>your.email@example.com</a><br>Copyright &copy; {datetime.datetime.now().year}</div>
    </div>
    """
    # --- Interactive Controls Injection ---
    # Prepare JS arrays for chains, ligands, mutations
    js_chains = pdb_meta['chains']
    js_ligands = pdb_meta['ligands']
    js_mutations = [
        {'chain': m.get('chain'), 'resnum': m.get('resnum'), 'mutation': m.get('mutation')}
        for m in mutation_list
    ]
    mutation_js_array = json.dumps([
        {"chain": m["chain"], "resi": m["resnum"]}
        for m in mutation_list
    ])
    # --- New: Interactive Controls Panel ---
    # Build button HTML outside the f-string to avoid curly brace issues
    def safe_viewer_call(js_code):
        return f"if(window.viewer && typeof window.viewer.setStyle==='function'){{{{{js_code}}}}}else{{{{alert('3D viewer not ready yet. Please reload the page.');}}}}"
    chain_buttons = ''
    for c in js_chains:
        js = f"window.viewer.setStyle({{chain: '{c}'}}, {{stick:{{colorscheme:'redCarbon'}}}});window.viewer.render();"
        chain_buttons += f'<button onclick="{safe_viewer_call(js)}" disabled>{c}</button>'
    ligand_buttons = ''
    for l in js_ligands:
        js = f"window.viewer.setStyle({{resn: '{l}'}}, {{stick:{{colorscheme:'orangeCarbon'}}}});window.viewer.render();"
        ligand_buttons += f'<button onclick="{safe_viewer_call(js)}" disabled>{l}</button>'
    reset_all_js = "window.viewer.setStyle({{}}, {{stick:{{}}}});window.viewer.render();"
    cartoon_js = "window.viewer.setStyle({{}}, {{cartoon:{{}}}});window.viewer.render();"
    stick_js = "window.viewer.setStyle({{}}, {{stick:{{}}}});window.viewer.render();"
    sphere_js = "window.viewer.setStyle({{}}, {{sphere:{{}}}});window.viewer.render();"
    line_js = "window.viewer.setStyle({{}}, {{line:{{}}}});window.viewer.render();"
    spectrum_js = "window.viewer.setStyle({{}}, {{cartoon: {{color: 'spectrum'}}}});window.viewer.render();"
    by_chain_js = "window.viewer.setStyle({{}}, {{cartoon: {{color: 'chain'}}}});window.viewer.render();"
    by_element_js = "window.viewer.setStyle({{}}, {{cartoon: {{color: 'element'}}}});window.viewer.render();"
    by_residue_js = "window.viewer.setStyle({{}}, {{cartoon: {{color: 'residue'}}}});window.viewer.render();"
    controls_panel = f'''
    <div id="bv-controls" style="margin:16px 0 24px 260px; padding:12px 16px; background:#f8f9fa; border-radius:8px; border:1px solid #e1e4e8;">
      <b>Interactive Controls:</b>
      <div style="margin-top:8px;">
        <b>Chains:</b>
        {chain_buttons}
        <button onclick="{safe_viewer_call(reset_all_js)}" disabled>Reset All</button>
      </div>
      <div style="margin-top:8px;">
        <b>Ligands:</b>
        {ligand_buttons}
      </div>
      <div style="margin-top:8px;">
        <button onclick="highlightMutations()" disabled>Highlight Mutations</button>
      </div>
      <div style="margin-top:8px;">
        <b>Style:</b>
        <button onclick="{safe_viewer_call(cartoon_js)}" disabled>Cartoon</button>
        <button onclick="{safe_viewer_call(stick_js)}" disabled>Stick</button>
        <button onclick="{safe_viewer_call(sphere_js)}" disabled>Sphere</button>
        <button onclick="{safe_viewer_call(line_js)}" disabled>Line</button>
      </div>
      <div style="margin-top:8px;">
        <b>Color:</b>
        <button onclick="{safe_viewer_call(spectrum_js)}" disabled>Spectrum</button>
        <button onclick="{safe_viewer_call(by_chain_js)}" disabled>By Chain</button>
        <button onclick="{safe_viewer_call(by_element_js)}" disabled>By Element</button>
        <button onclick="{safe_viewer_call(by_residue_js)}" disabled>By Residue</button>
      </div>
    </div>
    <script>
      function enableControls() {{{{
        var btns = document.querySelectorAll('#bv-controls button');
        btns.forEach(function(b) {{{{ b.disabled = false; }}}});
      }}}}
      var mutations = {mutation_js_array};
      function highlightMutations() {{{{
        if(window.viewer && typeof window.viewer.setStyle==='function'){{{{
          window.viewer.setStyle({{}}, {{}}); // clear
          mutations.forEach(function(m) {{{{
            window.viewer.setStyle({{chain: m.chain, resi: m.resi}}, {{stick: {{colorscheme: 'magentaCarbon'}}}});
          }}}});
          window.viewer.render();
        }}}} else {{{{
          alert('3D viewer not ready yet. Please reload the page.');
        }}}}
      }}}}
    </script>
    '''
    # --- End Interactive Controls Panel ---

    # Save the visualization to an HTML file, prepending info
    import re
    def inject_js_head(html, js_path="3Dmol-min.js"):
        js_tag = f'<script src="{js_path}"></script>'
        if "<head>" in html:
            # Only add if not already present
            if js_tag not in html:
                html = html.replace("<head>", "<head>\n" + js_tag)
        else:
            # No <head> found, inject at top
            html = js_tag + "\n" + html
        return html
    with open(output_html, 'w') as html_file:
        html = viewer._make_html()
        # --- Extract the 3Dmol viewer <div> block and script ---
        viewer_div_block = None
        viewer_div_id = None
        viewer_div_style = None
        viewer_script_block = None
        # Find the <div id="3dmolviewer_...">...</div> block
        div_pattern = re.compile(r'(<div\s+id="(3dmolviewer_[^"]+)"([^>]*)>.*?</div>)', re.DOTALL)
        div_match = div_pattern.search(html)
        if div_match:
            viewer_div_block = div_match.group(1)
            viewer_div_id = div_match.group(2)
            viewer_div_style = div_match.group(3)
            print(Fore.CYAN + f"[DIAGNOSTIC] Found viewer block id: {viewer_div_id}" + Style.RESET_ALL)
            print(Fore.CYAN + f"[DIAGNOSTIC] Viewer block HTML (first 100 chars): {viewer_div_block[:100].replace(chr(10),' ')}..." + Style.RESET_ALL)
            # --- Extract the viewer creation <script> block ---
            # The viewer variable is usually viewer_xxxxx, derived from the div id
            viewer_var = viewer_div_id.replace("3dmolviewer_", "viewer_")
            # Try to find the script block that creates the viewer
            script_pattern = re.compile(r'(<script>.*?var\s+' + re.escape(viewer_var) + r'\s*=.*?createViewer.*?</script>)', re.DOTALL)
            script_match = script_pattern.search(html)
            if script_match:
                viewer_script_block = script_match.group(1)
            else:
                # Fallback: try to find any script block that contains createViewer and the viewer_var
                script_pattern2 = re.compile(r'(<script>.*?' + re.escape(viewer_var) + r'.*?createViewer.*?</script>)', re.DOTALL)
                script_match2 = script_pattern2.search(html)
                if script_match2:
                    viewer_script_block = script_match2.group(1)
                else:
                    print(Fore.YELLOW + f"[WARNING] Could not find viewer creation <script> block for {viewer_var}. The viewer may not initialize correctly." + Style.RESET_ALL)
        else:
            # Diagnostic: print all div ids found
            all_ids = re.findall(r'<div\s+id="([^"]+)"', html)
            print(Fore.RED + "[DIAGNOSTIC] Viewer block not found! HTML ids found: " + str(all_ids) + Style.RESET_ALL)
        # --- Remove the viewer block and script from html ---
        html_wo_viewer = html
        if viewer_div_block:
            html_wo_viewer = html_wo_viewer.replace(viewer_div_block, '')
        if viewer_script_block:
            html_wo_viewer = html_wo_viewer.replace(viewer_script_block, '')
        # --- Ensure viewer div has margin-left:260px;min-height:480px; ---
        if viewer_div_block:
            # Ensure style attribute
            if "style=" in viewer_div_block:
                # Patch style string
                style_match = re.search(r'style="([^"]*)"', viewer_div_block)
                if style_match:
                    style_val = style_match.group(1)
                    # Add margin-left:260px;min-height:480px; if not present
                    if 'margin-left:' not in style_val:
                        style_val += ';margin-left:260px;'
                    if 'min-height:' not in style_val:
                        style_val += 'min-height:480px;'
                    # Patch style in div
                    viewer_div_block = re.sub(r'style="[^"]*"', f'style="{style_val}"', viewer_div_block)
            else:
                # Insert style attribute
                viewer_div_block = viewer_div_block.replace('>', ' style="margin-left:260px;min-height:480px;>', 1)
        # Ensure info/controls/panels have margin-left:260px; or are wrapped
        def ensure_margin_left(html_snip):
            # If already has margin-left:260px, return as is
            if "margin-left:260px" in html_snip:
                return html_snip
            # Add margin-left:260px to top div or wrap in div
            m = re.match(r'^(\s*<div[^>]*style=")([^"]*)"', html_snip)
            if m:
                style_str = m.group(2)
                if "margin-left:" not in style_str:
                    new_style = style_str + ";margin-left:260px;"
                    html_snip = html_snip.replace(m.group(0), m.group(1) + new_style + '"', 1)
                return html_snip
            else:
                # Wrap in a div
                return f"<div style='margin-left:260px;'>{html_snip}</div>"
        info_html_margin = ensure_margin_left(info_html)
        controls_panel_margin = ensure_margin_left(controls_panel)

        # --- Compose the HTML ---
        # Find <body> tag in html_wo_viewer
        if '<body>' in html_wo_viewer:
            body_split = html_wo_viewer.split('<body>', 1)
            before_body = body_split[0] + '<body>'
            after_body = body_split[1]
            content = sidebar_html
            if viewer_div_block:
                content += viewer_div_block
                if viewer_script_block:
                    # --- Inject window.viewer assignment script after viewer_script_block ---
                    # Find the viewer variable name
                    viewer_var_match = re.search(r'var (viewer_\\d+) =', viewer_script_block)
                    viewer_var = viewer_var_match.group(1) if viewer_var_match else None
                    assign_script = ''
                    if viewer_var:
                        assign_script = f"""
<script>
(function() {{
  // Assign the py3Dmol viewer variable to window.viewer for global access
  if (typeof window.{viewer_var} === 'object' && typeof window.{viewer_var}.setStyle === 'function') {{
    window.viewer = window.{viewer_var};
    console.log('[binding_visualizer] window.viewer assigned from {viewer_var}', window.viewer);
    enableControls();
  }} else {{
    // fallback for $3Dmol.viewers
    if (window.$3Dmol && $3Dmol.viewers) {{
      for (var id in $3Dmol.viewers) {{
        if ($3Dmol.viewers[id]) {{
          window.viewer = $3Dmol.viewers[id];
          console.log('[binding_visualizer] window.viewer assigned from $3Dmol.viewers', id, window.viewer);
          enableControls();
          break;
        }}
      }}
    }}
    if (!window.viewer) {{
      console.error('[binding_visualizer] FATAL: viewer not found for global assignment!');
    }}
  }}
}})();
</script>
"""
                    content += viewer_script_block + assign_script
                else:
                    content += viewer_script_block if viewer_script_block else ''
            content += info_html_margin + controls_panel_margin
            html_final = before_body + content + after_body
        else:
            # Compose content in the same way as the if branch
            content = sidebar_html
            if viewer_div_block:
                content += viewer_div_block
                if viewer_script_block:
                    # --- Inject window.viewer assignment script after viewer_script_block ---
                    viewer_var_match = re.search(r'var (viewer_\\d+) =', viewer_script_block)
                    viewer_var = viewer_var_match.group(1) if viewer_var_match else None
                    assign_script = ''
                    if viewer_var:
                        assign_script = f"""
<script>
(function() {{
  if (typeof window.{viewer_var} === 'object' && typeof window.{viewer_var}.setStyle === 'function') {{
    window.viewer = window.{viewer_var};
    console.log('[binding_visualizer] window.viewer assigned from {viewer_var}', window.viewer);
    enableControls();
  }} else {{
    if (window.$3Dmol && $3Dmol.viewers) {{
      for (var id in $3Dmol.viewers) {{
        if ($3Dmol.viewers[id]) {{
          window.viewer = $3Dmol.viewers[id];
          console.log('[binding_visualizer] window.viewer assigned from $3Dmol.viewers', id, window.viewer);
          enableControls();
          break;
        }}
      }}
    }}
    if (!window.viewer) {{
      console.error('[binding_visualizer] FATAL: viewer not found for global assignment!');
    }}
  }}
}})();
</script>
"""
                    content += viewer_script_block + assign_script
                else:
                    content += viewer_script_block if viewer_script_block else ''
            content += info_html_margin + controls_panel_margin
            html_final = content

        # --- Final HTML tweaks ---
        # Ensure <html> tag has lang attribute
        if '<html' in html_final and 'lang=' not in html_final:
            html_final = html_final.replace('<html', '<html lang="en"')
        # Ensure <meta> charset is set
        if '<meta' in html_final and 'charset=' not in html_final:
            html_final = html_final.replace('<meta', '<meta charset="UTF-8"')
        # Ensure <title> is set
        if '<title>' not in html_final:
            html_final = html_final.replace('<head>', '<head><title>PDB Structure Viewer</title>')
        # Ensure body has margin
        if '<body' in html_final and 'style=' not in html_final:
            html_final = html_final.replace('<body', '<body style="margin:0;"')

        # --- Debug: Save intermediate HTML for inspection ---
        debug_html_path = os.path.join(this_script_folder_path, f"debug_{pdb_id}.html")
        with open(debug_html_path, 'w') as debug_file:
            debug_file.write(html_final)
        print(Fore.YELLOW + f"[DEBUG] Intermediate HTML saved to {debug_html_path}. Inspect for troubleshooting." + Style.RESET_ALL)

        # Write final HTML to file
        html_file.write(html_final)

    # --- Always append robust assignment script at the end of the HTML file ---
    robust_assign_script = '''
<script>
(function() {
  function enableControls() {
    var btns = document.querySelectorAll('#bv-controls button');
    btns.forEach(function(b) { b.disabled = false; });
  }
  for (var k in window) {
    if (/^viewer_\d+$/.test(k) && typeof window[k] === "object" && typeof window[k].setStyle === "function") {
      window.viewer = window[k];
      console.log("[binding_visualizer] window.viewer assigned from", k, window.viewer);
      enableControls();
      break;
    }
  }
  if (!window.viewer && window.$3Dmol && $3Dmol.viewers) {
    for (var id in $3Dmol.viewers) {
      if ($3Dmol.viewers[id]) {
        window.viewer = $3Dmol.viewers[id];
        console.log("[binding_visualizer] window.viewer assigned from $3Dmol.viewers", id, window.viewer);
        enableControls();
        break;
      }
    }
  }
  if (!window.viewer) {
    console.error("[binding_visualizer] FATAL: viewer not found for global assignment!");
  }
})();
</script>
'''
    # Append to the HTML file
    with open(output_html, 'a') as f:
        f.write(robust_assign_script)

    print(Fore.GREEN + f"[SUCCESS] Visualization saved to {output_html}" + Style.RESET_ALL)
    print(Fore.CYAN + "[INFO] Opening visualization in default web browser..." + Style.RESET_ALL)
    # Open the saved HTML file in the default web browser
    import webbrowser
    webbrowser.open('file://' + os.path.realpath(output_html))

# --- Sidebar HTML (simple placeholder) ---
sidebar_html = f"""
<div id='bv-sidebar' style='position:fixed;top:0;left:0;width:240px;height:100%;background:#2a5298;color:#fff;padding:24px 12px 12px 12px;z-index:1000;overflow-y:auto;'>
  <h2 style='font-size:20px;margin-top:0;margin-bottom:16px;'>MM Structure<br>Visualizer</h2>
  <div style='font-size:14px;line-height:1.5;'>
    <b>Controls:</b><br>
    Use the interactive panel to explore chains, ligands, mutations, and styles.<br><br>
    <b>Tips:</b><br>
    - Drag to rotate<br>
    - Scroll to zoom<br>
    - Double-click to center<br>
    - Hover for atom info<br>
  </div>
  <div style='position:absolute;bottom:16px;left:12px;font-size:11px;color:#cfd8dc;'>
    &copy; {datetime.datetime.now().year} MM Lab
  </div>
</div>
"""
# --- End Sidebar HTML ---

# --- Debug: Save sidebar HTML to file ---
debug_sidebar_path = os.path.join(this_script_folder_path, "debug_sidebar.html")
with open(debug_sidebar_path, 'w') as debug_file:
    debug_file.write(sidebar_html)
print(Fore.YELLOW + f"[DEBUG] Sidebar HTML saved to {debug_sidebar_path}. Inspect for troubleshooting." + Style.RESET_ALL)

def main():
    print(Fore.CYAN + '[DEBUG] Entered main()' + Style.RESET_ALL)
    # Initialize logging configuration
    logging.basicConfig(
        filename=os.path.join(this_script_folder_path, "binding_visualizer.log"),
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s"
    )
    try:
        print(Fore.CYAN + f"[DEBUG] About to fetch PDB data for {config['pdb_id']}..." + Style.RESET_ALL)
        pdb_data = fetch_pdb_data(config['pdb_id'])
        
        # Parse PDB metadata for report generation
        print(Fore.CYAN + f"[DEBUG] Parsing PDB metadata for {config['pdb_id']}..." + Style.RESET_ALL)
        pdb_metadata = parse_pdb_header(pdb_data)
        
        print(Fore.CYAN + f"[DEBUG] About to visualize structure for {config['pdb_id']}..." + Style.RESET_ALL)
        output_html = os.path.join(this_script_folder_path, f"{config['pdb_id']}_structure_viewer.html")
        visualize_structure(
            pdb_data, 
            config['viewer']['width'], 
            config['viewer']['height'], 
            config['visualization']['chain_style'], 
            config['visualization']['residue_style'], 
            output_html
        )
        
        # Generate PDF report if enabled in config
        generate_pdf = config.get('generate_pdf', True)  # Default to True
        if generate_pdf:
            print(Fore.CYAN + f"[DEBUG] Generating PDF report for {config['pdb_id']}..." + Style.RESET_ALL)
            try:
                pdf_path = generate_pdf_report(config, pdb_data, pdb_metadata, this_script_folder_path)
                if pdf_path:
                    print(Fore.GREEN + f"[SUCCESS] PDF report generated: {pdf_path}" + Style.RESET_ALL)
                else:
                    print(Fore.YELLOW + "[WARNING] PDF generation failed, but LaTeX file is available for manual compilation" + Style.RESET_ALL)
            except Exception as pdf_error:
                print(Fore.YELLOW + f"[WARNING] PDF generation failed: {pdf_error}" + Style.RESET_ALL)
                logging.warning("PDF generation failed: %s", pdf_error)
        
        print(Fore.GREEN + f"[DEBUG] Workflow completed for {config['pdb_id']}!" + Style.RESET_ALL)
    except Exception as error:
        logging.error("An error occurred in the main function: %s", error)
        print(Fore.RED + "[DEBUG] An error occurred. Traceback is shown below:" + Style.RESET_ALL)
        print(Fore.YELLOW + traceback.format_exc() + Style.RESET_ALL)

if __name__ == "__main__":
    main()

def generate_structure_image(pdb_data, width=800, height=600, output_path=None):
    """
    Generate a static image of the molecular structure for PDF inclusion.
    
    Args:
        pdb_data (str): PDB structure data
        width (int): Image width in pixels
        height (int): Image height in pixels
        output_path (str): Path to save the image (optional)
    
    Returns:
        str: Base64 encoded image data or file path
    """
    print(Fore.CYAN + "[INFO] Generating static structure image for PDF..." + Style.RESET_ALL)
    
    try:
        # Create a py3Dmol viewer for image generation
        viewer = py3Dmol.view(width=width, height=height)
        viewer.addModel(pdb_data, 'pdb')
        
        # Apply basic styling
        viewer.setStyle({'chain': 'A'}, {'stick': {}})
        viewer.setStyle({'resn': 'BOR'}, {'stick': {'colorscheme': 'cyanCarbon'}})
        viewer.zoomTo()
        viewer.setBackgroundColor('white')
        viewer.render()
        
        # Get PNG data
        png_data = viewer.png()
        
        if output_path:
            # Save to file
            with open(output_path, 'wb') as f:
                f.write(png_data)
            print(Fore.GREEN + f"[SUCCESS] Structure image saved to {output_path}" + Style.RESET_ALL)
            return output_path
        else:
            # Return base64 encoded data
            encoded_data = base64.b64encode(png_data).decode('utf-8')
            print(Fore.GREEN + "[SUCCESS] Structure image generated as base64 data" + Style.RESET_ALL)
            return encoded_data
            
    except Exception as e:
        print(Fore.RED + f"[ERROR] Failed to generate structure image: {e}" + Style.RESET_ALL)
        logging.error("Failed to generate structure image: %s", e)
        return None

def create_latex_template():
    """
    Return the LaTeX template used for PDF reporting.

    The enhanced implementation available in ``modules/binding_visualizer/sources``
    already contains a rich, well-tested template.  If that function is available
    we reuse it directly; otherwise we fall back to a minimal local definition.
    """
    parent_sources = os.path.join(
        os.path.dirname(__file__),
        "sources",
    )
    try:
        sys.path.insert(0, parent_sources)
        from binding_visualizer import create_latex_template as sources_template  # type: ignore

        return sources_template()
    except Exception:
        pass
    finally:
        if parent_sources in sys.path:
            sys.path.remove(parent_sources)

    # Minimal fallback template.
    return r"""
\documentclass[11pt,a4paper]{article}
\usepackage[utf8]{inputenc}
\usepackage{geometry}
\usepackage{graphicx}
\usepackage{booktabs}
\usepackage{longtable}
\usepackage{hyperref}
\usepackage{xcolor}
\usepackage{fancyhdr}
\geometry{margin=2cm}

\pagestyle{fancy}
\fancyhf{}
\fancyhead[L]{\textbf{MM Drug Discovery Report}}
\fancyhead[R]{\today}
\fancyfoot[C]{\thepage}

\title{\textbf{\Large Structural Report: {{ pdb_id }}}}
\author{bmyCure4MM Binding Visualizer}
\date{\today}

\begin{document}
\maketitle
\tableofcontents
\newpage

\section{Structure Overview}
\begin{table}[h!]
\centering
\begin{tabular}{ll}
\toprule
\textbf{Property} & \textbf{Value} \\
\midrule
PDB ID & {{ pdb_id }} \\
Method & {{ method or 'N/A' }} \\
Resolution & {{ resolution or 'N/A' }} \\
Chains & {{ chains|length }} \\
Ligands & {{ ligands|join(', ') if ligands else 'None detected' }} \\
\bottomrule
\end{tabular}
\caption{Key structural characteristics}
\end{table}

\section{Visualization}
{% if structure_image_path %}
\begin{figure}[h!]
\centering
\includegraphics[width=0.8\textwidth]{{ "{" }}{{ structure_image_path }}{{ "}" }}
\caption{Rendered structure of {{ pdb_id }}{% if ligand %} with ligand {{ ligand }} highlighted{% endif %}.}
\end{figure}
{% else %}
\noindent\textit{Structure snapshot unavailable.}
{% endif %}

{% if mutations %}
\section{Mutation Highlights}
\begin{longtable}{llll}
\toprule
\textbf{Chain} & \textbf{Residue} & \textbf{Mutation} & \textbf{Effect} \\
\midrule
{% for mutation in mutations %}
{{ mutation.chain }} & {{ mutation.resnum }} & {{ mutation.mutation }} & {{ mutation.effect }} \\
{% endfor %}
\bottomrule
\end{longtable}
{% endif %}

{% if therapies %}
\section{Therapeutic Context}
{% for therapy in therapies %}
\subsection*{{ "{" }}{{ therapy.name }}{{ "}" }}
\begin{itemize}
\item \textbf{Clinical Phase:} {{ therapy.clinical_phase }}
\item \textbf{Mechanism:} {{ therapy.mechanism }}
\item \textbf{Relevance:} {{ therapy.mm_relevance }}
{% if therapy.resistance_mutations %}
\item \textbf{Resistance Mutations:} {{ therapy.resistance_mutations|join(', ') }}
{% endif %}
\end{itemize}
{% endfor %}
{% endif %}

\section{Analysis Notes}
\begin{itemize}
\item Binding site detection: {% if binding_site_detection.enabled %}Enabled ({{ binding_site_detection.cutoff_angstrom }}\,\AA{} cutoff){% else %}Disabled{% endif %}
\item Report generated by bmyCure4MM Binding Visualizer.
\end{itemize}

\section{References}
\begin{enumerate}
\item RCSB Protein Data Bank (\url{https://www.rcsb.org/})
\item py3Dmol: Interactive molecular visualisation
\item bmyCure4MM: Multiple Myeloma research toolkit
\end{enumerate}
\end{document}
"""
def generate_pdf_report(config, pdb_data, pdb_metadata, output_dir=None):
    """
    Generate a PDF report from LaTeX template with molecular structure analysis.
    
    Args:
        config (dict): Configuration dictionary
        pdb_data (str): PDB structure data
        pdb_metadata (dict): Parsed PDB metadata
        output_dir (str): Output directory for PDF (optional)
    
    Returns:
        str: Path to generated PDF file
    """
    print(Fore.CYAN + "[INFO] Generating PDF report from LaTeX..." + Style.RESET_ALL)
    
    if output_dir is None:
        output_dir = this_script_folder_path
    
    pdb_id = config['pdb_id']
    base_filename = f"{pdb_id}_structure_report";
    
    # Generate structure image
    image_path = os.path.join(output_dir, f"{pdb_id}_structure.png")
    structure_image = generate_structure_image(pdb_data, 800, 600, image_path)
    
    # Prepare template variables
    template_vars = {
        'title': f'Molecular Structure Analysis: {pdb_id}',
        'pdb_id': pdb_id,
        'method': pdb_metadata.get('method'),
        'resolution': pdb_metadata.get('resolution'),
        'chains': pdb_metadata.get('chains', []),
        'ligands': pdb_metadata.get('ligands', []),
        'ligand': config.get('ligand'),
        'mutations': config.get('mutations', []),
        'therapies': config.get('therapies', []),
        'binding_site_detection': config.get('binding_site_detection', {}),
        'structure_image_path': image_path if structure_image else None
    }
    
    # Create LaTeX content
    template = Template(create_latex_template())
    latex_content = template.render(**template_vars)
    
    # Write LaTeX file
    tex_path = os.path.join(output_dir, f"{base_filename}.tex")
    with open(tex_path, 'w', encoding='utf-8') as f:
        f.write(latex_content)
    
    print(Fore.GREEN + f"[SUCCESS] LaTeX file generated: {tex_path}" + Style.RESET_ALL)
    
    # Compile to PDF
    pdf_path = compile_latex_to_pdf(tex_path, output_dir)
    
    return pdf_path

def compile_latex_to_pdf(tex_path, output_dir):
    """
    Compile LaTeX file to PDF using pdflatex.
    
    Args:
        tex_path (str): Path to LaTeX file
        output_dir (str): Output directory
    
    Returns:
        str: Path to generated PDF file or None if compilation failed
    """
    print(Fore.CYAN + "[INFO] Compiling LaTeX to PDF..." + Style.RESET_ALL)
    
    # Check if pdflatex is available
    if not shutil.which('pdflatex'):
        print(Fore.YELLOW + "[WARNING] pdflatex not found. Please install LaTeX distribution (e.g., TeX Live, MiKTeX)" + Style.RESET_ALL)
        print(Fore.YELLOW + "[INFO] LaTeX file saved for manual compilation: " + tex_path + Style.RESET_ALL)
        return None
    
    try:
        # Change to output directory
        original_dir = os.getcwd()
        os.chdir(output_dir)
        
        # Get base filename
        base_name = os.path.splitext(os.path.basename(tex_path))[0]
        
        # Compile LaTeX (run twice for proper references)
        for run in range(2):
            print(Fore.CYAN + f"[INFO] Running pdflatex (pass {run + 1}/2)..." + Style.RESET_ALL)
            result = subprocess.run(
                ['pdflatex', '-interaction=nonstopmode', f"{base_name}.tex"],
                capture_output=True,
                text=True,
                timeout=60
            )
            
            if result.returncode != 0:
                print(Fore.RED + f"[ERROR] pdflatex compilation failed:" + Style.RESET_ALL)
                print(Fore.YELLOW + result.stdout + Style.RESET_ALL)
                print(Fore.RED + result.stderr + Style.RESET_ALL)
                os.chdir(original_dir)
                return None
        
        # Return to original directory
        os.chdir(original_dir)
        
        pdf_path = os.path.join(output_dir, f"{base_name}.pdf")
        
        if os.path.exists(pdf_path):
            print(Fore.GREEN + f"[SUCCESS] PDF report generated: {pdf_path}" + Style.RESET_ALL)
            
            # Clean up auxiliary files
            cleanup_extensions = ['.aux', '.log', '.out', '.toc', '.fdb_latexmk', '.fls']
            for ext in cleanup_extensions:
                aux_file = os.path.join(output_dir, f"{base_name}{ext}")
                if os.path.exists(aux_file):
                    os.remove(aux_file)
            
            return pdf_path
        else:
            print(Fore.RED + "[ERROR] PDF file was not generated" + Style.RESET_ALL)
            return None
            
    except subprocess.TimeoutExpired:
        print(Fore.RED + "[ERROR] LaTeX compilation timed out" + Style.RESET_ALL)
        os.chdir(original_dir)
        return None
    except Exception as e:
        print(Fore.RED + f"[ERROR] LaTeX compilation failed: {e}" + Style.RESET_ALL)
        os.chdir(original_dir)
        return None
