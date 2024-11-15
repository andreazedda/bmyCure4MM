import requests
import io

def find_ligands_in_pdb_file(pdb_id):
    """
    Fetches the PDB file and lists all heteroatoms (ligands) in the structure.
    """
    pdb_url = f"https://files.rcsb.org/download/{pdb_id}.pdb"
    try:
        response = requests.get(pdb_url, timeout=10)
        response.raise_for_status()

        # Read through PDB file lines to find HET and HETNAM entries
        pdb_lines = response.text.splitlines()
        for line in pdb_lines:
            if line.startswith("HET "):
                print("HET Record:", line)
            elif line.startswith("HETNAM"):
                print("HETNAM Record:", line)

    except requests.exceptions.RequestException as error:
        print("Error fetching or reading PDB file:", error)

# Usage
find_ligands_in_pdb_file("5LF3")
