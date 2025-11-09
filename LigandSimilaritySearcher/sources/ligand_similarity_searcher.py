"""Command-line interface for the LigandSimilaritySearcher module."""

import os
from pathlib import Path
import logging
import yaml
import pandas as pd
from lib import similarity
from lib import fallback_similarity


BASE_DIR = Path(__file__).resolve().parents[1]
CONFIG_FILE = BASE_DIR / "configs" / "ligand_similarity_searcher.yaml"
DATA_DIR = BASE_DIR / "data"
LOG_FILE = DATA_DIR / "ligand_similarity_searcher.log"


def load_config(path: Path) -> dict:
    """Load YAML configuration file."""
    with path.open() as f:
        return yaml.safe_load(f)


def main() -> None:
    """Run the ligand similarity search workflow."""
    config = load_config(CONFIG_FILE)
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    logging.basicConfig(
        filename=LOG_FILE,
        level=logging.INFO,
        format="%(asctime)s - %(levelname)s - %(message)s",
    )

    try:
        input_smiles = config["input_smiles"]
        fp_type = config.get("fingerprint_type", "morgan")
        radius = config.get("radius", 2)
        n_results = config.get("n_results", 10)
        output_csv = Path(config.get("output_csv", "ligand_similarity_results.csv"))
        if not output_csv.is_absolute():
            output_csv = BASE_DIR / output_csv
        output_csv.parent.mkdir(parents=True, exist_ok=True)

        logging.info("Searching similar compounds for %s", input_smiles)
        ranked = similarity.search_similar_compounds(
            input_smiles, fp_type, radius, n_results
        )

        if not ranked:
            logging.warning("No similar compounds found")
            print("No similar compounds found")
            # Create empty CSV with headers
            df = pd.DataFrame(columns=["CID", "SMILES", "Similarity", "MW", "LogP", "TPSA", "HBA", "HBD"])
            df.to_csv(output_csv, index=False)
            logging.info("Empty results saved to %s", output_csv)
            return

        rows = []
        for item in ranked:
            # Skip items without SMILES
            if not item.get("smiles"):
                logging.warning(f"Skipping CID {item['cid']} - no SMILES available")
                continue
                
            try:
                desc = similarity.compute_descriptors(item["smiles"])
                row = {
                    "CID": item["cid"],
                    "SMILES": item["smiles"],
                    "Similarity": item["similarity"],
                    **desc,
                }
                rows.append(row)
            except Exception as exc:
                logging.warning(f"Failed to compute descriptors for CID {item['cid']}: {exc}")
                continue

        if not rows:
            logging.warning("No valid results after descriptor calculation")
            print("No valid results after descriptor calculation")
            df = pd.DataFrame(columns=["CID", "SMILES", "Similarity", "MW", "LogP", "TPSA", "HBA", "HBD"])
            df.to_csv(output_csv, index=False)
            return

        df = pd.DataFrame(rows)
        df.to_csv(output_csv, index=False)
        logging.info("Results saved to %s", output_csv)
        logging.info(f"Saved {len(rows)} results with full details")
        print(f"Results written to {output_csv}")
        print(f"Found {len(rows)} similar compounds with full data")
    except Exception as exc:
        logging.exception("Search failed: %s", exc)
        print(f"An error occurred: {exc}")
        import traceback
        traceback.print_exc()
        raise  # Re-raise the exception so caller knows it failed


if __name__ == "__main__":
    main()
