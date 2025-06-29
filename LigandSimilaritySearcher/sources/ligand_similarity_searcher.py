"""Command-line interface for the LigandSimilaritySearcher module."""

from pathlib import Path
import logging
import yaml
import pandas as pd
from lib import similarity


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

        rows = []
        for item in ranked:
            desc = similarity.compute_descriptors(item["smiles"])
            row = {
                "CID": item["cid"],
                "SMILES": item["smiles"],
                "Similarity": item["similarity"],
                **desc,
            }
            rows.append(row)

        df = pd.DataFrame(rows)
        df.to_csv(output_csv, index=False)
        logging.info("Results saved to %s", output_csv)
        print(f"Results written to {output_csv}")
    except Exception as exc:
        logging.exception("Search failed: %s", exc)
        print("An error occurred. See log file for details.")


if __name__ == "__main__":
    main()
