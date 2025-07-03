# bmyCure4MM

A collection of utilities and experiments around drug discovery for multiple
myeloma.  The repository contains small Python scripts used in our research
pipeline together with a set of Jupyter notebooks found under `lab/`.

## Installation

1. Clone the repository
   ```bash
   git clone https://github.com/yourusername/bmyCure4MM.git
   cd bmyCure4MM
   ```
2. Install dependencies
   ```bash
   pip install -r lab/requirements.txt
   ```

## Quick Start

After installing the dependencies, generate the default configuration files:
```bash
python pipelines/processes/settings_generator.py
```
This creates the folder structure under `pipelines/data` and writes
`pipelines/configs/general_settings.json` with absolute paths.

You can then run any of the processing scripts, e.g. to visualise a PDB
structure:
```bash
python pipelines/processes/binding_visualizer.py
```
Outputs and logs are stored under `pipelines/data`.
You can also run the LigandSimilaritySearcher module:
```bash
python LigandSimilaritySearcher/sources/ligand_similarity_searcher.py
```

See the [docs](docs/README.md) directory for a more detailed overview of the
available scripts and notebooks.

## Contributing

We welcome contributions! Please fork the repository, create a branch for your
changes and open a pull request when ready.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE)
file for details.
