# bmyCure4MM

A collection of utilities and experiments around drug discovery for multiple
myeloma.  The repository contains small Python scripts used in our research
pipeline together with a set of Jupyter notebooks found under `lab/`, and a Django web portal for clinical decision support.

## Installation

1. Clone the repository
   ```bash
   git clone https://github.com/yourusername/bmyCure4MM.git
   cd bmyCure4MM
   ```
2. Install dependencies
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   pip install -r requirements.txt
   ```

## Quick Start

### Option 1: Simple Development Script (Recommended)

Start all services (Redis, Celery, Django) with one command:
```bash
./start_dev.sh
```
This will automatically:
- Start Redis server in background
- Start Celery worker for async job processing
- Start Django development server on http://127.0.0.1:8001
- Stop all services when you press Ctrl+C

**Requirements**: Redis must be installed (`brew install redis` on macOS)

### Option 2: Supervisord (Production-like)

For more control over services:
```bash
# Install supervisor
pip install supervisor

# Start all services
./manage_services.sh start

# Check status
./manage_services.sh status

# View logs
./manage_services.sh logs

# Stop all services
./manage_services.sh stop
```

### Option 3: Docker Compose (Most Portable)

If you have Docker installed:
```bash
docker-compose up -d
```
Access the portal at http://127.0.0.1:8001

### Option 4: Manual Start (Development without Celery)

If you don't want to run Redis/Celery during development:
```bash
export CELERY_TASK_ALWAYS_EAGER=True
python manage.py runserver 8001
```
Jobs will run synchronously (blocking) instead of in background.

### Standalone Scripts

You can also run processing scripts directly:
```bash
python pipelines/processes/settings_generator.py
python pipelines/processes/binding_visualizer.py
python LigandSimilaritySearcher/sources/ligand_similarity_searcher.py
```
Outputs and logs are stored under `pipelines/data`.

See the [docs](docs/README.md) directory for a more detailed overview of the
available scripts and notebooks.

## Contributing

We welcome contributions! Please fork the repository, create a branch for your
changes and open a pull request when ready.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE)
file for details.
