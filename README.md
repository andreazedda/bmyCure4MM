# bmyCure4MM 🧬

> A comprehensive Django-based platform for Multiple Myeloma research, clinical decision support, and drug discovery.

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![Django 4.2+](https://img.shields.io/badge/django-4.2+-green.svg)](https://www.djangoproject.com/)

<img width="1438" height="771" alt="image" src="https://github.com/user-attachments/assets/5b17abc2-8697-445f-9667-73ff2642b58e" />

---

## 🎯 Overview

**bmyCure4MM** is an integrated research platform that combines:

- **🔬 Drug Discovery Tools** - Molecular structure analysis, ADME property calculation, and similarity search
- **💊 PKPD Simulation** - Pharmacokinetic/pharmacodynamic modeling for treatment optimization
- **👤 Digital Patient Twins** - Virtual patient modeling for personalized treatment planning
- **📊 Clinical Decision Support** - Evidence-based scenario management with educational validation
- **🧪 Computational Biology** - Integration with PDB, molecular visualization, and cheminformatics

Perfect for researchers, clinicians, and students working on multiple myeloma treatment optimization and drug discovery.

---

## ✨ Key Features

### 🔬 **ChemTools - Drug Discovery Suite**
- **3D Molecular Visualization**: Interactive py3Dmol-based structure viewer
- **Similarity Search**: Find similar compounds using RDKit fingerprints
- **ADME Prediction**: Calculate drug-like properties (Lipinski's Rule of Five)
- **PDB Integration**: Fetch and analyze protein structures from Protein Data Bank
- **Structure Caching**: Fast retrieval of previously analyzed molecules

### 💊 **Simulator - PKPD Modeling**
- **Treatment Optimization**: Bayesian optimization with Optuna for finding optimal drug regimens
- **Virtual Patients**: 7 clinical archetypes (newly diagnosed, refractory, frail elderly, etc.)
- **Mathematical Models**: Differential equation-based tumor growth and drug response models
- **Clinical Scenarios**: Educational case studies with comprehensive validation
- **Difficulty Scoring**: Auto-calculated complexity for training scenarios
- **📊 NEW: Interactive Plot Embedding**: Inline Plotly visualizations with embed-safe rendering

### 👥 **Clinic - Patient Management**
- Patient record management with comprehensive history tracking
- Clinical assessments tracking (lab values, biomarkers)
- Treatment history documentation
- **🎯 NEW: AI-Powered Decision Support System**:
  - Automatic interpretation of simulation results
  - Priority-ranked actionable recommendations (Critical/High/Medium/Low)
  - Smart accordion UI that opens automatically on warnings
  - Direct links to modify scenarios and implement suggestions
  - Specific numerical guidance ("Reduce doses by 20-30%", "Extend horizon to 224 days")
  - 6 intelligent scenarios: high toxicity, low efficacy, tumor growth, early recurrence, moderate toxicity, favorable balance
  - Step-by-step implementation guides
  - Bilingual support (English/Italian)
- CRUD operations with role-based access

### 📚 **Docs Viewer**
- Multi-language documentation (English/Italian)
- Markdown rendering with syntax highlighting
- Search functionality
- Contextual help system

---

## 🚀 Quick Start (For Learners & New Users)

**Want to learn how to simulate treatments? Start here! 👇**

### 🎓 Option 1: Quick Start with Demo Data (Recommended)

Perfect for learning the platform without creating data manually:

```bash
# Clone the repository
git clone https://github.com/andreazedda/bmyCure4MM.git
cd bmyCure4MM

# Create and activate virtual environment
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# 🎉 ONE-COMMAND SETUP with demo patients!
./quick_start.sh

# Start the server
python manage.py runserver
```

**✅ You now have:**
- ✓ 3 demo patients (high/medium/low risk)
- ✓ 4 pre-filled clinical assessments
- ✓ Admin user (username: `admin`, password: `admin123`)
- ✓ Ready to simulate treatments immediately!

**📖 Next Step:** Visit `http://127.0.0.1:8000` and look for the blue **"🚀 New to the Platform?"** card on the dashboard. It will guide you through your first simulation step-by-step.

---

### 🔧 Option 2: Manual Setup

For production or custom setup:

```bash
# Clone the repository
git clone https://github.com/andreazedda/bmyCure4MM.git
cd bmyCure4MM

# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Set up environment variables (optional)
cp .env.example .env

# IMPORTANT: Generate and set DJANGO_SECRET_KEY
# Generate a secure key:
python -c 'from django.core.management.utils import get_random_secret_key; print(get_random_secret_key())'
# Copy the output and add to .env file:
# DJANGO_SECRET_KEY=<paste-generated-key-here>

# Run database migrations
python manage.py migrate

# Create a superuser (admin)
python manage.py createsuperuser

# Install Redis (if not already installed)
# macOS: brew install redis
# Linux: sudo apt-get install redis-server
# Windows: https://github.com/microsoftarchive/redis/releases
```

### Running the Application

**Option 1: All services with one command (Recommended)**
```bash
./start_dev.sh
```
This starts Redis, Celery, and Django development server automatically.

**Option 2: Without background tasks (Development)**
```bash
export CELERY_TASK_ALWAYS_EAGER=True
python manage.py runserver 8001
```

**Option 3: Docker Compose (Most Portable)**
```bash
docker-compose up -d
```

Visit **http://127.0.0.1:8001** in your browser.

---

## 📖 Documentation

Comprehensive documentation is available in the [`/docs`](docs/) directory:

- **[Installation Guide](docs/guides/installation.md)** - Detailed setup instructions
- **[Architecture Overview](docs/guides/architecture.md)** - System design and components
- **[Development Guide](docs/development/DEVELOPMENT.md)** - Development workflows
- **[API Documentation](docs/api.md)** - REST API reference
- **[Contributing](CONTRIBUTING.md)** - How to contribute to the project

### Quick Links

| Topic | Documentation |
|-------|---------------|
| Getting Started | [Quickstart Guide](docs/en/quickstart.md) |
| PKPD Simulation | [Simulator Docs](docs/en/simulator.md) |
| Patient Twins | [Digital Twin System](docs/en/patient_twin.md) |
| **🆕 Decision Support** | **[AI Assistant Guide](docs/features/DECISION_SUPPORT_SYSTEM.md)** |
| Drug Discovery | [ChemTools Guide](docs/features/CHEMTOOLS_TESTING_SUMMARY.md) |
| Mathematical Models | [PKPD Models](docs/features/MATHEMATICAL_MODELS_DOCUMENTATION.md) |
| Embed Debugging | [Iframe Troubleshooting](docs/development/EMBED_DEBUG_GUIDE.md) |
| Troubleshooting | [Common Issues](docs/development/TROUBLESHOOTING.md) |

---

## 🏗️ Project Structure

```
bmyCure4MM/
├── chemtools/          # Drug discovery and molecular analysis
├── simulator/          # PKPD simulation and optimization
├── clinic/             # Patient management
├── docs_viewer/        # Documentation system
├── modules/            # Standalone analysis modules
│   ├── binding_visualizer/
│   └── lipinski_analyzer/
├── pipelines/          # Data processing pipelines
├── lab/                # Research notebooks and experiments
├── docs/               # Comprehensive documentation
├── templates/          # HTML templates
├── static/             # CSS, JavaScript, images
└── tests/              # Test suites
```

---

## 🔬 Use Cases

### For Researchers
- Evaluate drug candidates using ADME predictions
- Run PKPD simulations to understand drug kinetics
- Optimize treatment protocols using virtual patients
- Analyze protein-ligand interactions

### For Clinicians
- Access evidence-based treatment scenarios
- Understand drug dosing guidelines and contraindications
- Learn about risk stratification and patient assessment
- Review clinical decision support recommendations

### For Educators
- Create educational clinical scenarios with difficulty scoring
- Teach pharmacology through interactive simulations
- Demonstrate drug-target interactions with 3D visualization
- Assess learner understanding through scenario challenges

### For Students
- Learn PKPD concepts through hands-on simulation
- Practice clinical decision-making with virtual patients
- Explore drug discovery through molecular analysis
- Understand mathematical modeling in medicine

---

## 🧪 Testing

```bash
# Run all tests
python manage.py test

# Run specific app tests
python manage.py test simulator
python manage.py test chemtools

# Run with coverage
coverage run --source='.' manage.py test
coverage report

# Run with verbose output
python manage.py test -v 2
```

---

## 🧬 Design / Simulation Report (MVP)

An early design/simulation engine prototype can generate a JSON report with:
- TP/FP tradeoffs + Pareto front
- Bulk vs reservoir dynamics (with subclone evolution)
- Relapse forecast proxy
- Mechanistic toxicity decomposition + rationale

### Generate report (CLI)

```bash
python manage.py generate_design_report --out artifacts/design_report.json
```

### Fetch report (HTTP)

- `GET /api/design-report/`
- Optional query params: `seed`, `steps`, `dt_days`

```bash
curl 'http://127.0.0.1:8001/api/design-report/?seed=123&steps=12&dt_days=7'
```

**Current Test Coverage**: 94% (30/32 tests passing)

---

## 🤝 Contributing

We welcome contributions from the community! Whether you're fixing bugs, adding features, improving documentation, or suggesting ideas, your help is appreciated.

### How to Contribute

1. **Fork the repository**
2. **Create a feature branch** (`git checkout -b feature/amazing-feature`)
3. **Make your changes** and add tests
4. **Commit your changes** (`git commit -m 'Add amazing feature'`)
5. **Push to your branch** (`git push origin feature/amazing-feature`)
6. **Open a Pull Request**

See [CONTRIBUTING.md](CONTRIBUTING.md) for detailed guidelines.

### Development Setup

```bash
# Install dev dependencies
pip install -r requirements.txt

# Run tests before committing
python manage.py test

# Check code style (optional)
ruff check .
```

---

## 📊 Technology Stack

### Backend
- **Django 4.2+** - Web framework
- **Django REST Framework** - API development
- **Celery** - Asynchronous task queue
- **Redis** - Message broker and caching

### Scientific Computing
- **NumPy & SciPy** - Numerical computing and ODE solving
- **RDKit** - Cheminformatics and molecular analysis
- **BioPython** - Biological structure parsing
- **Optuna** - Bayesian optimization
- **py3Dmol** - 3D molecular visualization

### Frontend
- **HTMX** - Dynamic HTML without JavaScript frameworks
- **Alpine.js** - Lightweight reactive components
- **Chart.js** - Interactive charts and visualizations
- **Custom CSS** - Responsive design

### Database
- **SQLite** (development) / **PostgreSQL** (production-ready)

---

## 🔐 Security & Privacy

This repository is prepared for public release with:
- ✅ No hardcoded credentials or API keys
- ✅ Environment variable configuration (see `.env.example`)
- ✅ No absolute filesystem paths
- ✅ Comprehensive `.gitignore` for sensitive files
- ✅ CSRF protection and secure session handling
- ✅ SQL injection prevention through Django ORM

**Important**: Never commit:
- `.env` files with real credentials
- Database files (`*.sqlite3`, `*.db`)
- API keys or secrets
- User data or media uploads

---

## 📄 License

This project is licensed under the **MIT License** - see the [LICENSE](LICENSE) file for details.

```
MIT License

Copyright (c) 2024 Andrea Zedda

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```

---

## 👤 Author

**Andrea Zedda**
- GitHub: [@andreazedda](https://github.com/andreazedda)
- Email: andreazedda@outlook.it

---

## 🙏 Acknowledgments

- Multiple Myeloma research community
- Open-source contributors
- Django and Python communities
- RDKit and BioPython developers
- Clinical collaborators and advisors

---

## 📞 Support & Contact

- **Issues**: Report bugs or request features via [GitHub Issues](https://github.com/andreazedda/bmyCure4MM/issues)
- **Discussions**: Join conversations in [GitHub Discussions](https://github.com/andreazedda/bmyCure4MM/discussions)
- **Documentation**: Check the [docs/](docs/) folder for detailed guides
- **Email**: For private inquiries, contact andreazedda@outlook.it

---

## 🗺️ Roadmap

### Upcoming Features
- [ ] PostgreSQL migration guide
- [ ] Enhanced API documentation with Swagger/OpenAPI
- [ ] Real-time collaboration features
- [ ] Mobile-responsive interface improvements
- [ ] Integration with external drug databases (DrugBank, ChEMBL)
- [ ] Machine learning models for response prediction
- [ ] Multi-user simulation sessions
- [ ] Export/import clinical scenarios

### Research Goals
- [ ] Validate PKPD models against clinical trial data
- [ ] Expand virtual patient library
- [ ] Add more multiple myeloma drugs to library
- [ ] Integrate pharmacogenomics data
- [ ] Develop resistance prediction models

---

## 📈 Project Status

**Active Development** - This project is actively maintained and welcomes contributions.

- **Latest Release**: v1.0.0
- **Last Updated**: November 2024
- **Test Coverage**: 94%
- **Documentation**: Comprehensive

---

**⭐ If you find this project useful, please consider giving it a star on GitHub!**


