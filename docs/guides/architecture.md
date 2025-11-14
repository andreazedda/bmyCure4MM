# bmyCure4MM Architecture Overview

## System Architecture

bmyCure4MM is a Django-based web application for multiple myeloma (MM) clinical decision support and drug discovery research.

## High-Level Architecture

```
┌─────────────────────────────────────────────────────────────┐
│                      Web Interface (Django)                  │
│                    Templates + HTMX + Alpine.js              │
└────────────┬────────────────────────────────┬────────────────┘
             │                                │
             ├────────────────┐              │
             │                │              │
┌────────────▼──────┐  ┌──────▼───────┐  ┌──▼──────────────┐
│   Simulator        │  │  ChemTools   │  │  Clinic         │
│   (PKPD Models)    │  │  (Drug       │  │  (Patient       │
│                    │  │   Discovery) │  │   Management)   │
└────────────┬───────┘  └──────┬───────┘  └──┬──────────────┘
             │                 │              │
             └─────────────────┴──────────────┘
                          │
             ┌────────────▼───────────────┐
             │   Celery Task Queue        │
             │   (Background Processing)  │
             └────────────┬───────────────┘
                          │
             ┌────────────▼───────────────┐
             │   Redis (Message Broker)   │
             └────────────────────────────┘
```

## Core Components

### 1. **Simulator App** 
**Purpose**: PKPD simulation and clinical scenario management

**Key Features**:
- Digital patient twins (virtual patients)
- Treatment regimen optimization
- Clinical scenario builder with validation
- Mathematical pharmacokinetic/pharmacodynamic models
- Difficulty scoring for educational scenarios

**Key Models**:
- `Scenario` - Clinical training scenarios with patient parameters
- `Regimen` - Treatment protocols with drug dosing
- `Simulation` - Individual simulation runs
- `SimulationParameter` - PKPD model parameters

**Technologies**:
- NumPy/SciPy for numerical computations
- Optuna for treatment optimization
- Mathematical models for drug kinetics

### 2. **ChemTools App**
**Purpose**: Drug discovery and molecular analysis tools

**Key Features**:
- Molecular structure visualization (py3Dmol)
- Drug similarity search
- ADME property calculation
- Lipinski's Rule of Five analysis
- PDB structure fetching and analysis
- Chemical structure caching

**Key Models**:
- `DrugParameter` - Pharmacological properties
- `SimilaritySearch` - Drug similarity results

**Technologies**:
- RDKit for cheminformatics
- py3Dmol for 3D visualization
- BioPython for PDB parsing

### 3. **Clinic App**
**Purpose**: Patient data management

**Key Features**:
- Patient records management
- Clinical assessments
- Treatment history tracking
- CRUD operations with validation

**Key Models**:
- `Patient` - Patient demographics and medical history
- `Assessment` - Clinical evaluations

### 4. **Docs Viewer App**
**Purpose**: Multi-language documentation system

**Key Features**:
- Markdown rendering
- Multi-language support (EN/IT)
- Search functionality
- Dynamic documentation loading

## Data Flow

### Simulation Workflow
```
User → Create Scenario → Define Regimen → Configure Parameters
    ↓
Run Simulation → Celery Task → PKPD Model Execution
    ↓
Store Results → Generate Visualizations → Display to User
```

### Drug Discovery Workflow
```
User → Input SMILES/PDB → ChemTools Processing
    ↓
Calculate Properties → Similarity Search → 3D Visualization
    ↓
Cache Results → Display Interactive Viewer
```

### Treatment Optimization
```
Define Patient Twin → Set Optimization Goals → Configure Constraints
    ↓
Optuna Trial → Evaluate Regimen → Score Response
    ↓
Iterate → Find Optimal Protocol → Present Results
```

## Technology Stack

### Backend
- **Framework**: Django 4.x
- **Database**: SQLite (dev), PostgreSQL (production-ready)
- **Task Queue**: Celery with Redis
- **API**: Django REST Framework

### Frontend
- **Templates**: Django Templates with Jinja2
- **Interactivity**: HTMX for dynamic updates
- **JavaScript**: Alpine.js for reactive components
- **CSS**: Custom styles + utility classes
- **Charts**: Chart.js for visualizations

### Scientific Computing
- **NumPy**: Array operations and numerical computing
- **SciPy**: Differential equation solving (odeint)
- **RDKit**: Molecular informatics
- **BioPython**: Biological data structures
- **Optuna**: Bayesian optimization

### Visualization
- **py3Dmol**: 3D molecular structures
- **Matplotlib**: Static plots (converted to base64)
- **Chart.js**: Interactive web charts

## Database Schema

### Core Relationships
```
User ─┬─ Patient ─── Assessment
      │
      ├─ Scenario ──┬─── Regimen
      │             │
      │             └─── Simulation ─── SimulationParameter
      │
      └─ DrugParameter ─── SimilaritySearch
```

## Configuration

### Environment Variables
```bash
# Django
DJANGO_SECRET_KEY          # Secret key (required in production)
DJANGO_DEBUG               # Debug mode (0 or 1)
ALLOWED_HOSTS              # Comma-separated hostnames
CSRF_TRUSTED_ORIGINS       # Comma-separated origins

# Celery
CELERY_BROKER_URL          # Redis connection (redis://localhost:6379/0)
CELERY_RESULT_BACKEND      # Result storage
CELERY_TASK_ALWAYS_EAGER   # Run tasks synchronously (dev)

# Features
PREDLAB_V2                 # Enable experimental features (0 or 1)
```

### File Structure
```
bmyCure4MM/
├── mmportal/              # Project settings
│   ├── settings.py        # Django configuration
│   ├── urls.py            # URL routing
│   └── celery.py          # Celery configuration
├── simulator/             # PKPD simulation app
│   ├── models.py          # Data models
│   ├── views.py           # View logic
│   ├── forms.py           # Form validation
│   ├── difficulty_scoring.py  # Mathematical models
│   └── virtual_patients.py    # Patient archetypes
├── chemtools/             # Drug discovery app
│   ├── models.py
│   ├── views.py
│   ├── tasks.py           # Celery tasks
│   └── utils.py           # Chemical calculations
├── clinic/                # Patient management app
├── docs_viewer/           # Documentation system
├── modules/               # Standalone modules
│   ├── binding_visualizer/
│   └── lipinski_analyzer/
├── pipelines/             # Data processing pipelines
├── lab/                   # Research notebooks
├── templates/             # HTML templates
├── static/                # CSS, JS, images
├── media/                 # User uploads
├── docs/                  # Documentation
└── tests/                 # Test suites
```

## Design Patterns

### Models
- **Fat Models, Thin Views**: Business logic in models
- **Form Validation**: Educational error messages with clinical context
- **Mixins**: Reusable model behaviors

### Views
- **Class-Based Views**: For CRUD operations
- **HTMX Partials**: For dynamic page updates
- **API ViewSets**: For REST endpoints

### Tasks
- **Asynchronous Processing**: Long-running computations
- **Progress Tracking**: Real-time status updates
- **Error Handling**: Graceful failure with logging

## Security Considerations

### Authentication & Authorization
- Django's built-in authentication
- Login required for sensitive operations
- User-scoped data access

### Input Validation
- Form-level validation with `clean()` methods
- Model-level validation
- API serializer validation
- XSS prevention through template escaping

### Secrets Management
- Environment variables for sensitive data
- `.env` files (not committed)
- No hardcoded credentials

## Performance Optimization

### Caching
- Structure caching (PDB files, molecules)
- Query result caching
- Static file compression

### Database
- Indexed fields for frequent queries
- Select/prefetch related for N+1 prevention
- Connection pooling

### Background Processing
- Celery for CPU-intensive tasks
- Redis for fast task queuing
- Result caching

## Monitoring & Logging

### Application Logging
```python
# Activity logging
LOGGING = {
    'loggers': {
        'activity': {...},        # User actions
        'ux.education': {...},    # Educational interactions
        'chemtools.tasks': {...}, # Background tasks
        'celery': {...},          # Task queue
    }
}
```

### Log Files
- `logs/activity.log` - User activity tracking
- `logs/django.log` - Application errors
- `logs/celery_tasks.log` - Background job status

## Deployment

### Development
```bash
./start_dev.sh  # Redis + Celery + Django
```

### Production Checklist
- [ ] Set `DJANGO_DEBUG=0`
- [ ] Configure strong `DJANGO_SECRET_KEY`
- [ ] Set `ALLOWED_HOSTS`
- [ ] Use PostgreSQL instead of SQLite
- [ ] Enable HTTPS
- [ ] Configure static file serving (WhiteNoise/CDN)
- [ ] Set up proper logging/monitoring
- [ ] Configure backup strategy
- [ ] Enable error tracking (Sentry)

## Future Architecture Considerations

### Scalability
- Migrate to PostgreSQL for production
- Add caching layer (Redis/Memcached)
- Load balancing for multiple workers
- Database read replicas

### Microservices (if needed)
- Separate computation service
- Dedicated API gateway
- Message queue (RabbitMQ/Kafka)

### Cloud Deployment
- Container orchestration (Kubernetes)
- Managed services (AWS RDS, ElastiCache)
- CDN for static assets
- Auto-scaling groups

---

For implementation details, see [Development Guide](development/DEVELOPMENT.md).
