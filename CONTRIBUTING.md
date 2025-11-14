# Contributing to bmyCure4MM

Thank you for your interest in contributing to bmyCure4MM! This document provides guidelines and instructions for contributing to the project.

## 🌟 Ways to Contribute

- 🐛 **Report bugs** - Submit detailed bug reports
- ✨ **Suggest features** - Propose new functionality
- 📖 **Improve documentation** - Fix typos, add examples, clarify concepts
- 🔧 **Submit code** - Fix bugs or implement features
- 🧪 **Add tests** - Improve test coverage
- 🎨 **Enhance UI/UX** - Improve the web interface

## 🚀 Getting Started

### 1. Fork and Clone

```bash
# Fork the repository on GitHub, then:
git clone https://github.com/YOUR_USERNAME/bmyCure4MM.git
cd bmyCure4MM
```

### 2. Set Up Development Environment

```bash
# Create virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install dependencies
pip install -r requirements.txt

# Apply migrations
python manage.py migrate

# Create superuser (optional)
python manage.py createsuperuser
```

### 3. Install Redis (for Celery tasks)

**macOS:**
```bash
brew install redis
brew services start redis
```

**Linux:**
```bash
sudo apt-get install redis-server
sudo systemctl start redis
```

**Windows:**
Download Redis from https://github.com/microsoftarchive/redis/releases

### 4. Run Development Server

```bash
# Option 1: All services with one command
./start_dev.sh

# Option 2: Manual start (without Celery)
export CELERY_TASK_ALWAYS_EAGER=True
python manage.py runserver 8001
```

Visit http://127.0.0.1:8001 to see the app running.

## 📝 Development Workflow

### 1. Create a Branch

```bash
git checkout -b feature/your-feature-name
# or
git checkout -b fix/bug-description
```

Use prefixes:
- `feature/` - New functionality
- `fix/` - Bug fixes
- `docs/` - Documentation changes
- `test/` - Test additions/improvements
- `refactor/` - Code refactoring

### 2. Make Your Changes

- Write clear, readable code
- Follow existing code style (PEP 8 for Python)
- Add docstrings to functions and classes
- Update tests if needed
- Update documentation if needed

### 3. Test Your Changes

```bash
# Run all tests
python manage.py test

# Run specific app tests
python manage.py test simulator
python manage.py test chemtools

# Run with verbose output
python manage.py test -v 2

# Check code style (optional)
ruff check .
```

### 4. Commit Your Changes

```bash
git add .
git commit -m "Brief description of changes"
```

**Good commit messages:**
- `Add drug dosing validation to RegimenForm`
- `Fix calculation error in PKPD model`
- `Update installation docs for macOS`
- `Add tests for patient twin functionality`

**Avoid:**
- `Update`
- `Fix stuff`
- `Changes`

### 5. Push and Create Pull Request

```bash
git push origin feature/your-feature-name
```

Then create a Pull Request on GitHub with:
- Clear title describing the change
- Description explaining what and why
- Reference any related issues (`Fixes #123`)
- Screenshots for UI changes

## 🧪 Testing Guidelines

### Writing Tests

```python
from django.test import TestCase
from simulator.models import Scenario

class ScenarioTests(TestCase):
    def setUp(self):
        """Set up test fixtures."""
        self.scenario = Scenario.objects.create(
            title="Test Scenario",
            patient_age=65,
            # ... other required fields
        )
    
    def test_difficulty_calculation(self):
        """Test that difficulty score is calculated correctly."""
        self.assertIsNotNone(self.scenario.difficulty_score)
        self.assertGreaterEqual(self.scenario.difficulty_score, 0)
        self.assertLessEqual(self.scenario.difficulty_score, 100)
```

### Test Coverage

- Aim for >80% coverage on new code
- Test edge cases and error conditions
- Include integration tests for complex features
- Test form validation thoroughly

## 📖 Documentation Guidelines

### Code Documentation

```python
def calculate_clearance(creatinine, age, weight, is_female=False):
    """
    Calculate creatinine clearance using Cockcroft-Gault equation.
    
    Args:
        creatinine (float): Serum creatinine in mg/dL
        age (int): Patient age in years
        weight (float): Body weight in kg
        is_female (bool): True if patient is female
        
    Returns:
        float: Creatinine clearance in mL/min
        
    Raises:
        ValueError: If any parameter is out of valid range
        
    Example:
        >>> calculate_clearance(1.2, 65, 70, is_female=False)
        72.5
    """
    if creatinine <= 0 or age <= 0 or weight <= 0:
        raise ValueError("Parameters must be positive")
    
    clearance = ((140 - age) * weight) / (72 * creatinine)
    if is_female:
        clearance *= 0.85
    
    return clearance
```

### Markdown Documentation

- Use clear headings and structure
- Include code examples
- Add screenshots for UI features
- Keep language simple and accessible
- Use bullet points and tables for clarity

## 🎨 Code Style

### Python Style (PEP 8)

```python
# Good
def calculate_auc(concentrations, times):
    """Calculate area under curve using trapezoidal rule."""
    return np.trapz(concentrations, times)

# Avoid
def calc_auc(c,t): return np.trapz(c,t)
```

### Django Best Practices

- Use model validation in `clean()` methods
- Keep views simple, move logic to models/utils
- Use class-based views when appropriate
- Follow REST conventions for API endpoints

### JavaScript Style

```javascript
// Good
const fetchSimulationResults = async (simulationId) => {
    const response = await fetch(`/api/simulations/${simulationId}/`);
    return response.json();
};

// Avoid
function fetch_sim_res(id){return fetch('/api/simulations/'+id+'/')}
```

## 🐛 Reporting Bugs

### Before Submitting

1. Check existing issues to avoid duplicates
2. Verify it's reproducible on latest version
3. Gather relevant information

### Bug Report Template

```markdown
**Description**
Clear description of the bug.

**To Reproduce**
Steps to reproduce:
1. Go to '...'
2. Click on '...'
3. See error

**Expected Behavior**
What should happen.

**Actual Behavior**
What actually happens.

**Environment**
- OS: [e.g., macOS 14.0]
- Python: [e.g., 3.11.5]
- Django: [e.g., 4.2.7]
- Browser: [e.g., Chrome 120]

**Additional Context**
Screenshots, logs, error messages, etc.
```

## ✨ Feature Requests

### Feature Request Template

```markdown
**Problem Statement**
What problem does this solve?

**Proposed Solution**
How should it work?

**Alternatives Considered**
What other approaches were considered?

**Use Cases**
Who would benefit and how?

**Additional Context**
Mockups, references, examples, etc.
```

## 📚 Resources

- **Documentation**: `/docs/` folder
- **Development Guide**: `/docs/development/DEVELOPMENT.md`
- **API Docs**: `/docs/api.md`
- **Mathematical Models**: `/docs/features/MATHEMATICAL_MODELS_DOCUMENTATION.md`

## 🤝 Code of Conduct

### Our Pledge

We are committed to providing a welcoming and inspiring community for all. We pledge to:
- Use welcoming and inclusive language
- Respect differing viewpoints and experiences
- Accept constructive criticism gracefully
- Focus on what's best for the community
- Show empathy towards others

### Unacceptable Behavior

- Harassment or discriminatory language
- Personal attacks or trolling
- Publishing others' private information
- Other unprofessional conduct

## 📧 Contact

For questions or discussions:
- **Issues**: Use GitHub Issues for bugs and features
- **Security**: Report security issues privately (see SECURITY.md)
- **General**: Start a discussion in GitHub Discussions

## 📄 License

By contributing, you agree that your contributions will be licensed under the MIT License.

---

**Thank you for contributing to bmyCure4MM!** 🎉
