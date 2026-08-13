# Testing Guide

## Overview

This project includes comprehensive testing at multiple levels:
- **Unit tests**: pytest for Python backend
- **Integration tests**: Django test client
- **E2E tests**: Playwright for browser automation
- **Accessibility tests**: axe-core integration

## Setup

### Python Tests

```bash
# Install test dependencies
uv sync --frozen --group test

# Run all Python tests
pytest -q

# Run with coverage
pytest --cov=simulator --cov=chemtools --cov-report=term-missing

# Run specific test file
pytest simulator/tests/test_api_help_item_cache.py -v

# Run tests in parallel
pytest -n auto
```

### E2E Tests (Playwright)

```bash
# Install Python Playwright dependencies
uv sync --frozen --group test --group screenshots

# Install Playwright browsers
playwright install chromium

# Run all E2E tests (headless)
pytest tests/e2e/

# Run with visible browser (headed mode)
pytest tests/e2e/ --headed

# Run specific test class
pytest tests/e2e/test_ui_improvements.py::TestDoseColorZones

# Run specific test method
pytest tests/e2e/test_ui_improvements.py::TestDoseColorZones::test_dose_input_shows_green_zone_for_safe_values

# Run with slow motion (500ms delays for debugging)
pytest tests/e2e/ --slowmo 500

# Generate video recordings
pytest tests/e2e/ --video on

# Run with screenshots on failure
pytest tests/e2e/ --screenshot on

# Debug mode with browser devtools
pytest tests/e2e/ --headed --slowmo 1000
```

**Note**: E2E tests use pytest-django's `live_server` fixture (automatically starts test server).

### Old E2E Tests (TypeScript/Node - Deprecated)

For legacy TypeScript E2E tests:

```bash
# Navigate to e2e directory
cd tests/e2e_legacy

# Install Node dependencies
npm install

# Install Playwright browsers
npx playwright install --with-deps

# Run E2E tests
npm run test:e2e
```

## Test Structure

### Python Tests

#### API Tests
- `test_api_help_item_cache.py`: ETag, 304 responses, cache headers
- `test_api_help_search_rank.py`: Search matching, deduplication, cutoff
- `test_api_ux_audit_limits.py`: JSON validation, body truncation, security

#### Integration Tests
- `test_badge_utils_integration.py`: BadgeUtils presence, data-testid attributes

### E2E Tests (Playwright)

#### Help System
- `test_help_and_focus_trap.spec.ts`: Focus management, keyboard navigation, ARIA
- `test_cmdk_debounce_and_search.spec.ts`: Command palette, search, debounce

#### UI Components
- `test_badges_status.spec.ts`: Dose badges, status changes, aria-live announcements

#### Accessibility
- `test_a11y.spec.ts`: axe-core integration, WCAG 2.1 AA compliance

## Key Features Tested

### API Layer
- ✅ ETag caching for help articles
- ✅ 304 Not Modified responses
- ✅ Search result ranking and deduplication
- ✅ 20-result cutoff
- ✅ JSON validation and body truncation (2048 bytes max)
- ✅ Event field truncation (64 chars max)

### JavaScript Utilities
- ✅ `BadgeUtils.computeDoseStatus()` with 'unknown' fallback
- ✅ `BadgeUtils.computeInteractionStatus()`
- ✅ `BadgeUtils.debounce()` utility (150ms default)

### Accessibility
- ✅ aria-live region for screen reader announcements
- ✅ Focus trap in help drawer
- ✅ Focus return after modal close
- ✅ Keyboard navigation (Tab, Shift+Tab, Esc, F1)
- ✅ ARIA attributes (role, aria-modal, aria-controls)
- ✅ data-testid for E2E stability

### UX Features
- ✅ Help drawer with language toggle
- ✅ Command-K search with debounce
- ✅ Inline help buttons with context
- ✅ Preset prioritization in search
- ✅ Type icons (📘 📈 🎯 🧪)

## CI/CD Integration

### GitHub Actions Example

```yaml
name: Tests

on: [push, pull_request]

jobs:
  python-tests:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - uses: actions/setup-python@v4
        with:
          python-version: '3.9'
      - uses: astral-sh/setup-uv@v7
        with:
          version: "0.12.3"
      - run: uv sync --frozen --group test
      - run: uv run pytest --cov --cov-report=xml
      - uses: codecov/codecov-action@v3

  e2e-tests:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3
      - uses: actions/setup-node@v3
      - uses: actions/setup-python@v4
      - uses: astral-sh/setup-uv@v7
        with:
          version: "0.12.3"
      - run: uv sync --frozen --group test --group screenshots
      - run: uv run python manage.py migrate
      - run: cd tests/e2e && npm install
      - run: npx playwright install --with-deps
      - run: npm run test:e2e
```

## Test Data

Tests create explicit model fixtures so scientific inputs remain visible in
the test source.

## Coverage Goals

- **Backend**: >85% coverage on critical paths
- **Accessibility**: Zero critical/serious axe violations
- **E2E**: Key user flows (help, search, form validation)

## Debugging Tips

### Pytest
```bash
# Stop on first failure
pytest -x

# Show print statements
pytest -s

# Verbose output
pytest -vv

# Run only tests matching pattern
pytest -k "test_help"
```

### Playwright
```bash
# Show browser
npx playwright test --headed

# Debug mode (pause execution)
npx playwright test --debug

# Generate test code
npx playwright codegen http://127.0.0.1:8000

# View test report
npx playwright show-report
```

## Known Issues

- TypeScript errors in E2E tests are expected without `npm install`
- Some tests require specific fixtures (Scenario with pk=1)
- Accessibility tests may flag minor contrast issues (non-blocking)

## Contributing

When adding new features:
1. Write unit tests first (TDD)
2. Add integration tests for API endpoints
3. Create E2E tests for user-facing changes
4. Run axe checks on new UI components
5. Update this README with new test categories
