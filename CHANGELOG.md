# Changelog

All notable changes to the bmyCure4MM project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added - January 2026

#### 🎯 AI-Powered Decision Support System
- **Automatic simulation interpretation** with qualitative grading (good/caution/bad)
- **Priority-ranked recommendations** system (critical/high/medium/low)
- **6 intelligent scenario detection patterns**:
  - High toxicity (≥30% healthy loss) → automatic critical alert with dose reduction guidance
  - Moderate toxicity (20-30%) → medium priority with monitoring recommendations
  - Poor efficacy (<30% tumor reduction) → high priority suggesting dose escalation or horizon extension
  - Tumor growth (negative reduction) → critical alert recommending regimen switch
  - Early recurrence (<90 days) → medium priority suggesting treatment duration extension
  - Favorable balance → low priority with fine-tuning suggestions
- **Actionable guidance** with specific numerical values (e.g., "Reduce doses by 20-30%", "Extend horizon to 224 days")
- **Smart accordion UI** that auto-opens on warnings/errors, stays closed for favorable results
- **Direct action buttons** ("🔧 Go to Scenario & Implement") linking directly to scenario editor
- **4-step implementation guide** integrated into patient page
- **Bilingual support** (English/Italian) for all recommendations
- **Color-coded assessment** banners (green/yellow/red) for instant status recognition
- Enhanced `_interpret_latest_simulation()` function in `clinic/views.py` (~150 lines)
- Full documentation: [Decision Support System Guide](docs/features/DECISION_SUPPORT_SYSTEM.md)

#### 📊 Interactive Plot Embedding
- **Embed-safe plot endpoint** (`/sim/attempts/<pk>/plot/embed/`) solving CDN-based Plotly.js iframe issues
- **Inline Plotly.js injection** (~3MB) eliminating external CDN dependencies that fail in iframe context
- **Security hardening** with path validation preventing directory traversal attacks
- **Cache-Control headers** (`no-store`) ensuring fresh plot data on each load
- **Lazy loading support** for improved page performance
- New view: `attempt_plot_embed()` in `simulator/views_manage.py`
- Updated URL routing in `simulator/urls.py`
- Modified template: `clinic/templates/clinic/patient_detail.html` iframe now uses embed endpoint

#### 🔍 Debug Logging Infrastructure
- **EmbedDebugMiddleware** for iframe request/response logging
- **Persistent log files** (`logs/embed_debug.log`) for troubleshooting
- **Detailed logging** capturing:
  - Request method, path, headers
  - Response status, headers, content length
  - Request duration timing
  - Automatic filtering for embed-related requests only
- Configuration in `mmportal/settings.py` with dedicated logger and file handler
- Full troubleshooting guide: [Embed Debug Guide](docs/development/EMBED_DEBUG_GUIDE.md)

#### 📚 Documentation Updates
- **README.md** enhanced with new feature descriptions:
  - Added "Interactive Plot Embedding" to Simulator features
  - Added comprehensive "AI-Powered Decision Support System" section with 8 feature bullets
  - Updated Quick Links table with references to new documentation
- **New comprehensive guides**:
  - `docs/features/DECISION_SUPPORT_SYSTEM.md` (500+ lines) covering:
    - System overview and capabilities
    - How interpretation works (thresholds, grading logic)
    - Detailed explanation of all 6 automatic scenarios with examples
    - Recommendation structure and priority system
    - UI behavior and user workflow
    - Implementation guide for users and developers
    - Technical customization options
    - Troubleshooting section
    - Best practices and medical disclaimers
  - `docs/development/EMBED_DEBUG_GUIDE.md` (450+ lines) covering:
    - Grey iframe problem diagnosis
    - Solution architecture (embed-safe endpoint)
    - Debug logging system usage
    - Comprehensive troubleshooting checklist
    - Security considerations (X-Frame-Options, CSP, path validation)
    - Performance optimization strategies
    - Testing procedures (manual and automated)
    - Common error messages reference table
    - Monitoring and metrics recommendations
- **Simulator documentation updates**:
  - `docs/en/simulator.md`: Added "AI-Powered Decision Support" section with usage guide
  - `docs/it/simulatore.md`: Added "Assistente Decisionale AI" section (Italian translation)
  - Both include 6 scenario examples, workflow demonstrations, and medical disclaimers

### Changed - January 2026

#### Clinic Module
- **Enhanced patient detail page** (`clinic/templates/clinic/patient_detail.html`):
  - Replaced simple alert banners with interactive accordion component
  - Conditional aria-expanded and show/collapse logic based on simulation status
  - Priority-based recommendation cards with color coding
  - Highlighted action boxes (blue background) for quick scanning
  - Priority badges (CRITICAL/HIGH) prominently displayed
  - Maintains backward compatibility with existing metrics table
- **Updated patient detail view** (`clinic/views.py`):
  - Modified to pass simulation parameters to interpretation function
  - Enhanced context with `latest_simulation_scenario_url` for direct navigation
  - Improved error handling for missing simulation data

#### Simulator Module
- **Iframe source updated** in patient detail template from CDN endpoint to embed-safe endpoint
- **Added Cache-Control headers** to prevent stale plot rendering
- **Path security validation** in plot serving logic

### Fixed - January 2026
- **Grey iframe issue**: Plots now render correctly in iframe context (no longer grey/blank)
- **CDN Plotly.js loading failure**: Solved by switching to inline JavaScript injection
- **Template syntax errors**: Fixed unclosed {% if %} blocks in accordion implementation
- **Path security vulnerability**: Added validation to prevent directory traversal in plot file serving

### Security - January 2026
- **Path validation** in `attempt_plot_embed()` ensures plot files are within MEDIA_ROOT
- **User authentication** verified for plot access (only owner can view their simulation plots)
- **X-Frame-Options** properly configured (SAMEORIGIN) for iframe security
- **Cache-Control** headers prevent sensitive data caching

---

## [Previous Releases]

### [1.0.0] - November 2024

#### Added
- Initial public release
- Comprehensive documentation structure
- CONTRIBUTING.md with full contributor guidelines
- Architecture guide with system diagrams
- Environment variable template (.env.example)
- Enhanced .gitignore for security

#### Security
- Removed all absolute filesystem paths
- DJANGO_SECRET_KEY now required (no default)
- Secret key validation at startup
- Production-ready security settings

#### Documentation
- Complete README.md rewrite (400+ lines)
- docs/ reorganization with logical structure
- Multiple installation guides
- Quick start documentation
- Use cases for different audiences

---

## Version History Summary

| Version | Date | Key Features |
|---------|------|--------------|
| Unreleased | Jan 2026 | AI Decision Support, Interactive Plot Embedding, Debug Infrastructure |
| 1.0.0 | Nov 2024 | Initial public release, Security hardening, Documentation overhaul |

---

## Breaking Changes

### None in Current Release
All changes in January 2026 update are backward compatible. Existing simulations, scenarios, and data remain fully functional.

---

## Upgrade Guide

### From 1.0.0 to Unreleased (Jan 2026)

1. **Pull Latest Code**:
   ```bash
   git pull origin main
   ```

2. **No Database Migrations Required**: No schema changes in this release

3. **Optional: Clear Browser Cache**: For best experience with new plot embedding
   ```
   Hard refresh: Ctrl+Shift+R (Windows/Linux) or Cmd+Shift+R (Mac)
   ```

4. **Optional: Review New Settings**:
   - No new required environment variables
   - Existing X_FRAME_OPTIONS should be 'SAMEORIGIN' (default)

5. **Restart Application**:
   ```bash
   # If using start_dev.sh
   ./start_dev.sh
   
   # If using Docker
   docker-compose restart
   ```

6. **Test New Features**:
   - Run a simulation and check patient page for decision support accordion
   - Verify plot renders in iframe (not grey box)
   - Check logs/embed_debug.log for iframe request logging

---

## Contributors

This release includes contributions from the development team focused on enhancing user experience, decision support capabilities, and troubleshooting infrastructure.

Special thanks to all contributors who provided feedback on iframe rendering issues and UX improvement suggestions.

---

## Links

- [Decision Support System Documentation](docs/features/DECISION_SUPPORT_SYSTEM.md)
- [Embed Debug Guide](docs/development/EMBED_DEBUG_GUIDE.md)
- [Simulator User Guide](docs/en/simulator.md)
- [Contributing Guidelines](CONTRIBUTING.md)
- [Security Policy](SECURITY.md)

---

**Note**: This project is intended for research and educational purposes. All decision support recommendations are heuristic guides based on mathematical models and should not replace clinical judgment or established medical guidelines.
