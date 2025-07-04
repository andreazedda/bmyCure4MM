# Binding Visualizer Module Documentation

## Overview

The **Binding Visualizer** is a comprehensive Python module designed for interactive visualization of protein-drug binding interactions using PDB structures. It leverages the py3Dmol library to create web-based 3D molecular visualizations with advanced interaction capabilities, specifically designed for multiple myeloma (MM) drug research.

## Features

- **3D Molecular Visualization**: Interactive protein structure visualization using py3Dmol
- **PDB Data Integration**: Automatic fetching and parsing of PDB structures from RCSB database
- **Drug-Target Analysis**: Specialized visualization of drug-protein binding sites
- **Mutation Mapping**: Visual representation of resistance mutations and their effects
- **Therapy Information**: Clinical phase and drug information integration
- **Interactive Controls**: Real-time style modifications and structure exploration
- **PDF Report Generation**: Automatic LaTeX-based PDF reports with structure analysis
- **Offline Support**: Local 3Dmol.js caching for offline use
- **Comprehensive Logging**: Detailed operation tracking and error reporting

## Architecture

### Core Components

1. **PDB Data Handler**: Fetches and processes protein structure data
2. **Visualization Engine**: Creates interactive 3D molecular representations
3. **Configuration Manager**: Handles YAML-based settings and parameters
4. **HTML Generator**: Creates standalone web-based visualizations
5. **Interactive Controls**: JavaScript-based UI for real-time manipulation
6. **PDF Report Generator**: LaTeX-based scientific report creation
7. **Image Generator**: Static molecular structure images for publications

### File Structure

```
binding_visualizer/
├── binding_visualizer.py          # Main module script
├── binding_visualizer.yaml        # Configuration file
├── requirements.txt              # Python dependencies
├── 3Dmol-min.js                 # Local 3Dmol.js library (auto-downloaded)
├── docs/                        # Documentation
│   └── README.md               # This documentation
├── sources/                     # Source code directory
│   └── binding_visualizer.py   # Core implementation
├── *.pdb                       # PDB structure files
├── *_structure_viewer.html     # Generated HTML visualizations
├── *_structure_report.pdf      # Generated PDF reports
├── *_structure_report.tex      # LaTeX source files
├── *_structure.png             # Structure images for PDF
├── debug_*.html               # Debug output files
└── binding_visualizer.log     # Application logs
```

## Configuration

The module is configured via `binding_visualizer.yaml`:

```yaml
# PDB Structure Configuration
pdb_id: "5LF3"                    # Target PDB structure ID

# Viewer Settings
viewer:
  width: 800                      # Visualization width (pixels)
  height: 600                     # Visualization height (pixels)

# Visualization Styles
visualization:
  chain_style:                    # Protein chain visualization
    stick: {}                     # Stick representation
  residue_style:                  # Specific residue styling
    stick: 
      colorscheme: "cyanCarbon"   # Color scheme

# Drug/Ligand Configuration
ligand: "BOR"                     # Main ligand identifier

# PDF Report Generation
generate_pdf: true               # Enable PDF report generation
pdf_options:
  include_structure_image: true  # Include 3D structure image
  image_width: 800              # Image width for PDF
  image_height: 600             # Image height for PDF
  cleanup_aux_files: true       # Remove LaTeX auxiliary files

# Binding Site Analysis
binding_site_detection:
  enabled: true                   # Enable binding site mapping
  cutoff_angstrom: 5.0           # Distance cutoff for binding site

# Mutation Analysis
mutations:
  - chain: K                      # Protein chain
    resnum: 45                    # Residue number
    mutation: A45T                # Mutation notation
    effect: resistance            # Clinical effect

# Therapy Information
therapies:
  - name: Bortezomib             # Drug name
    pdb_ligand: BOR              # PDB ligand code
    clinical_phase: Approved      # Development phase
    mechanism: Proteasome inhibitor # Mechanism of action
```

## Core Functions

### `fetch_pdb_data(pdb_id)`

Retrieves PDB structure data from the RCSB database.

**Parameters:**
- `pdb_id` (str): PDB identifier (e.g., "5LF3")

**Returns:**
- `str`: PDB file content as text

**Features:**
- Automatic timeout handling (10 seconds)
- Comprehensive error logging
- HTTP status code validation
- Network exception handling

### `parse_pdb_header(pdb_data)`

Extracts metadata from PDB file headers.

**Parameters:**
- `pdb_data` (str): Raw PDB file content

**Returns:**
- `dict`: Extracted metadata containing:
  - `method`: Experimental method (X-ray, NMR, etc.)
  - `resolution`: Structure resolution
  - `ligands`: List of non-water ligands
  - `chains`: Available protein chains

### `visualize_structure(pdb_data, width, height, chain_style, residue_style, output_html)`

Creates interactive 3D molecular visualization.

**Parameters:**
- `pdb_data` (str): PDB structure data
- `width` (int): Viewer width in pixels
- `height` (int): Viewer height in pixels
- `chain_style` (dict): Visualization style for protein chains
- `residue_style` (dict): Visualization style for specific residues
- `output_html` (str): Output HTML file path

**Features:**
- Interactive rotation, zoom, and pan controls
- Hover labels for atom identification
- Multiple representation styles (stick, cartoon, surface)
- Color scheme customization
- Automatic structure centering and scaling

### `ensure_local_3dmoljs(local_js_path)`

Manages local 3Dmol.js library for offline visualization.

**Parameters:**
- `local_js_path` (str): Path to local 3Dmol.js file

**Features:**
- Automatic download if library is missing
- File integrity verification
- Network timeout handling
- Fallback error management

### `file_hash(path)`

Calculates SHA256 hash for file integrity verification.

**Parameters:**
- `path` (str): File path

**Returns:**
- `str`: SHA256 hexadecimal hash

### `generate_structure_image(pdb_data, width=800, height=600, output_path=None)`

Generates a static PNG image of the molecular structure for PDF inclusion.

**Parameters:**
- `pdb_data` (str): PDB structure data
- `width` (int): Image width in pixels (default: 800)
- `height` (int): Image height in pixels (default: 600)
- `output_path` (str): Optional file path to save image

**Returns:**
- `str`: Base64 encoded image data or file path

**Features:**
- High-quality PNG generation
- Customizable dimensions
- Base64 encoding for embedding
- File output for external use

### `create_latex_template()`

Creates a comprehensive LaTeX template for molecular structure reports.

**Returns:**
- `str`: LaTeX template with placeholder variables

**Features:**
- Professional scientific document layout
- Automatic table generation for metadata
- Structure image inclusion
- Mutation analysis tables
- Therapy information sections
- Hyperlinked references

### `generate_pdf_report(config, pdb_data, pdb_metadata, output_dir=None)`

Generates a complete PDF report from LaTeX template with molecular analysis.

**Parameters:**
- `config` (dict): Configuration dictionary
- `pdb_data` (str): PDB structure data
- `pdb_metadata` (dict): Parsed PDB metadata
- `output_dir` (str): Output directory (optional)

**Returns:**
- `str`: Path to generated PDF file

**Features:**
- Automated LaTeX content generation
- Structure image embedding
- Mutation and therapy analysis
- Professional formatting
- Automatic PDF compilation

### `compile_latex_to_pdf(tex_path, output_dir)`

Compiles LaTeX source to PDF using pdflatex.

**Parameters:**
- `tex_path` (str): Path to LaTeX source file
- `output_dir` (str): Output directory for PDF

**Returns:**
- `str`: Path to generated PDF or None if failed

**Features:**
- Two-pass compilation for proper references
- Error handling and logging
- Automatic auxiliary file cleanup
- Timeout protection
- Cross-platform compatibility

## Generated HTML Structure

The module creates comprehensive HTML visualizations with:

### Interactive Components

1. **3D Viewer Panel**: Main molecular visualization area
2. **Sidebar Information**: Metadata and controls
3. **Interactive Controls**: Style and view manipulation
4. **Mutation Markers**: Visual indicators for resistance mutations
5. **Therapy Information**: Drug and clinical data display

### JavaScript Functionality

- **Real-time Style Changes**: Dynamic modification of visualization styles
- **Chain Selection**: Toggle visibility of specific protein chains
- **Ligand Highlighting**: Emphasis on drug binding sites
- **Mutation Visualization**: Color-coded resistance mutations
- **Info Panels**: Contextual information display

## Usage Examples

### Basic Structure Visualization

```python
import binding_visualizer as bv

# Load configuration
config = bv.load_config("binding_visualizer.yaml")

# Fetch PDB data
pdb_data = bv.fetch_pdb_data("5LF3")

# Create visualization
bv.visualize_structure(
    pdb_data=pdb_data,
    width=800,
    height=600,
    chain_style={"stick": {}},
    residue_style={"stick": {"colorscheme": "cyanCarbon"}},
    output_html="5LF3_viewer.html"
)
```

### Custom Drug-Target Analysis

```python
# Configure for specific drug analysis
config = {
    "pdb_id": "5LF3",
    "ligand": "BOR",  # Bortezomib
    "mutations": [
        {"chain": "K", "resnum": 45, "mutation": "A45T", "effect": "resistance"}
    ],
    "viewer": {"width": 1200, "height": 800}
}

# Generate enhanced visualization
bv.create_drug_target_visualization(config)
```

### Batch Processing Multiple Structures

```python
pdb_ids = ["5LF3", "4R3O", "5L5B"]

for pdb_id in pdb_ids:
    try:
        pdb_data = bv.fetch_pdb_data(pdb_id)
        output_file = f"{pdb_id}_analysis.html"
        bv.visualize_structure(pdb_data, 800, 600, {}, {}, output_file)
        print(f"✓ Generated visualization for {pdb_id}")
    except Exception as e:
        print(f"✗ Failed to process {pdb_id}: {e}")
```

### PDF Report Generation

```python
# Generate comprehensive PDF report
config = bv.load_config("binding_visualizer.yaml")
pdb_data = bv.fetch_pdb_data(config['pdb_id'])
pdb_metadata = bv.parse_pdb_header(pdb_data)

# Generate PDF report
pdf_path = bv.generate_pdf_report(
    config=config,
    pdb_data=pdb_data,
    pdb_metadata=pdb_metadata,
    output_dir="./reports/"
)

print(f"PDF report generated: {pdf_path}")
```

### Custom LaTeX Template

```python
# Create custom LaTeX template
custom_template = bv.create_latex_template()

# Modify template for specific needs
custom_template = custom_template.replace(
    "\\section{Summary}",
    "\\section{Custom Analysis Section}\\n\\section{Summary}"
)

# Generate report with custom template
# (requires manual template modification)
```

### Structure Image Generation

```python
# Generate high-quality structure image
pdb_data = bv.fetch_pdb_data("5LF3")

# Save as file
image_path = bv.generate_structure_image(
    pdb_data=pdb_data,
    width=1200,
    height=900,
    output_path="5LF3_publication.png"
)

# Get as base64 data for embedding
base64_image = bv.generate_structure_image(
    pdb_data=pdb_data,
    width=800,
    height=600
)
```

## Dependencies

### Python Packages

```bash
pip install -r requirements.txt
```

**Core Dependencies:**
- `py3Dmol==2.4.2`: 3D molecular visualization
- `requests==2.32.3`: HTTP requests for PDB data
- `PyYAML==6.0.2`: Configuration file parsing
- `colorama==0.4.6`: Colored console output
- `biopython==1.85`: Molecular data processing
- `numpy==2.0.2`: Numerical computations
- `jinja2==3.1.2`: LaTeX template rendering
- `matplotlib==3.7.2`: Image processing and generation
- `pillow==10.0.0`: Image manipulation

### System Dependencies

For PDF generation functionality:
- **LaTeX Distribution**: TeX Live (Linux/macOS) or MiKTeX (Windows)
  ```bash
  # Ubuntu/Debian
  sudo apt-get install texlive-latex-extra texlive-fonts-recommended
  
  # macOS (with Homebrew)
  brew install --cask mactex
  
  # Windows: Download and install MiKTeX from https://miktex.org/
  ```

### External Resources

- **3Dmol.js**: Downloaded automatically from https://3Dmol.csb.pitt.edu/
- **RCSB PDB**: Structure data from https://www.rcsb.org/
- **Web Browser**: For viewing generated HTML visualizations

## Output Files

### Generated Visualizations

1. **Main HTML File**: `{pdb_id}_structure_viewer.html`
   - Complete interactive visualization
   - Embedded 3Dmol.js functionality
   - Custom controls and information panels

2. **Debug Files**: `debug_{pdb_id}.html`
   - Intermediate HTML for troubleshooting
   - Useful for development and debugging

3. **Log Files**: `binding_visualizer.log`
   - Timestamped operation logs
   - Error tracking and performance metrics

4. **PDF Reports**: `{pdb_id}_structure_report.pdf`
   - Comprehensive molecular analysis reports
   - Professional LaTeX-formatted documents
   - Structure images and data tables
   - Mutation and therapy analysis

5. **LaTeX Source**: `{pdb_id}_structure_report.tex`
   - LaTeX source files for custom editing
   - Template-based generation
   - Reproducible document creation

6. **Structure Images**: `{pdb_id}_structure.png`
   - High-quality molecular structure images
   - Publication-ready resolution
   - Embedded in PDF reports

### HTML Features

- **Responsive Design**: Adapts to different screen sizes
- **Accessibility**: Screen reader compatible with ARIA labels
- **Performance Optimized**: Efficient rendering for large structures
- **Cross-browser Compatible**: Works in all modern web browsers

## Error Handling

### Network Errors

```python
try:
    pdb_data = fetch_pdb_data("5LF3")
except requests.exceptions.RequestException as e:
    print(f"Network error: {e}")
    # Fallback to local PDB file
    with open("5LF3.pdb", "r") as f:
        pdb_data = f.read()
```

### Configuration Errors

```python
try:
    config = load_config("binding_visualizer.yaml")
except FileNotFoundError:
    print("Configuration file not found. Using defaults.")
    config = get_default_config()
```

### Visualization Errors

```python
try:
    visualize_structure(pdb_data, 800, 600, {}, {}, "output.html")
except Exception as e:
    logging.error(f"Visualization failed: {e}")
    # Generate simplified fallback visualization
    create_fallback_visualization(pdb_data, "output.html")
```

## Best Practices

### Performance Optimization

1. **Local Caching**: Store frequently used PDB files locally
2. **Batch Processing**: Process multiple structures in sequence
3. **Memory Management**: Clear large datasets after processing
4. **Compression**: Use compressed PDB formats when available

### Configuration Management

1. **Version Control**: Track configuration changes
2. **Validation**: Verify configuration parameters before processing
3. **Documentation**: Comment complex visualization settings
4. **Backup**: Maintain backup configurations for critical analyses

### Quality Assurance

1. **Validation**: Verify PDB data integrity
2. **Testing**: Test visualizations across different browsers
3. **Logging**: Monitor performance and error rates
4. **Debugging**: Use debug output files for troubleshooting

## Troubleshooting

### Common Issues

**Issue: 3Dmol.js fails to load**
```
Solution: Check internet connection or use local copy
- Verify 3Dmol-min.js exists in module directory
- Check file permissions and integrity
```

**Issue: PDB fetch timeout**
```
Solution: Increase timeout or use local PDB files
- Download PDB files manually for offline use
- Adjust timeout settings in configuration
```

**Issue: Visualization not rendering**
```
Solution: Check JavaScript console for errors
- Verify browser JavaScript is enabled
- Check for conflicting CSS or JavaScript
```

**Issue: PDF generation fails**
```
Solution: Check LaTeX installation and dependencies
- Install complete LaTeX distribution (TeX Live/MiKTeX)
- Verify pdflatex is in system PATH
- Check LaTeX compilation logs for specific errors
```

**Issue: Structure image generation fails**
```
Solution: Check py3Dmol functionality
- Verify py3Dmol installation and version
- Check if headless display is available
- Try reducing image resolution
```

**Issue: LaTeX compilation timeout**
```
Solution: Optimize LaTeX compilation
- Reduce image size and complexity
- Check for missing LaTeX packages
- Increase timeout in configuration
```

### Debug Mode

Enable detailed debugging by setting environment variable:

```bash
export BINDING_VISUALIZER_DEBUG=1
python binding_visualizer.py
```

This enables:
- Verbose console output
- Additional debug files
- Performance timing information
- Memory usage tracking

## Contributing

### Development Setup

1. Clone the repository
2. Install development dependencies
3. Run tests to verify functionality
4. Follow coding standards and documentation guidelines

### Code Style

- Follow PEP 8 conventions
- Use type hints for function parameters
- Include comprehensive docstrings
- Add unit tests for new features

### Testing

```bash
# Run unit tests
python -m pytest tests/

# Run integration tests
python -m pytest tests/integration/

# Run performance tests
python -m pytest tests/performance/
```

## License

This module is part of the bmyCure4MM project. See the project root for license information.

## Support

For issues, questions, or contributions:

1. Check existing documentation
2. Review troubleshooting section
3. Examine log files for error details
4. Create detailed issue reports with:
   - Configuration used
   - Error messages
   - System information
   - Steps to reproduce

---

*Last updated: July 3, 2025*
*Version: 1.0.0*
*Author: bmyCure4MM Development Team*
