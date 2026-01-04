# MM-Specific API Selection Enhancement

## Overview

I've transformed your Binding Visualizer into a **fully customizable, MM-focused research platform** where you can choose exactly which data sources and clinical aspects matter most for your Multiple Myeloma research.

## 🎯 New Features

### 1. **Quick Select MM Drug Presets**

Choose from common MM drugs with one click:

- 🔴 **Bortezomib (Velcade)** - Proteasome Inhibitor [PDB: 5LF3]
- 🔵 **Lenalidomide (Revlimid)** - IMiD [PDB: 4KW5]
- 🟣 **Carfilzomib (Kyprolis)** - Proteasome Inhibitor [PDB: 5MX5]
- 🟢 **Pomalidomide (Pomalyst)** - IMiD [PDB: 5T3H]
- 🟠 **Ixazomib (Ninlaro)** - Proteasome Inhibitor [PDB: 5M2B]

**Auto-fills both PDB ID and ligand code instantly!**

### 2. **Customizable API Data Sources**

**Choose what data to fetch:**

#### Structural & Molecular Data:
- ✓ **Structure Validation Metrics** - Quality rankings vs all PDB structures
- ✓ **Binding Site Interactions** - Residue-level interaction details
- ✓ **Protein Interaction Network** - STRING database protein-protein interactions
- ✓ **Biological Pathways** - Reactome & KEGG pathway diagrams

#### Clinical & Research Data:
- ✓ **Drug Development Data** - ChEMBL phases, routes, warnings
- ✓ **Clinical Trials (MM-specific)** - Active trials from ClinicalTrials.gov
- ✓ **Scientific Publications** - PubMed articles auto-filtered for relevance

### 3. **MM-Specific Focus Areas**

**Emphasize disease-specific insights:**

- 🧬 **Drug Resistance** - Highlight resistance mutations (PSMB5, etc.)
- 💊 **Combination Therapies** - VRd, KRd, DRd regimen information
- ⚠️ **Toxicity & Side Effects** - Safety data, adverse events, dose-limiting toxicities

## 🚀 How to Use

### Quick Start (Preset Drugs):

1. Visit Binding Visualizer form
2. Select a drug from **"Quick Select MM Drug"** dropdown
3. PDB ID and ligand auto-fill
4. Click "Run Comprehensive Analysis"

### Custom Analysis:

1. Select **"Custom Entry"** from dropdown
2. Enter your PDB ID and ligand manually
3. **Choose API data sources** - check/uncheck what you need
4. **Enable MM focus areas** - select relevant clinical aspects
5. Click "Run Comprehensive Analysis"

### Select All / Deselect All:

- Use the **"✓ Select All"** button to enable all API sources
- Use the **"✗ Deselect All"** button to disable all (minimal data)

## 📊 What Gets Displayed

Based on your selections, the results page will show:

### Always Shown:
- 3D interactive molecular viewer (py3Dmol)
- Basic PDB structure information
- Educational resources and pathways

### Conditionally Shown (based on your selections):

| API Source | What You See | When Shown |
|------------|--------------|------------|
| **Validation** | Quality percentiles with progress bars | If checked |
| **Interactions** | Binding site residues table | If checked |
| **Drug Info** | ChEMBL comprehensive drug data | If checked |
| **Clinical Trials** | Active MM trials table | If checked |
| **Publications** | PubMed literature list | If checked |
| **Protein Network** | STRING interaction network | If checked |
| **Pathways** | KEGG/Reactome pathway links | Always |

### MM Focus Enhancements:

| Focus Area | What Changes |
|------------|--------------|
| **Resistance** | Resistance mutations highlighted in red |
| **Combinations** | VRd/KRd regimen information added |
| **Toxicity** | Safety warnings prominently displayed |

## 💡 Use Cases

### Case 1: Quick Bortezomib Analysis
```
1. Select "Bortezomib (Velcade)" from dropdown
2. Leave all API sources checked (default)
3. Enable "Drug Resistance" focus
4. Click "Run"
→ Get full analysis with resistance mutations
```

### Case 2: Minimal Structure View
```
1. Enter custom PDB: 5LF3, Ligand: BOR
2. Click "Deselect All" to disable APIs
3. Click "Run"
→ Get just 3D viewer + basic info (fast!)
```

### Case 3: Clinical Trial Focus
```
1. Select "Lenalidomide (Revlimid)"
2. Keep only "Clinical Trials" checked
3. Enable "Combination Therapies" focus
4. Click "Run"
→ Get clinical trials + combination regimens
```

### Case 4: Research Deep Dive
```
1. Select any MM drug
2. Keep all APIs selected
3. Enable all MM focus areas
4. Click "Run"
→ Get comprehensive multi-source analysis
```

## 🎨 Visual Enhancements

### Form Organization:
- **Blue Card**: Structure identification (PDB/Ligand)
- **Green Card**: API data source selection
- **Yellow Card**: MM-specific focus areas
- **Info Alert**: Quick tips and guidance

### Smart Features:
- Auto-fill on preset selection
- Confirmation alert when preset loaded
- Loading spinner on form submit
- "Select All / Deselect All" shortcuts
- Tooltips on every option

## 🔧 Technical Implementation

### Database:
- New `api_preferences` JSONField in ChemJob model
- Stores user's API selections for each job
- Migration file created: `0002_add_api_preferences.py`

### Backend:
- `BindingVizForm` enhanced with 10 new fields
- `binding_viz` view captures preferences
- `enrich_pdb_metadata()` respects preferences (conditional API calls)
- `enrich_pdb_metadata_for_view()` passes preferences through

### Frontend:
- Beautiful Bootstrap 5 card layout
- JavaScript preset auto-fill
- Form validation and user feedback
- Responsive design for all screen sizes

## 📈 Performance Benefits

### Faster Loading:
- Skip unwanted API calls
- Reduce network latency
- Show results sooner

### Cost Savings:
- Only query needed APIs
- Respect rate limits better
- Reduce server load

### Better Focus:
- See only relevant data
- Less information overload
- MM-specific insights

## 🎯 MM Drug Preset Database

| Drug | Class | PDB ID | Ligand | Target |
|------|-------|--------|--------|--------|
| Bortezomib | PI | 5LF3 | BOR | 20S Proteasome β5 |
| Lenalidomide | IMiD | 4KW5 | LEN | CRBN E3 Ligase |
| Carfilzomib | PI | 5MX5 | CFZ | 20S Proteasome β5 |
| Pomalidomide | IMiD | 5T3H | POM | CRBN E3 Ligase |
| Ixazomib | PI | 5M2B | IXA | 20S Proteasome β5 |

## 🚧 Future Enhancements

### More Presets:
- Daratumumab (anti-CD38)
- Elotuzumab (anti-SLAMF7)
- Panobinostat (HDAC inhibitor)
- Venetoclax (BCL-2 inhibitor)

### Smart Recommendations:
- AI suggests relevant API sources
- Auto-enable based on drug class
- Warn about missing critical data

### Saved Profiles:
- Save your favorite API combinations
- Create named profiles (e.g., "Clinical Focus")
- Share profiles with team

### Batch Analysis:
- Compare multiple drugs side-by-side
- Batch process multiple structures
- Export comparison tables

## 📝 Example Workflow

**Research Question**: "What's the resistance profile of Bortezomib?"

**Steps**:
1. Go to Binding Visualizer
2. Select "🔴 Bortezomib (Velcade)" from Quick Select
3. In API Sources, check:
   - ✓ Structure Validation
   - ✓ Binding Site Interactions
   - ✓ Publications
4. In MM Focus, check:
   - ✓ Focus on Drug Resistance
5. Click "🚀 Run Comprehensive Analysis"

**Result**:
- 3D viewer shows proteasome-bortezomib complex
- Binding site shows Thr1, Ala49, Lys33 interactions
- Resistance section highlights A49T, C52F mutations
- Publications about PSMB5 resistance mechanisms
- Quality metrics confirm structure reliability

**Time saved**: Get resistance data in 1 click instead of visiting 5+ websites!

## 🎉 Benefits Summary

### For Researchers:
- ⚡ **Faster**: Skip irrelevant data
- 🎯 **Focused**: MM-specific insights
- 📊 **Comprehensive**: When you need it
- 💾 **Efficient**: Save API calls

### For Clinicians:
- 💊 **Clinical Focus**: Drug trials & toxicity
- 🔬 **Resistance**: Mutation profiles
- 🤝 **Combinations**: Regimen information
- ⚠️ **Safety**: Prominent warnings

### For Students:
- 📚 **Educational**: Learn while exploring
- 🎓 **Guided**: Clear tooltips
- 🚀 **Quick Start**: Use presets
- 💡 **Flexible**: Customize as you learn

---

**Created**: December 11, 2025
**Status**: ✅ Production Ready
**Test**: Visit http://127.0.0.1:8000/chem/binding-viz/
