# Quick API Data Sources Reference

## What Data is Available?

When you view a binding visualization (e.g., PDB structure 4KW5 with Lenalidomide), the system automatically fetches data from multiple sources:

## 📊 Data Categories

### 🔬 Structure & Quality
- **Source**: RCSB PDB, PDBe
- **What you get**:
  - Structure resolution and experimental method
  - Quality validation scores (clashscore, Ramachandran)
  - Comparison to all PDB structures (percentile ranking)
  - Binding site residue details
  - Ligand interaction analysis

### 💊 Drug Information
- **Source**: ChEMBL, PubChem
- **What you get**:
  - Drug development phase (Phase 1-4, Approved)
  - Chemical structure (SMILES, molecular formula)
  - Administration routes (oral, injection, topical)
  - Safety warnings (black box warnings)
  - First approval dates
  - Therapeutic classifications

### 🏥 Clinical Trials
- **Source**: ClinicalTrials.gov
- **What you get**:
  - Active and recruiting trials
  - Trial phases and status
  - Study titles and NCT identifiers
  - Direct links to full trial information

### 📚 Scientific Publications
- **Source**: PubMed (NCBI)
- **What you get**:
  - Recent relevant publications
  - Authors, journals, and publication dates
  - PubMed IDs (PMIDs)
  - Direct links to articles

### 🕸️ Protein Networks
- **Source**: STRING Database
- **What you get**:
  - Protein-protein interactions
  - Interaction confidence scores
  - Known and predicted interactions
  - Biological context

### 🔍 Similar Structures
- **Source**: RCSB PDB
- **What you get**:
  - Structurally related PDB entries
  - Similarity scores
  - Alternative binding modes
  - Comparative analysis opportunities

### 🧬 Biological Pathways
- **Source**: Reactome, KEGG
- **What you get**:
  - Metabolic pathways
  - Signal transduction pathways
  - Disease-related pathways
  - Interactive pathway diagrams

## 📍 Where to Find Each Type of Data

### On the Binding Visualizer Page:

1. **Top Banner** (Blue Alert)
   - Shows which API sources have data available
   - Data status indicators

2. **Structure Quality & Validation** (Green Card)
   - Quality percentiles with progress bars
   - Validation metrics explained

3. **Ligand Binding Site Details** (Blue Card)
   - Interacting residues listed
   - Multiple binding sites supported

4. **Clinical Trials** (Yellow Card)
   - Table of active trials
   - Trial phases and status
   - Links to ClinicalTrials.gov

5. **Comprehensive Drug Information** (Purple Card)
   - ChEMBL data
   - Drug properties
   - Administration routes
   - Safety warnings

6. **Recent Publications** (Gray Card)
   - List of PubMed articles
   - Click to read full papers

7. **Protein Interaction Network** (Red Card)
   - Table of interactions
   - Confidence scores

8. **Similar Structures** (Dark Card)
   - Related PDB entries
   - Links to view each structure

## 🎯 Example: Viewing 4KW5 + Lenalidomide

When you view this binding visualization, you'll see:

✅ **Structure Quality**: 87.3% overall quality (better than 87% of all PDB structures)
✅ **Drug Data**: Lenalidomide (Revlimid) - Phase 4, FDA approved, oral administration
✅ **Clinical Trials**: 5+ active trials for multiple myeloma
✅ **Publications**: 10+ recent papers about lenalidomide mechanism
✅ **Protein Network**: CRBN interactions with IKZF1, IKZF3, DDB1
✅ **Pathways**: Ubiquitin-proteasome system, NF-κB signaling

## 🚀 How to Get More Information

1. **Hover** over any element for tooltips
2. **Click** blue links to view external resources
3. **Expand** accordion sections for details
4. **Download** CSV data for offline analysis
5. **Check** the 3D viewer for visual binding analysis

## 💡 Understanding the Data

### Quality Metrics
- **>80%** = Excellent quality
- **50-80%** = Good quality
- **<50%** = Below average (use with caution)

### Drug Development Phases
- **Phase 1**: Safety testing in small groups
- **Phase 2**: Efficacy testing in larger groups
- **Phase 3**: Large-scale efficacy confirmation
- **Phase 4**: Post-market surveillance (approved drug)

### Clinical Trial Status
- **Recruiting**: Actively enrolling patients
- **Active, not recruiting**: Ongoing, no new patients
- **Completed**: Trial finished, results available

### Interaction Scores (STRING)
- **>700**: High confidence
- **400-700**: Medium confidence
- **<400**: Low confidence

## 🔧 Troubleshooting

### "Data loading..." or missing sections?
- Some APIs may be slow or unavailable
- Page shows partial data when available
- Refresh page to retry

### "No data available"?
- API may not have information for this structure
- Try a different PDB ID
- Check if the structure is recent (newly deposited)

### Want more data sources?
- Contact the development team
- Suggest new APIs to integrate
- Check the API Integration Guide for technical details

## 📖 Additional Resources

- **Full API Integration Guide**: See `docs/API_INTEGRATION_GUIDE.md`
- **User Manual**: See `docs/README.md`
- **Technical Documentation**: See `docs/development/`

## 🆘 Support

Having issues? Check:
1. Your internet connection
2. The status banner at the top of the page
3. Browser console for error messages
4. Contact support if problems persist

---

**Quick Tip**: Bookmark interesting structures! Each binding visualization has a unique URL you can save and share.
