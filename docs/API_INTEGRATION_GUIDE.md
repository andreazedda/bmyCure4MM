# Comprehensive API Integration Guide

## Overview

The bmyCure4MM Binding Visualizer now integrates data from **10+ authoritative external APIs** to provide comprehensive, multi-source information for protein-ligand complexes and drug discovery research.

## Integrated API Sources

### 1. **RCSB PDB API** (Protein Data Bank)
- **Base URL**: `https://data.rcsb.org/rest/v1`
- **Data Provided**:
  - Structure metadata and annotations
  - Experimental method and resolution
  - Primary citations and publications
  - Chemical component information
- **Methods**:
  - `get_pdb_summary()` - Comprehensive entry metadata
  - `get_pdb_molecules()` - Ligand/molecule information
  - `get_pdb_citations()` - Publication citations

### 2. **PDBe API** (Protein Data Bank in Europe)
- **Base URL**: `https://www.ebi.ac.uk/pdbe/api`
- **Data Provided**:
  - Structure validation metrics
  - Quality assessment percentiles
  - Binding site details
  - Ligand-protein interactions
- **Methods**:
  - `get_pdbe_validation()` - Quality metrics (clashscore, Ramachandran, etc.)
  - `get_pdbe_ligand_interactions()` - Detailed binding site residues

### 3. **UniProt API**
- **Base URL**: `https://rest.uniprot.org/uniprotkb`
- **Data Provided**:
  - Protein names and functions
  - Gene information
  - Organism details
  - Biological pathways
- **Methods**:
  - `get_uniprot_info()` - Protein information by UniProt ID

### 4. **ChEMBL API** (European Bioinformatics Institute)
- **Base URL**: `https://www.ebi.ac.uk/chembl/api/data`
- **Data Provided**:
  - Drug development status (max phase)
  - Therapeutic classifications
  - Molecular structures (SMILES)
  - Administration routes (oral, parenteral, topical)
  - Safety warnings (black box)
  - First approval dates
  - ATC classifications
- **Methods**:
  - `search_chembl_by_name()` - Basic drug search
  - `get_chembl_drug_details()` - Comprehensive drug information

### 5. **PubChem API** (NCBI)
- **Base URL**: `https://pubchem.ncbi.nlm.nih.gov/rest/pug`
- **Data Provided**:
  - Chemical identifiers (CID)
  - IUPAC names
  - Molecular formulas
  - Molecular weights
  - Chemical structures
- **Methods**:
  - `get_pubchem_compound()` - Compound information by name

### 6. **Reactome API**
- **Base URL**: `https://reactome.org/ContentService`
- **Data Provided**:
  - Biological pathways
  - Pathway hierarchies
  - Species-specific pathways
  - Pathway diagrams and links
- **Methods**:
  - `get_reactome_pathways()` - Pathways for genes/proteins

### 7. **PubMed/NCBI E-utilities**
- **Base URL**: `https://eutils.ncbi.nlm.nih.gov/entrez/eutils`
- **Data Provided**:
  - Scientific publications
  - Article titles and abstracts
  - Authors and journals
  - Publication dates and PMIDs
- **Methods**:
  - `search_pubmed()` - Search publications with automatic relevance ranking

### 8. **ClinicalTrials.gov API v2**
- **Base URL**: `https://clinicaltrials.gov/api/v2`
- **Data Provided**:
  - Active and completed clinical trials
  - Trial phases and status
  - Study titles and descriptions
  - NCT identifiers
- **Methods**:
  - `search_clinical_trials()` - Find trials for drugs/compounds

### 9. **STRING Database** (Protein Interactions)
- **Base URL**: `https://string-db.org/api`
- **Data Provided**:
  - Protein-protein interaction networks
  - Interaction confidence scores
  - Evidence types
  - Known and predicted interactions
- **Methods**:
  - `get_protein_interactions()` - Network data for proteins

### 10. **Additional Resources Linked**
- **DrugBank** - Drug encyclopedic data (web links)
- **KEGG** - Metabolic pathway diagrams (web links)
- **PDB-REDO** - Optimized structure models (web links)
- **PDBsum** - Structural summaries (web links)
- **AlphaFold** - AI-predicted structures (web links)

## Data Flow Architecture

```
User Views PDB Structure (e.g., 4KW5 + Lenalidomide)
                  ↓
    chemtools/views.py: job_detail()
                  ↓
    enrich_pdb_metadata_for_view()
                  ↓
    ┌─────────────────────────────────┐
    │   PDBAPIClient Instance         │
    └─────────────────────────────────┘
                  ↓
    ┌─────────────────────────────────┐
    │  Parallel API Calls (async-safe)│
    ├─────────────────────────────────┤
    │ 1. RCSB PDB: Structure data     │
    │ 2. PDBe: Validation metrics     │
    │ 3. PDBe: Binding interactions   │
    │ 4. ChEMBL: Drug information     │
    │ 5. PubChem: Chemical data       │
    │ 6. PubMed: Publications         │
    │ 7. ClinicalTrials: Active trials│
    │ 8. STRING: Protein network      │
    │ 9. Reactome: Pathways           │
    └─────────────────────────────────┘
                  ↓
    Enriched metadata dictionary
                  ↓
    chemtools/templates/chemtools/job_detail.html
                  ↓
    ┌─────────────────────────────────┐
    │   Display Components            │
    ├─────────────────────────────────┤
    │ • API Sources Banner            │
    │ • Structure Validation Metrics  │
    │ • Ligand Interaction Details    │
    │ • Clinical Trials Table         │
    │ • ChEMBL Drug Information       │
    │ • PubMed Literature             │
    │ • Protein Interaction Network   │
    │ • Similar Structures            │
    │ • Educational Resources         │
    │ • 3D Molecular Viewer           │
    └─────────────────────────────────┘
```

## Key Features

### 1. **Structure Quality Assessment**
- Validation percentiles (compared to all PDB structures)
- Clashscore analysis
- Ramachandran plot outliers
- Overall quality ranking
- Visual progress bars for easy interpretation

### 2. **Binding Site Analysis**
- Detailed residue-level interactions
- Ligand binding site identification
- Key residue highlighting
- Multiple binding site support

### 3. **Drug Development Intelligence**
- Clinical development phase
- Approval status and dates
- Administration routes
- Safety warnings (black box)
- Therapeutic classifications

### 4. **Clinical Context**
- Active clinical trials
- Trial phases and status
- Multiple myeloma specific trials
- Direct links to ClinicalTrials.gov

### 5. **Scientific Literature**
- Relevant PubMed publications
- Automatic relevance ranking
- Author and journal information
- Direct PMID links

### 6. **Protein Network Context**
- Protein-protein interactions
- Interaction confidence scores
- STRING database integration
- Network visualization data

### 7. **Similar Structures**
- Structurally related PDB entries
- Similarity scores
- Alternative binding modes
- Comparative analysis support

## Usage Examples

### Viewing Structure 4KW5 with Lenalidomide

When you view a binding visualization job with PDB ID `4KW5` and ligand `LEN`, the system automatically:

1. Fetches structure metadata from RCSB PDB
2. Retrieves validation metrics from PDBe
3. Gets ligand interaction details
4. Searches ChEMBL for Lenalidomide drug data
5. Finds active clinical trials
6. Retrieves recent PubMed publications
7. Fetches protein interaction network for CRBN
8. Displays all data in organized, expandable sections

### Example API Response Data

```python
{
    'pdb_id': '4KW5',
    'validation_metrics': {
        'overall_quality': 87.3,
        'clashscore': 92.1,
        'ramachandran_outliers': 89.5
    },
    'chembl_drug_details': {
        'chembl_id': 'CHEMBL1336',
        'name': 'LENALIDOMIDE',
        'max_phase': 4,
        'first_approval': 2005,
        'oral': True,
        'black_box_warning': True
    },
    'clinical_trials': [
        {
            'nct_id': 'NCT04876092',
            'title': 'Lenalidomide in Multiple Myeloma',
            'phase': 'PHASE3',
            'status': 'RECRUITING'
        }
    ],
    'pubmed_literature': [
        {
            'pmid': '34567890',
            'title': 'Mechanisms of lenalidomide...',
            'authors': 'Smith J et al.',
            'journal': 'Blood',
            'year': '2024'
        }
    ],
    'protein_network': {
        'count': 10,
        'interactions': [...]
    }
}
```

## Error Handling

All API calls include:
- Timeout protection (10 seconds default)
- Exception handling with logging
- Graceful degradation (missing data doesn't break page)
- Partial data display when some APIs fail
- User-friendly error messages

## Performance Considerations

- **Caching**: Consider implementing Redis/Memcached for API responses
- **Async Calls**: APIs are called sequentially but could be parallelized
- **Rate Limiting**: Some APIs have rate limits (respect them)
- **Timeouts**: All requests have 10-second timeouts
- **Retry Logic**: Could be added for transient failures

## Template Display Logic

The template uses conditional rendering:

```django
{% if pdb_metadata.validation_metrics %}
  <!-- Show validation section -->
{% endif %}

{% if pdb_metadata.clinical_trials %}
  <!-- Show clinical trials table -->
{% endif %}
```

This ensures that:
- Only available data is displayed
- Page doesn't break if API fails
- User sees partial results even if some sources fail
- Loading states can be shown

## Future Enhancements

### Planned Additions:
1. **DrugBank API** - More comprehensive drug data (requires API key)
2. **OpenTargets** - Disease-target associations
3. **ProteinsAPI** - Additional protein annotations
4. **ChemSpider** - Chemical structure search
5. **Caching Layer** - Redis for API response caching
6. **Async Processing** - Celery tasks for slow API calls
7. **User Preferences** - Choose which APIs to query
8. **Export Options** - Download all API data as JSON/CSV

### Performance Improvements:
1. Implement concurrent API calls using `asyncio`
2. Add response caching with TTL
3. Queue slow API calls for background processing
4. Implement API health monitoring
5. Add fallback data sources

## API Rate Limits & Best Practices

| API | Rate Limit | Notes |
|-----|------------|-------|
| RCSB PDB | None specified | Be respectful |
| PDBe | None specified | Be respectful |
| ChEMBL | ~5 req/sec | Monitor usage |
| PubChem | ~5 req/sec | Use POST for bulk |
| PubMed | 3 req/sec (no key) | Register for higher |
| ClinicalTrials | 1000 req/hr | Monitor usage |
| STRING | None specified | Be respectful |

### Best Practices:
1. Always set User-Agent header
2. Implement exponential backoff for retries
3. Cache responses when possible
4. Respect rate limits
5. Handle errors gracefully
6. Log all API failures
7. Monitor API status pages

## Troubleshooting

### Common Issues:

**Issue**: API timeout errors
- **Solution**: Check network connectivity, increase timeout, implement retries

**Issue**: Missing data in display
- **Solution**: Check API logs, verify API endpoints still work, test with curl

**Issue**: Slow page loads
- **Solution**: Implement caching, async calls, or background tasks

**Issue**: Rate limit exceeded
- **Solution**: Implement request queuing, caching, or API key authentication

## Testing

### Manual Testing:
```bash
# Test API client directly
python manage.py shell
>>> from chemtools.pdb_api_client import PDBAPIClient
>>> client = PDBAPIClient()
>>> client.get_pdb_summary('4KW5')
>>> client.get_pdbe_validation('4KW5')
```

### Unit Tests:
```python
# tests/test_pdb_api_client.py
def test_get_pdb_summary():
    client = PDBAPIClient()
    result = client.get_pdb_summary('4KW5')
    assert result is not None
    assert 'struct' in result
```

## Support & Documentation

- **RCSB PDB**: https://data.rcsb.org/
- **PDBe**: https://www.ebi.ac.uk/pdbe/api/doc/
- **ChEMBL**: https://chembl.gitbook.io/chembl-interface-documentation/
- **PubChem**: https://pubchemdocs.ncbi.nlm.nih.gov/pug-rest
- **ClinicalTrials**: https://clinicaltrials.gov/data-api/about-api
- **STRING**: https://string-db.org/help/api/

## License & Attribution

All API data sources are properly attributed in the UI. Users should cite original sources when using this data in publications.

## Change Log

**Version 2.0** (December 2025)
- Added 10+ API integrations
- Structure validation metrics
- Clinical trials data
- Protein interaction networks
- Enhanced literature search
- Comprehensive drug information
- Improved error handling
- Better data visualization

**Version 1.0** (November 2025)
- Initial binding visualizer
- Basic PDB data
- 3D molecular viewer
- Educational resources

---

**Last Updated**: December 11, 2025
**Maintainer**: bmyCure4MM Development Team
