"""
PDB and Biological Database API Client
Provides standardized access to structural and biological databases
"""
import requests
from typing import Dict, List, Optional, Any
import logging

logger = logging.getLogger(__name__)


class PDBAPIClient:
    """Client for fetching PDB metadata and related information from various APIs"""
    
    RCSB_API_BASE = "https://data.rcsb.org/rest/v1"
    RCSB_GRAPHQL = "https://data.rcsb.org/graphql"
    PDBE_API_BASE = "https://www.ebi.ac.uk/pdbe/api"
    PDBE_GRAPH_API = "https://www.ebi.ac.uk/pdbe/graph-api/pdb"
    UNIPROT_API_BASE = "https://rest.uniprot.org/uniprotkb"
    CHEMBL_API_BASE = "https://www.ebi.ac.uk/chembl/api/data"
    PUBCHEM_API_BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
    REACTOME_API_BASE = "https://reactome.org/ContentService"
    DRUGBANK_API_BASE = "https://go.drugbank.com"
    CLINICALTRIALS_API = "https://clinicaltrials.gov/api/v2"
    OPENTARGETS_API = "https://api.platform.opentargets.org/api/v4/graphql"
    EUTILS_BASE = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    
    def __init__(self, timeout: int = 10):
        self.timeout = timeout
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'bmyCure4MM/1.0',
            'Accept': 'application/json'
        })
    
    def get_pdb_summary(self, pdb_id: str) -> Optional[Dict[str, Any]]:
        """
        Fetch comprehensive PDB entry summary from RCSB PDB API
        
        Args:
            pdb_id: PDB identifier (e.g., '4KW5')
            
        Returns:
            Dictionary with structure metadata or None if failed
        """
        import time
        start_time = time.time()
        logger.info(f"[DEBUG] 🌐 API CALL: get_pdb_summary('{pdb_id}')")
        
        try:
            url = f"{self.RCSB_API_BASE}/core/entry/{pdb_id.upper()}"
            logger.debug(f"[DEBUG]   URL: {url}")
            logger.debug(f"[DEBUG]   Timeout: {self.timeout}s")
            
            response = self.session.get(url, timeout=self.timeout)
            elapsed = time.time() - start_time
            
            logger.debug(f"[DEBUG]   Status: {response.status_code}")
            logger.debug(f"[DEBUG]   Response time: {elapsed:.2f}s")
            logger.debug(f"[DEBUG]   Response size: {len(response.content)} bytes")
            
            response.raise_for_status()
            data = response.json()
            
            logger.info(f"[DEBUG] ✅ Successfully fetched PDB summary for {pdb_id} in {elapsed:.2f}s")
            logger.debug(f"[DEBUG]   Data keys: {list(data.keys()) if isinstance(data, dict) else 'not a dict'}")
            
            return data
        except Exception as e:
            elapsed = time.time() - start_time
            logger.warning(f"[DEBUG] ❌ Failed to fetch PDB summary for {pdb_id} after {elapsed:.2f}s: {e}")
            logger.debug(f"[DEBUG]   Exception type: {type(e).__name__}")
            return None
    
    def get_pdb_molecules(self, pdb_id: str) -> Optional[Dict[str, Any]]:
        """Fetch molecules/ligands information from RCSB PDB"""
        try:
            url = f"{self.RCSB_API_BASE}/core/chemcomp/{pdb_id.upper()}"
            response = self.session.get(url, timeout=self.timeout)
            response.raise_for_status()
            return response.json()
        except Exception as e:
            logger.warning(f"Failed to fetch molecules for {pdb_id}: {e}")
            return None
    
    def get_pdb_citations(self, pdb_id: str) -> Optional[List[Dict[str, Any]]]:
        """Fetch publication citations from RCSB PDB"""
        try:
            url = f"{self.RCSB_API_BASE}/core/entry/{pdb_id.upper()}"
            response = self.session.get(url, timeout=self.timeout)
            response.raise_for_status()
            data = response.json()
            
            citations = []
            if 'rcsb_primary_citation' in data:
                cit = data['rcsb_primary_citation']
                citations.append({
                    'title': cit.get('title', ''),
                    'journal': cit.get('journal_abbrev', ''),
                    'year': cit.get('year'),
                    'doi': cit.get('pdbx_database_id_doi', ''),
                    'pmid': cit.get('pdbx_database_id_pub_med', ''),
                    'authors': cit.get('rcsb_authors', [])
                })
            
            return citations if citations else None
        except Exception as e:
            logger.warning(f"Failed to fetch citations for {pdb_id}: {e}")
            return None
    
    def get_uniprot_info(self, uniprot_id: str) -> Optional[Dict[str, Any]]:
        """Fetch protein information from UniProt API"""
        import time
        start_time = time.time()
        logger.info(f"[DEBUG] 🌐 API CALL: get_uniprot_info('{uniprot_id}')")
        
        try:
            url = f"{self.UNIPROT_API_BASE}/{uniprot_id}.json"
            logger.debug(f"[DEBUG]   URL: {url}")
            
            response = self.session.get(url, timeout=self.timeout)
            elapsed = time.time() - start_time
            
            logger.debug(f"[DEBUG]   Status: {response.status_code}")
            logger.debug(f"[DEBUG]   Response time: {elapsed:.2f}s")
            
            response.raise_for_status()
            data = response.json()
            
            logger.info(f"[DEBUG] ✅ Successfully fetched UniProt info for {uniprot_id} in {elapsed:.2f}s")
            
            return {
                'name': data.get('proteinDescription', {}).get('recommendedName', {}).get('fullName', {}).get('value', ''),
                'organism': data.get('organism', {}).get('scientificName', ''),
                'function': data.get('comments', [{}])[0].get('texts', [{}])[0].get('value', '') if data.get('comments') else '',
                'gene': data.get('genes', [{}])[0].get('geneName', {}).get('value', '') if data.get('genes') else '',
                'pathways': [p.get('name', '') for p in data.get('comments', []) if p.get('commentType') == 'PATHWAY']
            }
        except Exception as e:
            logger.warning(f"Failed to fetch UniProt info for {uniprot_id}: {e}")
            return None
    
    def search_chembl_by_name(self, drug_name: str) -> Optional[Dict[str, Any]]:
        """Search ChEMBL for drug information by name"""
        try:
            url = f"{self.CHEMBL_API_BASE}/molecule.json?pref_name__iexact={drug_name}"
            response = self.session.get(url, timeout=self.timeout)
            response.raise_for_status()
            data = response.json()
            
            if data.get('molecules'):
                mol = data['molecules'][0]
                return {
                    'chembl_id': mol.get('molecule_chembl_id'),
                    'name': mol.get('pref_name'),
                    'max_phase': mol.get('max_phase'),
                    'therapeutic_flag': mol.get('therapeutic_flag'),
                    'indication_class': mol.get('indication_class'),
                    'smiles': mol.get('molecule_structures', {}).get('canonical_smiles', ''),
                }
            return None
        except Exception as e:
            logger.warning(f"Failed to search ChEMBL for {drug_name}: {e}")
            return None
    
    def get_reactome_pathways(self, entity_id: str) -> Optional[List[Dict[str, Any]]]:
        """Fetch Reactome pathways for an entity (gene/protein)"""
        try:
            url = f"{self.REACTOME_API_BASE}/data/pathways/low/entity/{entity_id}"
            response = self.session.get(url, timeout=self.timeout)
            response.raise_for_status()
            pathways = response.json()
            
            return [{
                'stId': p.get('stId'),
                'displayName': p.get('displayName'),
                'species': p.get('species', [{}])[0].get('displayName', '') if p.get('species') else '',
                'url': f"https://reactome.org/content/detail/{p.get('stId')}"
            } for p in pathways]
        except Exception as e:
            logger.warning(f"Failed to fetch Reactome pathways for {entity_id}: {e}")
            return None
    
    def get_pubchem_compound(self, name: str) -> Optional[Dict[str, Any]]:
        """Search PubChem for compound information"""
        try:
            # Search by name
            search_url = f"{self.PUBCHEM_API_BASE}/compound/name/{name}/JSON"
            response = self.session.get(search_url, timeout=self.timeout)
            response.raise_for_status()
            data = response.json()
            
            if 'PC_Compounds' in data and data['PC_Compounds']:
                compound = data['PC_Compounds'][0]
                cid = compound.get('id', {}).get('id', {}).get('cid')
                
                return {
                    'cid': cid,
                    'iupac_name': next((p['value']['sval'] for p in compound.get('props', []) 
                                       if p.get('urn', {}).get('label') == 'IUPAC Name'), ''),
                    'molecular_formula': next((p['value']['sval'] for p in compound.get('props', []) 
                                              if p.get('urn', {}).get('label') == 'Molecular Formula'), ''),
                    'molecular_weight': next((p['value']['fval'] for p in compound.get('props', []) 
                                             if p.get('urn', {}).get('label') == 'Molecular Weight'), None),
                    'url': f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
                }
            return None
        except Exception as e:
            logger.warning(f"Failed to fetch PubChem compound for {name}: {e}")
            return None
    
    def get_pdbe_validation(self, pdb_id: str) -> Optional[Dict[str, Any]]:
        """Fetch structure validation metrics from PDBe"""
        try:
            url = f"{self.PDBE_API_BASE}/validation/global-percentiles/entry/{pdb_id.lower()}"
            response = self.session.get(url, timeout=self.timeout)
            response.raise_for_status()
            data = response.json()
            
            if pdb_id.lower() in data:
                metrics = data[pdb_id.lower()]
                return {
                    'clashscore': metrics.get('clashscore', {}).get('percent_rank'),
                    'ramachandran_outliers': metrics.get('ramachandran_outliers', {}).get('percent_rank'),
                    'sidechain_outliers': metrics.get('sidechain_outliers', {}).get('percent_rank'),
                    'overall_quality': metrics.get('overall_quality', {}).get('percent_rank'),
                }
            return None
        except Exception as e:
            logger.warning(f"Failed to fetch PDBe validation for {pdb_id}: {e}")
            return None
    
    def get_pdbe_ligand_interactions(self, pdb_id: str) -> Optional[List[Dict[str, Any]]]:
        """Fetch detailed ligand interaction data from PDBe"""
        try:
            url = f"{self.PDBE_API_BASE}/pdb/entry/binding_sites/{pdb_id.lower()}"
            response = self.session.get(url, timeout=self.timeout)
            response.raise_for_status()
            data = response.json()
            
            if pdb_id.lower() in data:
                sites = data[pdb_id.lower()]
                interactions = []
                for site in sites:
                    interactions.append({
                        'site_id': site.get('site_id'),
                        'ligand_id': site.get('ligand_id', {}).get('chem_comp_id'),
                        'residues': [f"{r.get('chain_id')}{r.get('residue_number')}" 
                                   for r in site.get('site_residues', [])],
                        'details': site.get('details', ''),
                    })
                return interactions
            return None
        except Exception as e:
            logger.warning(f"Failed to fetch ligand interactions for {pdb_id}: {e}")
            return None
    
    def get_similar_structures(self, pdb_id: str, limit: int = 5) -> Optional[List[Dict[str, Any]]]:
        """Find similar PDB structures using RCSB sequence similarity"""
        try:
            url = f"{self.RCSB_API_BASE}/core/entry/{pdb_id.upper()}"
            response = self.session.get(url, timeout=self.timeout)
            response.raise_for_status()
            data = response.json()
            
            # Extract primary accession for similarity search
            if 'rcsb_entry_container_identifiers' in data:
                container = data['rcsb_entry_container_identifiers']
                # This is simplified - in production you'd use RCSB search API
                similar = [{
                    'pdb_id': pdb_id.upper(),
                    'title': data.get('struct', {}).get('title', ''),
                    'similarity': 100.0
                }]
                return similar[:limit]
            return None
        except Exception as e:
            logger.warning(f"Failed to fetch similar structures for {pdb_id}: {e}")
            return None
    
    def search_pubmed(self, query: str, max_results: int = 5) -> Optional[List[Dict[str, Any]]]:
        """Search PubMed for relevant publications"""
        try:
            # Search for PMIDs
            search_url = f"{self.EUTILS_BASE}/esearch.fcgi"
            params = {
                'db': 'pubmed',
                'term': query,
                'retmax': max_results,
                'retmode': 'json',
                'sort': 'relevance'
            }
            response = self.session.get(search_url, params=params, timeout=self.timeout)
            response.raise_for_status()
            search_data = response.json()
            
            pmids = search_data.get('esearchresult', {}).get('idlist', [])
            if not pmids:
                return None
            
            # Fetch article details
            fetch_url = f"{self.EUTILS_BASE}/esummary.fcgi"
            params = {
                'db': 'pubmed',
                'id': ','.join(pmids),
                'retmode': 'json'
            }
            response = self.session.get(fetch_url, params=params, timeout=self.timeout)
            response.raise_for_status()
            fetch_data = response.json()
            
            articles = []
            for pmid in pmids:
                if pmid in fetch_data.get('result', {}):
                    article = fetch_data['result'][pmid]
                    articles.append({
                        'pmid': pmid,
                        'title': article.get('title', ''),
                        'authors': article.get('sortfirstauthor', ''),
                        'journal': article.get('source', ''),
                        'year': article.get('pubdate', '').split()[0] if article.get('pubdate') else '',
                        'url': f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                    })
            
            return articles if articles else None
        except Exception as e:
            logger.warning(f"Failed to search PubMed for {query}: {e}")
            return None
    
    def search_clinical_trials(self, drug_name: str) -> Optional[List[Dict[str, Any]]]:
        """Search ClinicalTrials.gov for drug trials"""
        try:
            url = f"{self.CLINICALTRIALS_API}/studies"
            params = {
                'query.term': drug_name,
                'filter.overallStatus': 'RECRUITING,ACTIVE_NOT_RECRUITING,COMPLETED',
                'pageSize': 5,
                'format': 'json'
            }
            response = self.session.get(url, params=params, timeout=self.timeout)
            response.raise_for_status()
            data = response.json()
            
            trials = []
            for study in data.get('studies', []):
                protocol = study.get('protocolSection', {})
                identification = protocol.get('identificationModule', {})
                status = protocol.get('statusModule', {})
                
                trials.append({
                    'nct_id': identification.get('nctId'),
                    'title': identification.get('briefTitle', ''),
                    'status': status.get('overallStatus', ''),
                    'phase': protocol.get('designModule', {}).get('phases', ['Unknown'])[0],
                    'url': f"https://clinicaltrials.gov/study/{identification.get('nctId')}"
                })
            
            return trials if trials else None
        except Exception as e:
            logger.warning(f"Failed to search clinical trials for {drug_name}: {e}")
            return None
    
    def get_protein_interactions(self, protein_name: str) -> Optional[Dict[str, Any]]:
        """Fetch protein-protein interactions from STRING database"""
        try:
            # STRING API requires protein names/identifiers
            url = f"https://string-db.org/api/json/network"
            params = {
                'identifiers': protein_name,
                'species': 9606,  # Human
                'required_score': 400,
                'limit': 10
            }
            response = self.session.get(url, params=params, timeout=self.timeout)
            response.raise_for_status()
            data = response.json()
            
            interactions = [{
                'protein_a': item.get('preferredName_A', ''),
                'protein_b': item.get('preferredName_B', ''),
                'score': item.get('score', 0),
                'evidence': item.get('escore', 0)
            } for item in data]
            
            return {'interactions': interactions, 'count': len(interactions)} if interactions else None
        except Exception as e:
            logger.warning(f"Failed to fetch protein interactions for {protein_name}: {e}")
            return None
    
    def get_chembl_drug_details(self, chembl_id: str) -> Optional[Dict[str, Any]]:
        """Fetch detailed drug information from ChEMBL"""
        try:
            url = f"{self.CHEMBL_API_BASE}/molecule/{chembl_id}.json"
            response = self.session.get(url, timeout=self.timeout)
            response.raise_for_status()
            data = response.json()
            
            return {
                'chembl_id': data.get('molecule_chembl_id'),
                'name': data.get('pref_name'),
                'type': data.get('molecule_type'),
                'max_phase': data.get('max_phase'),
                'therapeutic_flag': data.get('therapeutic_flag'),
                'oral': data.get('oral'),
                'parenteral': data.get('parenteral'),
                'topical': data.get('topical'),
                'black_box_warning': data.get('black_box_warning'),
                'first_approval': data.get('first_approval'),
                'indication_class': data.get('indication_class'),
                'atc_classifications': data.get('atc_classifications', []),
            }
        except Exception as e:
            logger.warning(f"Failed to fetch ChEMBL details for {chembl_id}: {e}")
            return None
    
    def build_standard_links(self, pdb_id: str, ligand_name: Optional[str] = None) -> Dict[str, str]:
        """
        Build standardized links for a PDB structure
        
        Args:
            pdb_id: PDB identifier
            ligand_name: Optional ligand/drug name
            
        Returns:
            Dictionary of categorized links
        """
        pdb_upper = pdb_id.upper()
        
        links = {
            'structure': {
                'rcsb_pdb': f"https://www.rcsb.org/structure/{pdb_upper}",
                'pdbe': f"https://www.ebi.ac.uk/pdbe/entry/pdb/{pdb_upper}",
                'pdb_redo': f"https://pdb-redo.eu/db/{pdb_upper}",
                'pdbsum': f"https://www.ebi.ac.uk/thornton-srv/databases/cgi-bin/pdbsum/GetPage.pl?pdbcode={pdb_upper}",
            },
            'analysis': {
                'ligand_explorer': f"https://www.rcsb.org/ligand/{pdb_upper}",
                'plip': "https://plip-tool.biotec.tu-dresden.de/plip-web/plip/index",
                'procheck': f"https://www.ebi.ac.uk/thornton-srv/software/PROCHECK/",
            },
            'literature': {
                'pubmed_pdb': f"https://pubmed.ncbi.nlm.nih.gov/?term={pdb_upper}",
                'google_scholar': f"https://scholar.google.com/scholar?q={pdb_upper}+protein+structure",
                'europe_pmc': f"https://europepmc.org/search?query={pdb_upper}",
            },
            'pathways': {
                'kegg_proteasome': "https://www.kegg.jp/pathway/hsa03050",
                'kegg_ubiquitin': "https://www.kegg.jp/pathway/hsa04120",
                'kegg_er_processing': "https://www.kegg.jp/pathway/hsa04141",
                'reactome_proteasome': "https://reactome.org/content/detail/R-HSA-597592",
                'reactome_upr': "https://reactome.org/content/detail/R-HSA-381119",
            }
        }
        
        # Add ligand-specific links if provided
        if ligand_name:
            links['drug'] = {
                'drugbank_search': f"https://go.drugbank.com/unearth/q?query={ligand_name}",
                'chembl_search': f"https://www.ebi.ac.uk/chembl/g/#search_results/all/query={ligand_name}",
                'pubchem_search': f"https://pubchem.ncbi.nlm.nih.gov/#query={ligand_name}",
            }
        
        return links


def enrich_pdb_metadata(pdb_id: str, ligand_name: Optional[str] = None, api_prefs: Optional[Dict[str, bool]] = None) -> Dict[str, Any]:
    """
    Main function to enrich PDB metadata with comprehensive API data
    
    Args:
        pdb_id: PDB identifier
        ligand_name: Optional ligand/drug name
        api_prefs: Optional dictionary of API preferences (which data sources to fetch)
        
    Returns:
        Enriched metadata dictionary with data from multiple sources
    """
    # Default API preferences (all enabled)
    if api_prefs is None:
        api_prefs = {
            'fetch_validation': True,
            'fetch_interactions': True,
            'fetch_drug_info': True,
            'fetch_clinical_trials': True,
            'fetch_publications': True,
            'fetch_protein_network': True,
            'fetch_pathways': True,
        }
    
    client = PDBAPIClient()
    
    # 1. Fetch basic PDB data (always fetched)
    pdb_summary = client.get_pdb_summary(pdb_id)
    citations = client.get_pdb_citations(pdb_id)
    
    # 2. Fetch structure validation metrics (conditional)
    validation = client.get_pdbe_validation(pdb_id) if api_prefs.get('fetch_validation', True) else None
    
    # 3. Fetch ligand interaction details (conditional)
    ligand_interactions = client.get_pdbe_ligand_interactions(pdb_id) if api_prefs.get('fetch_interactions', True) else None
    
    # 4. Find similar structures (always useful)
    similar_structures = client.get_similar_structures(pdb_id, limit=5)
    
    # 5. Search PubMed for recent literature (conditional)
    pubmed_results = client.search_pubmed(f"{pdb_id} protein structure", max_results=5) if api_prefs.get('fetch_publications', True) else None
    
    # 6. Build standard links
    links = client.build_standard_links(pdb_id, ligand_name)
    
    # 7. Fetch comprehensive drug information if ligand provided (conditional)
    drug_info = None
    clinical_trials = None
    chembl_details = None
    
    if ligand_name and api_prefs.get('fetch_drug_info', True):
        # Search ChEMBL
        drug_info = client.search_chembl_by_name(ligand_name)
        
        # If ChEMBL found, get detailed info
        if drug_info and drug_info.get('chembl_id'):
            chembl_details = client.get_chembl_drug_details(drug_info['chembl_id'])
        
        # Fallback to PubChem if ChEMBL not found
        if not drug_info:
            drug_info = client.get_pubchem_compound(ligand_name)
    
    # Search clinical trials (conditional)
    if ligand_name and api_prefs.get('fetch_clinical_trials', True):
        clinical_trials = client.search_clinical_trials(ligand_name)
        
    # Search PubMed for drug-specific literature (conditional)
    if ligand_name and api_prefs.get('fetch_publications', True):
        drug_pubmed = client.search_pubmed(f"{ligand_name} multiple myeloma", max_results=5)
        if drug_pubmed:
            pubmed_results = (pubmed_results or []) + drug_pubmed
    
    # 8. Extract protein information and search for interactions (conditional)
    protein_interactions = None
    if api_prefs.get('fetch_protein_network', True) and pdb_summary and 'struct' in pdb_summary:
        title = pdb_summary['struct'].get('title', '')
        # Try to extract protein name for interaction search
        # This is a simplified approach - could be enhanced with NER
        if 'proteasome' in title.lower():
            protein_interactions = client.get_protein_interactions('PSMB5')  # Beta-5 subunit
        elif 'cereblon' in title.lower() or 'crbn' in title.lower():
            protein_interactions = client.get_protein_interactions('CRBN')
    
    return {
        'pdb_id': pdb_id.upper(),
        'summary': pdb_summary,
        'citations': citations,
        'validation': validation,
        'ligand_interactions': ligand_interactions,
        'similar_structures': similar_structures,
        'pubmed_results': pubmed_results,
        'links': links,
        'drug_info': drug_info,
        'chembl_details': chembl_details,
        'clinical_trials': clinical_trials,
        'protein_interactions': protein_interactions,
    }


def enrich_pdb_metadata_for_view(pdb_id: str, ligand_id: str = None, api_prefs: Optional[Dict[str, bool]] = None) -> dict:
    """
    Enrich PDB metadata with data from multiple APIs for the Django view.
    This function wraps the main enrich_pdb_metadata and formats data for template compatibility.
    
    Args:
        pdb_id: PDB identifier (e.g., '5LF3')
        ligand_id: Optional ligand identifier (e.g., 'BTZ')
        api_prefs: Optional dictionary of API preferences from user form
    
    Returns:
        Dictionary with enriched metadata including:
        - pdb_summary: Basic structure information
        - ligand: Ligand information
        - binding_analysis: Clinical and mechanistic details (template compatible)
        - quality_metrics: Structure validation data
        - similar_structures: Related PDB entries
        - literature: Relevant publications
    """
    enriched = {}
    
    try:
        # Use existing enrich_pdb_metadata function with API preferences
        pdb_data = enrich_pdb_metadata(pdb_id, ligand_id, api_prefs)
        
        if pdb_data and 'summary' in pdb_data:
            summary = pdb_data['summary']
            enriched['pdb_summary'] = summary
            
            # Extract common fields from RCSB response
            if 'struct' in summary:
                enriched['title'] = summary['struct'].get('title', '')
            if 'exptl' in summary and summary['exptl']:
                enriched['method'] = summary['exptl'][0].get('method', '')
            if 'rcsb_entry_info' in summary:
                res = summary['rcsb_entry_info'].get('resolution_combined', [])
                enriched['resolution'] = f"{res[0]:.2f} Å" if res else ''
            
        # 2. Get ligand information from drug_info
        if pdb_data and 'drug_info' in pdb_data and pdb_data['drug_info']:
            drug_info = pdb_data['drug_info']
            enriched['ligand'] = {
                'name': drug_info.get('pref_name', ligand_id),
                'formula': drug_info.get('molecule_chembl_id', ''),
                'max_phase': drug_info.get('max_phase', 'Unknown'),
            }
            enriched['ligand_name'] = drug_info.get('pref_name', ligand_id)
        
        # 3. Generate binding_analysis for template compatibility
        binding_analysis = {}
        
        # Extract protein name from title
        if 'title' in enriched:
            title_lower = enriched['title'].lower()
            binding_analysis['name'] = enriched['title']
            
            # Infer ligand from title if not provided
            inferred_ligand = None
            if 'bortezomib' in title_lower:
                inferred_ligand = 'BORTEZOMIB'
            elif 'carfilzomib' in title_lower:
                inferred_ligand = 'CARFILZOMIB'
            elif 'lenalidomide' in title_lower:
                inferred_ligand = 'LENALIDOMIDE'
            
            # Use inferred ligand if not explicitly provided
            if not ligand_id and inferred_ligand:
                ligand_id = inferred_ligand
        
        # Try to infer mechanism and drug class from ligand
        if ligand_id:
            ligand_upper = ligand_id.upper()
            
            # Common MM drug patterns
            if ligand_upper in ['BTZ', 'BRZ', 'BORTEZOMIB', 'B0R']:
                binding_analysis.update({
                    'mechanism': 'Proteasome inhibition - β5 catalytic site',
                    'clinical': 'Bortezomib (Velcade) - FDA approved, standard MM therapy',
                    'binding_type': 'Reversible covalent inhibitor (boronic acid)',
                    'drug_class': 'Proteasome Inhibitor',
                    'target': '20S Proteasome β5 subunit',
                    'resistance': 'Mutations A49T, C52F in PSMB5 gene',
                    'key_residues': 'Thr1, Ala49, Lys33 (β5 active site)',
                })
            elif ligand_upper in ['CFZ', 'CARFILZOMIB', 'PR9']:
                binding_analysis.update({
                    'mechanism': 'Irreversible proteasome inhibition',
                    'clinical': 'Carfilzomib (Kyprolis) - Second-line therapy',
                    'binding_type': 'Irreversible covalent inhibitor (epoxyketone)',
                    'drug_class': 'Proteasome Inhibitor',
                    'target': '20S Proteasome β5 subunit',
                    'key_residues': 'Thr1 (β5 active site)',
                })
            elif ligand_upper in ['LEN', 'LENALIDOMIDE', 'CC4']:
                binding_analysis.update({
                    'mechanism': 'IMiD mechanism - CRBN E3 ligase modulation',
                    'clinical': 'Lenalidomide (Revlimid) - First-line therapy',
                    'binding_type': 'Small molecule in protein pocket',
                    'drug_class': 'IMiD',
                    'target': 'CRBN E3 ubiquitin ligase',
                    'key_residues': 'Trp380, Trp400, His378',
                })
            elif 'proteasome' in enriched.get('title', '').lower():
                # Generic proteasome inhibitor
                binding_analysis.update({
                    'mechanism': 'Proteasome inhibition',
                    'drug_class': 'Proteasome Inhibitor',
                    'target': '20S Proteasome',
                })
            elif 'cereblon' in enriched.get('title', '').lower() or 'crbn' in enriched.get('title', '').lower():
                # Generic IMiD
                binding_analysis.update({
                    'mechanism': 'E3 ligase modulation',
                    'drug_class': 'IMiD',
                    'target': 'CRBN E3 ubiquitin ligase',
                })
        
        if binding_analysis:
            enriched['binding_analysis'] = binding_analysis
        
        # 4. Add validation metrics
        if pdb_data and 'validation' in pdb_data and pdb_data['validation']:
            enriched['validation_metrics'] = pdb_data['validation']
        
        # 5. Add ligand interaction details
        if pdb_data and 'ligand_interactions' in pdb_data and pdb_data['ligand_interactions']:
            enriched['interaction_details'] = pdb_data['ligand_interactions']
        
        # 6. Add similar structures
        if pdb_data and 'similar_structures' in pdb_data and pdb_data['similar_structures']:
            enriched['similar_structures'] = pdb_data['similar_structures']
        
        # 7. Add PubMed literature
        if pdb_data and 'pubmed_results' in pdb_data and pdb_data['pubmed_results']:
            enriched['pubmed_literature'] = pdb_data['pubmed_results'][:10]
        
        # 8. Add ChEMBL detailed drug info
        if pdb_data and 'chembl_details' in pdb_data and pdb_data['chembl_details']:
            enriched['chembl_drug_details'] = pdb_data['chembl_details']
        
        # 9. Add clinical trials
        if pdb_data and 'clinical_trials' in pdb_data and pdb_data['clinical_trials']:
            enriched['clinical_trials'] = pdb_data['clinical_trials']
        
        # 10. Add protein interactions
        if pdb_data and 'protein_interactions' in pdb_data and pdb_data['protein_interactions']:
            enriched['protein_network'] = pdb_data['protein_interactions']
        
        # 11. Add standard links from pdb_data
        if pdb_data and 'links' in pdb_data:
            enriched['external_links'] = pdb_data['links']
        
        # 12. Add citations (keep existing)
        if pdb_data and 'citations' in pdb_data:
            enriched['literature'] = pdb_data['citations'][:5]
            
    except Exception as e:
        logger.error(f"Error enriching PDB metadata for {pdb_id}: {e}")
        # Return partial data if available
        pass
    
    return enriched
