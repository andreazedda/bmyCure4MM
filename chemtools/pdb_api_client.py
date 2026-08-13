"""
PDB and Biological Database API Client
Provides standardized access to structural and biological databases
"""
import requests
from typing import Dict, List, Optional, Any
import logging

from mmportal.governance import EpistemicLabel, governance_metadata

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
        """Fetch detailed drug information from ChEMBL including mechanisms and targets"""
        try:
            url = f"{self.CHEMBL_API_BASE}/molecule/{chembl_id}.json"
            response = self.session.get(url, timeout=self.timeout)
            response.raise_for_status()
            data = response.json()
            
            result = {
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
            
            # Fetch mechanism of action data
            mech_url = f"{self.CHEMBL_API_BASE}/mechanism.json"
            mech_params = {'molecule_chembl_id': chembl_id}
            mech_response = self.session.get(mech_url, params=mech_params, timeout=self.timeout)
            if mech_response.status_code == 200:
                mech_data = mech_response.json()
                if mech_data.get('mechanisms'):
                    result['mechanisms'] = [{
                        'action_type': m.get('action_type'),
                        'mechanism': m.get('mechanism_of_action'),
                        'target_chembl_id': m.get('target_chembl_id'),
                        'disease_efficacy': m.get('disease_efficacy'),
                    } for m in mech_data['mechanisms']]
                    
                    # Fetch target details for gene symbols and UniProt IDs
                    targets = []
                    for mech in result['mechanisms']:
                        target_id = mech.get('target_chembl_id')
                        if target_id:
                            target_detail = self._get_chembl_target_details(target_id)
                            if target_detail:
                                targets.append(target_detail)
                    
                    if targets:
                        result['targets'] = targets
            
            return result
        except Exception as e:
            logger.warning(f"Failed to fetch ChEMBL details for {chembl_id}: {e}")
            return None
    
    def _get_chembl_target_details(self, target_chembl_id: str) -> Optional[Dict[str, Any]]:
        """Fetch target details including gene symbols and UniProt IDs"""
        try:
            url = f"{self.CHEMBL_API_BASE}/target/{target_chembl_id}.json"
            response = self.session.get(url, timeout=self.timeout)
            response.raise_for_status()
            data = response.json()
            
            components = []
            for comp in data.get('target_components', []):
                components.append({
                    'gene_symbol': comp.get('component_synonym'),
                    'accession': comp.get('accession'),  # UniProt ID
                    'protein_name': comp.get('component_description'),
                })
            
            return {
                'target_chembl_id': target_chembl_id,
                'target_type': data.get('target_type'),
                'organism': data.get('organism'),
                'pref_name': data.get('pref_name'),
                'target_components': components,
            }
        except Exception as e:
            logger.debug(f"Failed to fetch target details for {target_chembl_id}: {e}")
            return None
    
    def estimate_mm_efficacy(
        self,
        drug_name: str,
        chembl_details: Optional[Dict],
        clinical_trials: Optional[List],
        pubmed_results: Optional[List],
    ) -> Dict[str, Any]:
        """Return an explicit gate instead of the unvalidated legacy score."""
        del drug_name, chembl_details, clinical_trials, pubmed_results
        return {
            "status": "INSUFFICIENT_EVIDENCE",
            "overall_score": None,
            "confidence": "UNKNOWN",
            "interpretation": "No efficacy estimate emitted.",
            "evidence_sources": [],
            "mm_relevance": {},
            "clinical_status": None,
            "key_findings": [],
            "limitations": [
                "The historical relevance score is not a validated efficacy endpoint.",
                "Source-level evidence binding is required before this output can be restored.",
            ],
            "governance": governance_metadata(
                epistemic_label=EpistemicLabel.UNKNOWN,
                output_kind="efficacy_validation_gate",
            ),
        }

    def estimate_survival_impact(
        self,
        efficacy_profile: Dict[str, Any],
        chembl_details: Optional[Dict],
    ) -> Dict[str, Any]:
        """Return an explicit gate because no validated survival model exists."""
        del efficacy_profile, chembl_details
        return {
            "status": "PATIENT_SPECIFIC_PREDICTION_NOT_VALIDATED",
            "median_pfs_months": None,
            "median_os_months": None,
            "response_rate_percent": None,
            "hazard_ratio_pfs": None,
            "hazard_ratio_os": None,
            "survival_benefit": None,
            "model_basis": [],
            "confidence": "UNKNOWN",
            "governance": governance_metadata(
                epistemic_label=EpistemicLabel.UNKNOWN,
                output_kind="survival_validation_gate",
            ),
        }

    def estimate_toxicity_profile(
        self,
        chembl_details: Optional[Dict],
        efficacy_profile: Dict[str, Any],
    ) -> Dict[str, Any]:
        """Return an explicit gate because no validated safety model exists."""
        del chembl_details, efficacy_profile
        return {
            "status": "INSUFFICIENT_EVIDENCE",
            "overall_risk": "UNKNOWN",
            "risk_score": None,
            "common_adverse_events": [],
            "serious_adverse_events": [],
            "dose_limiting_toxicities": [],
            "black_box_warnings": [],
            "mechanism_related_toxicities": [],
            "management_strategies": [],
            "risk_benefit_ratio": None,
            "governance": governance_metadata(
                epistemic_label=EpistemicLabel.UNKNOWN,
                output_kind="toxicity_validation_gate",
            ),
        }

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
    # Default API preferences. Unsupported scientific estimators are disabled.
    if api_prefs is None:
        api_prefs = {
            'fetch_validation': True,
            'fetch_interactions': True,
            'fetch_drug_info': True,
            'fetch_clinical_trials': True,
            'fetch_publications': True,
            'fetch_protein_network': True,
            'fetch_pathways': True,
            'estimate_mm_efficacy': False,
            'estimate_survival_impact': False,
            'estimate_toxicity': False,
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
    if api_prefs.get('fetch_protein_network', True):
        # First try to get gene symbol from ChEMBL targets
        gene_symbol = None
        if chembl_details and 'targets' in chembl_details:
            # Use the first target's gene symbol
            for target in chembl_details['targets']:
                components = target.get('target_components', [])
                if components and components[0].get('gene_symbol'):
                    gene_symbol = components[0]['gene_symbol']
                    break
        
        # Fallback: extract from PDB title using common protein mappings
        if not gene_symbol and pdb_summary and 'struct' in pdb_summary:
            title = pdb_summary['struct'].get('title', '')
            protein_gene_map = {
                'proteasome': 'PSMB5',
                'cereblon': 'CRBN',
                'trypsin': 'PRSS1',
                'thrombin': 'F2',
                'kinase': None,  # Too generic, skip
            }
            for protein_name, gene in protein_gene_map.items():
                if protein_name in title.lower() and gene:
                    gene_symbol = gene
                    break
        
        # Fetch interactions if we have a gene symbol
        if gene_symbol:
            protein_interactions = client.get_protein_interactions(gene_symbol)
    
    # 9. Fetch pathway information (conditional) - based on ligand or protein
    pathways = None
    if api_prefs.get('fetch_pathways', True):
        # Try to get pathways from ChEMBL target data if available
        if chembl_details and 'targets' in chembl_details:
            # Get UniProt IDs from ChEMBL targets
            for target in chembl_details['targets'][:3]:  # Limit to first 3 targets
                uniprot_id = target.get('target_components', [{}])[0].get('accession')
                if uniprot_id:
                    pathways = client.get_reactome_pathways(uniprot_id)
                    if pathways:
                        break  # Use first successful pathway fetch
        
        # Fallback: try to extract protein from title
        if not pathways and pdb_summary and 'struct' in pdb_summary:
            title = pdb_summary['struct'].get('title', '')
            # Map common proteins to UniProt IDs for pathway lookup
            protein_uniprot_map = {
                'proteasome': 'P28074',  # PSMB5
                'cereblon': 'Q96SW2',    # CRBN
                'trypsin': 'P07477',     # PRSS1
                'thrombin': 'P00734',    # F2
            }
            for protein_name, uniprot_id in protein_uniprot_map.items():
                if protein_name in title.lower():
                    pathways = client.get_reactome_pathways(uniprot_id)
                    break
    
    # Legacy efficacy, survival, and toxicity estimators are intentionally
    # unavailable. They did not bind outputs to governed evidence and could be
    # misread as clinically validated conclusions.
    mm_efficacy = None
    survival_impact = None
    toxicity_profile = None
    unsupported_governance = governance_metadata(
        epistemic_label=EpistemicLabel.UNKNOWN,
        output_kind='unsupported_chemistry_estimators',
    )
    unsupported_scientific_estimators = {
        'status': 'INSUFFICIENT_EVIDENCE',
        'governance': unsupported_governance,
        'requested': {
            name: bool(api_prefs.get(name, False))
            for name in (
                'estimate_mm_efficacy',
                'estimate_survival_impact',
                'estimate_toxicity',
            )
        },
        'limitations': [
            'No validated efficacy, survival, toxicity, or risk-benefit estimate is emitted.',
            'Source-level evidence and governed validation are required.',
        ],
    }
    
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
        'pathways': pathways,
        'mm_efficacy': mm_efficacy,
        'survival_impact': survival_impact,
        'toxicity_profile': toxicity_profile,
        'unsupported_scientific_estimators': unsupported_scientific_estimators,
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
        
        # 3. Generate binding_analysis dynamically from ChEMBL data
        binding_analysis = {}
        
        # Extract title for context
        if 'title' in enriched:
            binding_analysis['name'] = enriched['title']
        
        # Use ChEMBL data if available
        chembl_details = pdb_data.get('chembl_details') if pdb_data else None
        if chembl_details:
            # Extract mechanism of action from ChEMBL
            if chembl_details.get('mechanisms'):
                primary_mech = chembl_details['mechanisms'][0]
                binding_analysis['mechanism'] = primary_mech.get('mechanism') or primary_mech.get('action_type', 'Unknown mechanism')
                if primary_mech.get('target_name'):
                    binding_analysis['target'] = primary_mech['target_name']
            
            # Extract drug class from ATC classifications
            if chembl_details.get('atc_classifications'):
                atc = chembl_details['atc_classifications'][0]
                binding_analysis['drug_class'] = atc.split(' - ')[-1] if ' - ' in atc else atc
            
            # Build clinical info from approval data
            clinical_info = []
            if chembl_details.get('name'):
                clinical_info.append(chembl_details['name'])
            
            max_phase = chembl_details.get('max_phase')
            if max_phase == 4:
                clinical_info.append('FDA approved')
            elif max_phase == 3:
                clinical_info.append('Phase III trials')
            elif max_phase == 2:
                clinical_info.append('Phase II trials')
            elif max_phase == 1:
                clinical_info.append('Phase I trials')
            
            if chembl_details.get('first_approval'):
                clinical_info.append(f"First approved: {chembl_details['first_approval']}")
            
            if clinical_info:
                binding_analysis['clinical'] = ' - '.join(clinical_info)
            
            # Add administration routes
            routes = []
            if chembl_details.get('oral'):
                routes.append('oral')
            if chembl_details.get('parenteral'):
                routes.append('parenteral')
            if chembl_details.get('topical'):
                routes.append('topical')
            if routes:
                binding_analysis['administration'] = ', '.join(routes)
            
            # Add safety warnings
            if chembl_details.get('black_box_warning'):
                binding_analysis['warning'] = 'Black box warning'
        
        # Fallback: infer basic info from PDB title if ChEMBL data unavailable
        if not binding_analysis.get('mechanism') and 'title' in enriched:
            title_lower = enriched['title'].lower()
            
            # Generic protein-based inference
            if 'proteasome' in title_lower:
                binding_analysis['target'] = '20S Proteasome'
                binding_analysis['drug_class'] = 'Proteasome Inhibitor'
            elif 'cereblon' in title_lower or 'crbn' in title_lower:
                binding_analysis['target'] = 'CRBN E3 ubiquitin ligase'
                binding_analysis['drug_class'] = 'IMiD'
            elif 'kinase' in title_lower:
                binding_analysis['drug_class'] = 'Kinase Inhibitor'
            elif 'antibody' in title_lower or 'fab' in title_lower:
                binding_analysis['drug_class'] = 'Therapeutic Antibody'
        
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
        
        # 11. Add pathways from Reactome
        if pdb_data and 'pathways' in pdb_data and pdb_data['pathways']:
            enriched['reactome_pathways'] = pdb_data['pathways']
        
        # 12. Add MM efficacy estimation
        if pdb_data and 'mm_efficacy' in pdb_data and pdb_data['mm_efficacy']:
            enriched['mm_efficacy_profile'] = pdb_data['mm_efficacy']
        
        # 13. Add survival impact estimation
        if pdb_data and 'survival_impact' in pdb_data and pdb_data['survival_impact']:
            enriched['survival_impact'] = pdb_data['survival_impact']
        
        # 14. Add toxicity profile
        if pdb_data and 'toxicity_profile' in pdb_data and pdb_data['toxicity_profile']:
            enriched['toxicity_profile'] = pdb_data['toxicity_profile']

        if pdb_data and 'unsupported_scientific_estimators' in pdb_data:
            enriched['unsupported_scientific_estimators'] = pdb_data['unsupported_scientific_estimators']
        
        # 15. Add standard links from pdb_data
        if pdb_data and 'links' in pdb_data:
            enriched['external_links'] = pdb_data['links']
        
        # 14. Add citations (keep existing)
        if pdb_data and 'citations' in pdb_data:
            enriched['literature'] = pdb_data['citations'][:5]
            
    except Exception as e:
        logger.error(f"Error enriching PDB metadata for {pdb_id}: {e}")
        # Return partial data if available
        pass
    
    return enriched
