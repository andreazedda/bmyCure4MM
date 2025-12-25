"""
Test script to verify dynamic ChEMBL data fetching
"""
import requests
import json

class SimpleChEMBLClient:
    """Simplified ChEMBL client for testing"""
    CHEMBL_API_BASE = "https://www.ebi.ac.uk/chembl/api/data"
    
    def __init__(self):
        self.session = requests.Session()
        self.timeout = 10
    
    def search_chembl_by_name(self, drug_name):
        try:
            url = f"{self.CHEMBL_API_BASE}/molecule.json"
            params = {'pref_name__iexact': drug_name, 'limit': 1}
            response = self.session.get(url, params=params, timeout=self.timeout)
            response.raise_for_status()
            data = response.json()
            if data.get('molecules'):
                mol = data['molecules'][0]
                return {'chembl_id': mol.get('molecule_chembl_id')}
            return None
        except Exception as e:
            print(f"Error: {e}")
            return None
    
    def get_chembl_drug_details(self, chembl_id):
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
                'atc_classifications': data.get('atc_classifications', []),
            }
            
            # Fetch mechanisms
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
                    } for m in mech_data['mechanisms']]
                    
                    # Fetch target details
                    targets = []
                    for mech in result['mechanisms'][:3]:
                        target_id = mech.get('target_chembl_id')
                        if target_id:
                            target_detail = self._get_target_details(target_id)
                            if target_detail:
                                targets.append(target_detail)
                    if targets:
                        result['targets'] = targets
            
            return result
        except Exception as e:
            print(f"Error: {e}")
            return None
    
    def _get_target_details(self, target_chembl_id):
        try:
            url = f"{self.CHEMBL_API_BASE}/target/{target_chembl_id}.json"
            response = self.session.get(url, timeout=self.timeout)
            response.raise_for_status()
            data = response.json()
            
            components = []
            for comp in data.get('target_components', []):
                components.append({
                    'gene_symbol': comp.get('component_synonym'),
                    'accession': comp.get('accession'),
                    'protein_name': comp.get('component_description'),
                })
            
            return {
                'target_chembl_id': target_chembl_id,
                'pref_name': data.get('pref_name'),
                'target_components': components,
            }
        except Exception as e:
            return None

def test_dynamic_drug_fetching():
    """Test dynamic drug mechanism and target fetching"""
    client = SimpleChEMBLClient()
    
    # Test 1: MM drug (Bortezomib)
    print("=" * 60)
    print("Test 1: Bortezomib (MM drug)")
    print("=" * 60)
    
    drug_info = client.search_chembl_by_name('Bortezomib')
    if drug_info:
        print(f"✓ Found ChEMBL ID: {drug_info.get('chembl_id')}")
        
        details = client.get_chembl_drug_details(drug_info['chembl_id'])
        if details:
            print(f"✓ Drug Name: {details.get('name')}")
            print(f"✓ Max Phase: {details.get('max_phase')}")
            print(f"✓ Drug Type: {details.get('type')}")
            print(f"✓ ATC Classifications: {details.get('atc_classifications', [])[:2]}")
            
            if details.get('mechanisms'):
                print(f"\n✓ Mechanisms of Action:")
                for i, mech in enumerate(details['mechanisms'][:2], 1):
                    print(f"  {i}. {mech.get('mechanism', 'N/A')}")
                    print(f"     Action: {mech.get('action_type', 'N/A')}")
                    print(f"     Target: {mech.get('target_chembl_id', 'N/A')}")
            
            if details.get('targets'):
                print(f"\n✓ Targets:")
                for i, target in enumerate(details['targets'][:2], 1):
                    print(f"  {i}. {target.get('pref_name', 'N/A')}")
                    if target.get('target_components'):
                        comp = target['target_components'][0]
                        print(f"     Gene: {comp.get('gene_symbol', 'N/A')}")
                        print(f"     UniProt: {comp.get('accession', 'N/A')}")
        else:
            print("✗ Failed to fetch drug details")
    else:
        print("✗ Drug not found in ChEMBL")
    
    # Test 2: Kinase inhibitor (Imatinib)
    print("\n" + "=" * 60)
    print("Test 2: Imatinib (Kinase inhibitor)")
    print("=" * 60)
    
    drug_info = client.search_chembl_by_name('Imatinib')
    if drug_info:
        print(f"✓ Found ChEMBL ID: {drug_info.get('chembl_id')}")
        
        details = client.get_chembl_drug_details(drug_info['chembl_id'])
        if details:
            print(f"✓ Drug Name: {details.get('name')}")
            print(f"✓ Max Phase: {details.get('max_phase')}")
            print(f"✓ ATC Classifications: {details.get('atc_classifications', [])[:2]}")
            
            if details.get('mechanisms'):
                print(f"\n✓ Mechanisms ({len(details['mechanisms'])} found):")
                for i, mech in enumerate(details['mechanisms'][:3], 1):
                    print(f"  {i}. {mech.get('mechanism', 'N/A')}")
            
            if details.get('targets'):
                print(f"\n✓ Targets ({len(details['targets'])} found):")
                for i, target in enumerate(details['targets'][:3], 1):
                    print(f"  {i}. {target.get('pref_name', 'N/A')}")
        else:
            print("✗ Failed to fetch drug details")
    else:
        print("✗ Drug not found in ChEMBL")
    
    print("\n" + "=" * 60)
    print("✓ Dynamic data fetching test complete!")
    print("=" * 60)

if __name__ == '__main__':
    test_dynamic_drug_fetching()
