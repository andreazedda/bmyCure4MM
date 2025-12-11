#!/usr/bin/env python
"""
Quick test of the PDB API integration
"""
import os
import django

os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'mmportal.settings')
django.setup()

from chemtools.pdb_api_client import enrich_pdb_metadata_for_view

def test_pdb_api():
    """Test API enrichment for known PDB structures"""
    
    test_cases = [
        ('5LF3', None, 'Bortezomib + 20S Proteasome'),
        ('4KW5', 'CC4', 'Lenalidomide + Cereblon'),
        ('2F16', None, 'Bortezomib + Proteasome (older)'),
    ]
    
    for pdb_id, ligand_id, description in test_cases:
        print(f"\n{'='*60}")
        print(f"Testing: {pdb_id} ({description})")
        print(f"Ligand: {ligand_id or 'None specified'}")
        print('='*60)
        
        try:
            result = enrich_pdb_metadata_for_view(pdb_id, ligand_id)
            
            print(f"\n✅ API call successful for {pdb_id}")
            print(f"   Keys returned: {list(result.keys())}")
            
            if 'title' in result:
                print(f"   Title: {result['title'][:80]}...")
            
            if 'method' in result:
                print(f"   Method: {result['method']}")
            
            if 'resolution' in result:
                print(f"   Resolution: {result['resolution']}")
            
            if 'binding_analysis' in result:
                ba = result['binding_analysis']
                print(f"\n   Binding Analysis:")
                print(f"     - Drug Class: {ba.get('drug_class', 'N/A')}")
                print(f"     - Target: {ba.get('target', 'N/A')}")
                print(f"     - Mechanism: {ba.get('mechanism', 'N/A')[:60]}...")
                print(f"     - Clinical: {ba.get('clinical', 'N/A')[:60]}...")
            
            if 'ligand' in result:
                ligand = result['ligand']
                print(f"\n   Ligand Info:")
                print(f"     - Name: {ligand.get('name', 'N/A')}")
                if 'max_phase' in ligand:
                    print(f"     - Max Phase: {ligand['max_phase']}")
            
            if 'external_links' in result:
                print(f"\n   External links generated: {len(result['external_links'])} categories")
            
        except Exception as e:
            print(f"\n❌ Error processing {pdb_id}: {e}")
            import traceback
            traceback.print_exc()

if __name__ == '__main__':
    print("Testing PDB API Integration")
    print("This will make real API calls to RCSB PDB, ChEMBL, etc.")
    print("Please wait, APIs may be slow...\n")
    test_pdb_api()
    print("\n" + "="*60)
    print("Test completed!")
