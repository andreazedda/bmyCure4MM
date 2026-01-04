#!/usr/bin/env python
"""
Quick test to verify MM efficacy, survival impact, and toxicity features are working.
Tests with PDB 5T3H + ligand POM (Pomalidomide).
"""

import os
import django
import logging

# Setup Django
os.environ.setdefault('DJANGO_SECRET_KEY', 'test-secret-key-for-testing')
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'mmportal.settings')
django.setup()

# Now import Django modules
from chemtools.pdb_api_client import enrich_pdb_metadata_for_view

# Enable detailed logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def test_mm_features():
    """Test MM efficacy, survival, and toxicity estimation for Pomalidomide"""
    
    print("=" * 80)
    print("Testing MM Features with Pomalidomide (PDB: 5T3H, Ligand: POM)")
    print("=" * 80)
    
    # Test with default API preferences (should enable all features)
    print("\n1. Fetching enriched metadata with default preferences...")
    enriched = enrich_pdb_metadata_for_view('5T3H', ligand_id='POM', api_prefs=None)
    
    if not enriched:
        print("❌ ERROR: No enriched data returned!")
        return False
    
    print(f"✅ Enriched data keys: {list(enriched.keys())}")
    
    # Check for MM efficacy
    print("\n2. Checking MM Efficacy Profile...")
    if 'mm_efficacy_profile' in enriched:
        efficacy = enriched['mm_efficacy_profile']
        print(f"✅ MM Efficacy Profile Found!")
        print(f"   Overall Score: {efficacy.get('overall_score')}/100")
        print(f"   Confidence: {efficacy.get('confidence')}")
        print(f"   Target Relevance: {efficacy.get('target_relevance_score')}")
        print(f"   Mechanism Relevance: {efficacy.get('mechanism_relevance_score')}")
        print(f"   Clinical Evidence: {efficacy.get('clinical_evidence_score')}")
        print(f"   Literature Evidence: {efficacy.get('literature_evidence_score')}")
    else:
        print("❌ MM Efficacy Profile NOT FOUND")
        return False
    
    # Check for Survival Impact
    print("\n3. Checking Survival Impact...")
    if 'survival_impact' in enriched:
        survival = enriched['survival_impact']
        print(f"✅ Survival Impact Found!")
        print(f"   Median PFS: {survival.get('median_pfs_months')} months")
        print(f"   Median OS: {survival.get('median_os_months')} months")
        print(f"   Response Rate: {survival.get('response_rate_percent')}%")
        print(f"   HR PFS: {survival.get('hazard_ratio_pfs')}")
        print(f"   HR OS: {survival.get('hazard_ratio_os')}")
        print(f"   Confidence: {survival.get('confidence_level')}")
    else:
        print("❌ Survival Impact NOT FOUND")
        return False
    
    # Check for Toxicity Profile
    print("\n4. Checking Toxicity Profile...")
    if 'toxicity_profile' in enriched:
        toxicity = enriched['toxicity_profile']
        print(f"✅ Toxicity Profile Found!")
        print(f"   Overall Risk: {toxicity.get('overall_risk')}")
        print(f"   Risk Score: {toxicity.get('risk_score')}/100")
        print(f"   Risk-Benefit Ratio: {toxicity.get('risk_benefit_ratio')}")
        print(f"   Common AEs: {len(toxicity.get('common_adverse_events', []))} listed")
        print(f"   Serious AEs: {len(toxicity.get('serious_adverse_events', []))} listed")
        print(f"   Black Box Warnings: {len(toxicity.get('black_box_warnings', []))} listed")
        print(f"   Management Strategies: {len(toxicity.get('management_strategies', []))} listed")
        
        # Show black box warnings for Pomalidomide (IMiD)
        if toxicity.get('black_box_warnings'):
            print("\n   ⚠️  Black Box Warnings:")
            for warning in toxicity['black_box_warnings']:
                print(f"      - {warning}")
    else:
        print("❌ Toxicity Profile NOT FOUND")
        return False
    
    print("\n" + "=" * 80)
    print("✅ ALL FEATURES WORKING! Data will display in job_detail.html")
    print("=" * 80)
    return True

if __name__ == '__main__':
    success = test_mm_features()
    exit(0 if success else 1)
