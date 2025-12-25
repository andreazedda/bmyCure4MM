#!/usr/bin/env python
"""
Test script to verify async loading endpoint works correctly.
"""

import os
import django

os.environ.setdefault('DJANGO_SECRET_KEY', 'test-secret-key')
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'mmportal.settings')
django.setup()

from django.test import RequestFactory
from django.contrib.auth import get_user_model
from chemtools import models, views

def test_async_endpoint():
    """Test the new job_enriched_data endpoint"""
    
    print("=" * 80)
    print("Testing Async Loading Endpoint")
    print("=" * 80)
    
    # Create a test user
    User = get_user_model()
    user, _ = User.objects.get_or_create(username='testuser', defaults={'email': 'test@test.com'})
    
    # Create a test binding job
    job, created = models.ChemJob.objects.get_or_create(
        user=user,
        kind=models.ChemJob.BIND,
        input_a='5T3H',  # PDB ID
        input_b='POM',   # Pomalidomide
        defaults={'status': 'done'}
    )
    
    print(f"\n✅ Test job created: ID={job.pk}, PDB={job.input_a}, Ligand={job.input_b}")
    
    # Create a mock request
    factory = RequestFactory()
    request = factory.get(f'/chemtools/job/{job.pk}/enriched-data.json')
    request.user = user
    
    print(f"\n🌐 Calling job_enriched_data(request, pk={job.pk})...")
    print("   This will fetch data from ChEMBL, PubMed, ClinicalTrials.gov...")
    print("   Expected time: 3-5 seconds\n")
    
    # Call the view
    import time
    start = time.time()
    
    try:
        response = views.job_enriched_data(request, job.pk)
        elapsed = time.time() - start
        
        print(f"✅ Response received in {elapsed:.2f}s")
        print(f"   Status code: {response.status_code}")
        
        if response.status_code == 200:
            import json
            data = json.loads(response.content)
            
            print(f"\n📊 Data keys returned: {list(data.keys())}")
            
            # Check MM Efficacy
            if data.get('mm_efficacy_profile'):
                print(f"\n✅ MM Efficacy Profile:")
                print(f"   Score: {data['mm_efficacy_profile'].get('overall_score')}/100")
                print(f"   Confidence: {data['mm_efficacy_profile'].get('confidence')}")
            else:
                print(f"\n⚠️  No MM Efficacy data")
            
            # Check Survival Impact
            if data.get('survival_impact'):
                print(f"\n✅ Survival Impact:")
                print(f"   PFS: {data['survival_impact'].get('median_pfs_months')} months")
                print(f"   OS: {data['survival_impact'].get('median_os_months')} months")
            else:
                print(f"\n⚠️  No Survival Impact data")
            
            # Check Toxicity Profile
            if data.get('toxicity_profile'):
                print(f"\n✅ Toxicity Profile:")
                print(f"   Risk: {data['toxicity_profile'].get('overall_risk')}")
                print(f"   Score: {data['toxicity_profile'].get('risk_score')}/100")
            else:
                print(f"\n⚠️  No Toxicity data")
            
            print("\n" + "=" * 80)
            print("✅ ASYNC ENDPOINT WORKING! JavaScript will load this data.")
            print("=" * 80)
            
        else:
            print(f"\n❌ ERROR: Status code {response.status_code}")
            
    except Exception as e:
        print(f"\n❌ ERROR: {e}")
        import traceback
        traceback.print_exc()

if __name__ == '__main__':
    test_async_endpoint()
