#!/usr/bin/env python3
"""
Test script to verify MM drug preset and API preference functionality
"""

import os
import sys
import django

# Setup Django
sys.path.insert(0, '/Volumes/nvme/Github/bmyCure4MM')
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'mmportal.settings')
django.setup()

from chemtools.forms import BindingVizForm
from chemtools.models import ChemJob
from django.contrib.auth.models import User

def test_form_validation():
    """Test form with various inputs"""
    print("=" * 60)
    print("TEST 1: Form Validation")
    print("=" * 60)
    
    # Test with valid data
    data = {
        'mm_drug_preset': '5MX5_CFZ',
        'pdb_id': '5MX5',
        'ligand': 'CFZ',
        'fetch_validation': True,
        'fetch_interactions': True,
        'fetch_drug_info': True,
        'fetch_clinical_trials': True,
        'fetch_publications': True,
        'fetch_protein_network': True,
        'fetch_pathways': True,
        'mm_focus_resistance': True,
        'mm_focus_combinations': False,
        'mm_focus_toxicity': True,
    }
    
    form = BindingVizForm(data)
    if form.is_valid():
        print("✓ Form validation passed")
        print(f"  PDB ID: {form.cleaned_data['pdb_id']}")
        print(f"  Ligand: {form.cleaned_data['ligand']}")
        print(f"  Preset: {form.cleaned_data['mm_drug_preset']}")
        print(f"  Fetch validation: {form.cleaned_data['fetch_validation']}")
        print(f"  MM focus resistance: {form.cleaned_data['mm_focus_resistance']}")
    else:
        print("✗ Form validation failed:")
        for field, errors in form.errors.items():
            print(f"  {field}: {errors}")
    print()

def test_mm_drug_presets():
    """Test MM drug preset values"""
    print("=" * 60)
    print("TEST 2: MM Drug Preset Values")
    print("=" * 60)
    
    presets = {
        '5LF3_BOR': ('5LF3', 'BOR', 'Bortezomib'),
        '4KW5_LEN': ('4KW5', 'LEN', 'Lenalidomide'),
        '5MX5_CFZ': ('5MX5', 'CFZ', 'Carfilzomib'),
        '5T3H_POM': ('5T3H', 'POM', 'Pomalidomide'),
        '5M2B_IXA': ('5M2B', 'IXA', 'Ixazomib'),
    }
    
    for preset_key, (expected_pdb, expected_ligand, drug_name) in presets.items():
        data = {
            'mm_drug_preset': preset_key,
            'pdb_id': expected_pdb,
            'ligand': expected_ligand,
        }
        form = BindingVizForm(data)
        if form.is_valid():
            pdb = form.cleaned_data['pdb_id']
            lig = form.cleaned_data['ligand']
            if pdb == expected_pdb and lig == expected_ligand:
                print(f"✓ {drug_name}: PDB={pdb}, Ligand={lig}")
            else:
                print(f"✗ {drug_name}: Expected PDB={expected_pdb}, got {pdb}; Expected Ligand={expected_ligand}, got {lig}")
        else:
            print(f"✗ {drug_name}: Form validation failed - {form.errors}")
    print()

def test_api_preferences_storage():
    """Test that API preferences are stored correctly"""
    print("=" * 60)
    print("TEST 3: API Preferences Storage")
    print("=" * 60)
    
    # Get or create test user
    user, _ = User.objects.get_or_create(username='test_user')
    
    api_prefs = {
        'fetch_validation': True,
        'fetch_interactions': False,
        'fetch_drug_info': True,
        'fetch_clinical_trials': False,
        'fetch_publications': True,
        'fetch_protein_network': False,
        'fetch_pathways': True,
        'mm_focus_resistance': True,
        'mm_focus_combinations': True,
        'mm_focus_toxicity': False,
    }
    
    # Create test job
    job = ChemJob.objects.create(
        kind=ChemJob.BIND,
        input_a='5MX5',
        input_b='CFZ',
        api_preferences=api_prefs,
        user=user
    )
    
    # Retrieve and verify
    retrieved_job = ChemJob.objects.get(pk=job.pk)
    stored_prefs = retrieved_job.api_preferences
    
    all_match = True
    for key, expected_value in api_prefs.items():
        actual_value = stored_prefs.get(key)
        if actual_value == expected_value:
            print(f"✓ {key}: {actual_value}")
        else:
            print(f"✗ {key}: Expected {expected_value}, got {actual_value}")
            all_match = False
    
    if all_match:
        print("\n✓ All API preferences stored and retrieved correctly")
    else:
        print("\n✗ Some API preferences did not match")
    
    # Cleanup
    job.delete()
    print()

def test_pdb_id_validation():
    """Test PDB ID validation"""
    print("=" * 60)
    print("TEST 4: PDB ID Validation")
    print("=" * 60)
    
    test_cases = [
        ('5MX5', True, 'Valid PDB ID'),
        ('4KW5', True, 'Valid PDB ID'),
        ('1ABC', True, 'Valid PDB ID'),
        ('INVALID', False, 'Invalid format'),
        ('12345', False, 'Too long'),
        ('', False, 'Empty'),
    ]
    
    for pdb_id, should_be_valid, description in test_cases:
        data = {
            'pdb_id': pdb_id,
            'ligand': 'LIG',
        }
        form = BindingVizForm(data)
        is_valid = form.is_valid()
        
        if is_valid == should_be_valid:
            status = "✓"
        else:
            status = "✗"
        
        print(f"{status} {description}: '{pdb_id}' - Valid={is_valid}")
    print()

def test_checkbox_defaults():
    """Test that checkboxes have correct default values"""
    print("=" * 60)
    print("TEST 5: Checkbox Default Values")
    print("=" * 60)
    
    # Submit form with minimal data (no checkboxes)
    data = {
        'pdb_id': '5MX5',
        'ligand': 'CFZ',
    }
    
    form = BindingVizForm(data)
    if form.is_valid():
        # API sources should default to True (checked)
        api_fields = [
            'fetch_validation',
            'fetch_interactions',
            'fetch_drug_info',
            'fetch_clinical_trials',
            'fetch_publications',
            'fetch_protein_network',
            'fetch_pathways',
        ]
        
        # MM focus should default to False (unchecked)
        mm_fields = [
            'mm_focus_resistance',
            'mm_focus_combinations',
            'mm_focus_toxicity',
        ]
        
        all_correct = True
        for field in api_fields:
            value = form.cleaned_data.get(field, None)
            expected = True
            if value == expected:
                print(f"✓ {field}: {value} (expected {expected})")
            else:
                print(f"✗ {field}: {value} (expected {expected})")
                all_correct = False
        
        for field in mm_fields:
            value = form.cleaned_data.get(field, None)
            expected = False
            if value == expected:
                print(f"✓ {field}: {value} (expected {expected})")
            else:
                print(f"✗ {field}: {value} (expected {expected})")
                all_correct = False
        
        if all_correct:
            print("\n✓ All checkbox defaults are correct")
        else:
            print("\n✗ Some checkbox defaults are incorrect")
    else:
        print("✗ Form validation failed")
        print(form.errors)
    print()

def main():
    """Run all tests"""
    print("\n" + "=" * 60)
    print("BMYCURE4MM - API PREFERENCES & MM PRESETS TEST SUITE")
    print("=" * 60 + "\n")
    
    try:
        test_form_validation()
        test_mm_drug_presets()
        test_api_preferences_storage()
        test_pdb_id_validation()
        test_checkbox_defaults()
        
        print("=" * 60)
        print("ALL TESTS COMPLETED")
        print("=" * 60)
        
    except Exception as e:
        print(f"\n✗ Error during testing: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0

if __name__ == '__main__':
    sys.exit(main())
