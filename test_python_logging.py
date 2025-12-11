#!/usr/bin/env python3
"""
Test script to verify Python debug logging is present in all files
"""

import re
from pathlib import Path

def test_python_logging():
    """Check that debug logging has been added to Python files"""
    
    print("🔍 Checking Python files for debug logging...\n")
    
    files_to_check = [
        {
            'path': Path('/Volumes/nvme/Github/bmyCure4MM/chemtools/views.py'),
            'name': 'views.py (Django Views)',
            'patterns': [
                (r'\[DEBUG\].*binding_viz view called', 'binding_viz entry log'),
                (r'\[DEBUG\].*Form validated', 'Form validation log'),
                (r'\[DEBUG\].*API Preferences captured', 'API preferences log'),
                (r'\[DEBUG\].*Job created', 'Job creation log'),
                (r'\[DEBUG\].*Attempting to enqueue', 'Job enqueueing log'),
                (r'\[DEBUG\].*job_detail view called', 'job_detail entry log'),
                (r'\[DEBUG\].*Job loaded', 'Job loading log'),
                (r'\[DEBUG\].*Loading HTML', 'HTML loading log'),
                (r'\[DEBUG\].*Enriching metadata', 'Metadata enrichment log'),
                (r'logger\.info', 'Logger info calls'),
                (r'logger\.debug', 'Logger debug calls'),
            ]
        },
        {
            'path': Path('/Volumes/nvme/Github/bmyCure4MM/chemtools/pdb_api_client.py'),
            'name': 'pdb_api_client.py (API Client)',
            'patterns': [
                (r'\[DEBUG\].*API CALL:', 'API call logging'),
                (r'\[DEBUG\].*URL:', 'URL logging'),
                (r'\[DEBUG\].*Status:', 'Status code logging'),
                (r'\[DEBUG\].*Response time:', 'Response time logging'),
                (r'\[DEBUG\].*Response size:', 'Response size logging'),
                (r'logger\.info', 'Logger info calls'),
                (r'logger\.debug', 'Logger debug calls'),
                (r'import time', 'Time import for timing'),
            ]
        },
        {
            'path': Path('/Volumes/nvme/Github/bmyCure4MM/modules/binding_visualizer/binding_visualizer.py'),
            'name': 'binding_visualizer.py (Core Script)',
            'patterns': [
                (r'DEBUG MODE ACTIVE', 'Debug mode banner'),
                (r'\[DEBUG\].*Starting script', 'Script start log'),
                (r'\[DEBUG\].*Timestamp:', 'Timestamp log'),
                (r'\[DEBUG\].*Python version:', 'Python version log'),
                (r'\[DEBUG\].*Platform:', 'Platform log'),
                (r'\[DEBUG\].*Entered main', 'Main entry log'),
                (r'\[DEBUG\].*Fetching PDB data', 'PDB fetch log'),
                (r'\[DEBUG\].*PDB data fetched in', 'PDB fetch timing'),
                (r'logging\.info', 'Logging info calls'),
            ]
        }
    ]
    
    total_passed = 0
    total_failed = 0
    
    for file_info in files_to_check:
        file_path = file_info['path']
        file_name = file_info['name']
        patterns = file_info['patterns']
        
        print(f"{'='*70}")
        print(f"📁 {file_name}")
        print(f"{'='*70}")
        
        if not file_path.exists():
            print(f"❌ File not found: {file_path}\n")
            total_failed += len(patterns)
            continue
        
        content = file_path.read_text()
        
        file_passed = 0
        file_failed = 0
        
        for pattern, description in patterns:
            if re.search(pattern, content, re.IGNORECASE):
                print(f"✅ {description}")
                file_passed += 1
            else:
                print(f"❌ {description} - Pattern not found: {pattern}")
                file_failed += 1
        
        # Count logging statements
        info_count = len(re.findall(r'logger\.info', content))
        debug_count = len(re.findall(r'logger\.debug', content))
        warning_count = len(re.findall(r'logger\.warning', content))
        error_count = len(re.findall(r'logger\.error', content))
        print_debug_count = len(re.findall(r'print\(Fore\.\w+.*\[DEBUG\]', content))
        
        print(f"\n📊 Logging Statistics:")
        print(f"   logger.info: {info_count}")
        print(f"   logger.debug: {debug_count}")
        print(f"   logger.warning: {warning_count}")
        print(f"   logger.error: {error_count}")
        print(f"   print([DEBUG]): {print_debug_count}")
        print(f"   Total: {info_count + debug_count + warning_count + error_count + print_debug_count}")
        
        print(f"\n✅ Passed: {file_passed}/{len(patterns)}")
        if file_failed > 0:
            print(f"❌ Failed: {file_failed}/{len(patterns)}")
        print()
        
        total_passed += file_passed
        total_failed += file_failed
    
    print(f"{'='*70}")
    print(f"OVERALL RESULTS")
    print(f"{'='*70}")
    print(f"✅ Passed: {total_passed}")
    print(f"❌ Failed: {total_failed}")
    print(f"{'='*70}\n")
    
    if total_failed == 0:
        print("🎉 All Python debug logging checks passed!")
        return True
    else:
        print(f"⚠️ {total_failed} checks failed. Review the output above.")
        return False

if __name__ == '__main__':
    success = test_python_logging()
    exit(0 if success else 1)
