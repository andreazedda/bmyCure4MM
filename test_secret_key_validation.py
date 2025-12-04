#!/usr/bin/env python
"""
Test script to validate Django SECRET_KEY security enforcement.
Run with: python test_secret_key_validation.py
"""

import os
import sys
from pathlib import Path

# Add project to path
BASE_DIR = Path(__file__).resolve().parent
sys.path.insert(0, str(BASE_DIR))

from django.core.exceptions import ImproperlyConfigured

print('=' * 70)
print('DJANGO SECRET_KEY SECURITY VALIDATION TESTS')
print('=' * 70)

tests_passed = 0
tests_total = 4

# Test 1: Missing DJANGO_SECRET_KEY
print('\n[Test 1/4] Missing DJANGO_SECRET_KEY')
os.environ.pop('DJANGO_SECRET_KEY', None)
try:
    if 'mmportal.settings' in sys.modules:
        del sys.modules['mmportal.settings']
    import mmportal.settings as settings
    print('❌ FAIL: Should have raised ImproperlyConfigured')
except ImproperlyConfigured as e:
    print('✅ PASS: Correctly raised ImproperlyConfigured')
    print(f'   Message: {str(e)[:100]}...')
    tests_passed += 1

# Test 2: Insecure DJANGO_SECRET_KEY (contains "insecure")
print('\n[Test 2/4] Insecure DJANGO_SECRET_KEY (contains "insecure")')
os.environ['DJANGO_SECRET_KEY'] = 'django-insecure-test'
try:
    if 'mmportal.settings' in sys.modules:
        del sys.modules['mmportal.settings']
    import mmportal.settings as settings
    print('❌ FAIL: Should have raised ImproperlyConfigured')
except ImproperlyConfigured as e:
    print('✅ PASS: Correctly raised ImproperlyConfigured')
    print(f'   Message: {str(e)[:100]}...')
    tests_passed += 1

# Test 3: Default DJANGO_SECRET_KEY
print('\n[Test 3/4] Default DJANGO_SECRET_KEY')
os.environ['DJANGO_SECRET_KEY'] = 'django-insecure-please-change-me'
try:
    if 'mmportal.settings' in sys.modules:
        del sys.modules['mmportal.settings']
    import mmportal.settings as settings
    print('❌ FAIL: Should have raised ImproperlyConfigured')
except ImproperlyConfigured as e:
    print('✅ PASS: Correctly raised ImproperlyConfigured')
    print(f'   Message: {str(e)[:100]}...')
    tests_passed += 1

# Test 4: Valid DJANGO_SECRET_KEY
print('\n[Test 4/4] Valid strong DJANGO_SECRET_KEY')
os.environ['DJANGO_SECRET_KEY'] = 'valid-production-key-xh8_p9k_2m_5n_qwerty-xyz-123'
try:
    if 'mmportal.settings' in sys.modules:
        del sys.modules['mmportal.settings']
    import mmportal.settings as settings
    print('✅ PASS: Settings loaded successfully')
    print(f'   SECRET_KEY length: {len(settings.SECRET_KEY)} chars')
    tests_passed += 1
except ImproperlyConfigured as e:
    print('❌ FAIL: Should not raise ImproperlyConfigured with valid key')
    print(f'   Error: {e}')

print('\n' + '=' * 70)
if tests_passed == tests_total:
    print(f'✅ SUCCESS: All {tests_passed}/{tests_total} tests passed!')
    print('\nSecurity enforcement is working correctly:')
    print('  • App refuses to start without DJANGO_SECRET_KEY')
    print('  • App refuses keys containing "insecure"')
    print('  • App refuses known default/weak keys')
    print('  • App accepts valid strong keys')
    sys.exit(0)
else:
    print(f'❌ FAILURE: Only {tests_passed}/{tests_total} tests passed')
    sys.exit(1)
print('=' * 70)
