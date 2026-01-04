#!/usr/bin/env python3
"""
Diagnostic script to check if interactive controls are properly generated
"""

import os
import sys
import django

sys.path.insert(0, '/Volumes/nvme/Github/bmyCure4MM')
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'mmportal.settings')
django.setup()

from chemtools.models import ChemJob

def check_html_content():
    """Check the most recent BIND job's HTML content"""
    print("=" * 60)
    print("CHECKING INTERACTIVE CONTROLS IN RECENT JOB")
    print("=" * 60)
    
    # Get most recent BIND job with HTML content
    jobs = ChemJob.objects.filter(kind='BIND', out_html__isnull=False).order_by('-created')
    
    if not jobs.exists():
        print("✗ No BIND jobs with HTML content found")
        return
    
    job = jobs.first()
    print(f"✓ Found job #{job.pk}: {job.input_a}")
    print(f"  Created: {job.created}")
    print()
    
    # Read HTML content
    try:
        with open(job.out_html.path, 'r') as f:
            html = f.read()
        
        print(f"HTML file size: {len(html):,} bytes")
        print()
        
        # Check for key elements
        checks = {
            '#bv-controls': html.count('id="bv-controls"'),
            'button elements': html.count('<button'),
            'onclick attributes': html.count('onclick='),
            'disabled buttons': html.count('disabled>'),
            'enableControls function': html.count('function enableControls'),
            'window.viewer': html.count('window.viewer'),
            'highlightMutations': html.count('highlightMutations'),
            '3dmolviewer div': html.count('3dmolviewer'),
        }
        
        print("Element Check:")
        for element, count in checks.items():
            status = "✓" if count > 0 else "✗"
            print(f"  {status} {element}: {count} occurrence(s)")
        print()
        
        # Extract a sample button
        import re
        button_match = re.search(r'<button[^>]*onclick="([^"]+)"[^>]*>([^<]+)</button>', html)
        if button_match:
            onclick = button_match.group(1)
            text = button_match.group(2)
            print(f"Sample button:")
            print(f"  Text: {text}")
            print(f"  Onclick: {onclick[:100]}...")
            print()
        
        # Check if buttons are in the controls panel
        controls_match = re.search(r'<div id="bv-controls"[^>]*>(.*?)</div>\s*<script>', html, re.DOTALL)
        if controls_match:
            controls_html = controls_match.group(1)
            controls_buttons = controls_html.count('<button')
            print(f"✓ Found #bv-controls with {controls_buttons} buttons")
        else:
            print("✗ Could not find #bv-controls structure")
        
    except Exception as e:
        print(f"✗ Error reading HTML: {e}")
        import traceback
        traceback.print_exc()

if __name__ == '__main__':
    check_html_content()
