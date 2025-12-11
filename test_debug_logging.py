#!/usr/bin/env python3
"""
Test script to verify comprehensive debug logging is present in job_detail.html
"""

import re
from pathlib import Path

def test_debug_logging():
    """Check that all debug logging has been added to the template"""
    
    template_path = Path('/Volumes/nvme/Github/bmyCure4MM/chemtools/templates/chemtools/job_detail.html')
    
    if not template_path.exists():
        print(f"❌ Template not found: {template_path}")
        return False
    
    content = template_path.read_text()
    
    print("🔍 Checking for debug logging in job_detail.html...\n")
    
    checks = [
        # Banner
        ("Debug Mode Banner", r"DEBUG MODE ACTIVE"),
        
        # Function group headers
        ("Content Reorganization Group", r"\[DEBUG\].*Starting Content Reorganization"),
        ("Button Processing Group", r"\[DEBUG\].*Processing Control Buttons"),
        ("Viewer Availability Group", r"\[DEBUG\].*Checking Viewer Availability"),
        ("Initialization Group", r"\[DEBUG\].*Starting Viewer UI Initialization"),
        
        # Detailed logging
        ("Container Found Log", r"Container found:"),
        ("Info Div Log", r"Info div found:"),
        ("Table Count Log", r"Total tables found:"),
        ("Button Count Log", r"Total buttons:"),
        ("Viewer Ready Log", r"Viewer ready:"),
        
        # Performance tracking
        ("Performance Timer Start", r"console\.time\('reorganizeBindingContent'\)"),
        ("Performance Timer End", r"console\.timeEnd\('reorganizeBindingContent'\)"),
        
        # State checks
        ("DOM Check Object", r"DOM Check:"),
        ("Final State Summary", r"Final State Summary"),
        
        # Click tracking
        ("Click Event Listener", r"addEventListener\('click'"),
        ("Button Click Debug", r"\[DEBUG\].*Button Click"),
        
        # Extract field logging
        ("Extract Field Logging", r"extractField.*\$\{label\}"),
        
        # Summary tables
        ("Console Table", r"console\.table"),
    ]
    
    passed = 0
    failed = 0
    
    for name, pattern in checks:
        if re.search(pattern, content, re.IGNORECASE):
            print(f"✅ {name}")
            passed += 1
        else:
            print(f"❌ {name} - Pattern not found: {pattern}")
            failed += 1
    
    print(f"\n{'='*60}")
    print(f"Results: {passed} passed, {failed} failed")
    print(f"{'='*60}\n")
    
    # Count total console.log statements
    console_logs = len(re.findall(r'console\.(log|group|warn|error|time)', content))
    print(f"📊 Total console logging statements: {console_logs}")
    
    # Count console.group statements
    console_groups = len(re.findall(r'console\.group', content))
    print(f"📦 Console groups: {console_groups}")
    
    # Count debug markers
    debug_markers = len(re.findall(r'\[DEBUG\]', content))
    print(f"🔍 Debug markers: {debug_markers}")
    
    # Count emojis in logs
    emoji_logs = len(re.findall(r'console\.(log|group).*[🔄🎮🔍🚀📊🖱️✅❌⏳]', content))
    print(f"✨ Logs with emojis: {emoji_logs}")
    
    print()
    
    if failed == 0:
        print("🎉 All debug logging checks passed!")
        return True
    else:
        print(f"⚠️ {failed} checks failed. Review the output above.")
        return False

if __name__ == '__main__':
    success = test_debug_logging()
    exit(0 if success else 1)
