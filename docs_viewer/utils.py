"""
Utility functions for documentation viewer with security controls.
"""
import os
import re
from pathlib import Path
from typing import List, Tuple, Optional
from django.conf import settings
import markdown
from markdown.extensions import toc, fenced_code, tables, codehilite


# Whitelist of allowed documentation directories
ALLOWED_DOC_PATHS = [
    'docs/',
    'tests/',
    'README.md',
    'IMPLEMENTATION_LOG.md',
    'IMPLEMENTATION_SUMMARY.md',
]


def is_safe_path(path: str) -> bool:
    """
    Validate that path is safe and within allowed directories.
    
    Security checks:
    - No path traversal (../)
    - Must be in whitelist
    - Must exist
    - Must be a file (not directory)
    """
    # Normalize path
    normalized = os.path.normpath(path)
    
    # Check for path traversal
    if '..' in normalized or normalized.startswith('/'):
        return False
    
    # Check whitelist
    allowed = any(
        normalized.startswith(allowed_path) or normalized == allowed_path
        for allowed_path in ALLOWED_DOC_PATHS
    )
    if not allowed:
        return False
    
    # Check file exists
    full_path = os.path.join(settings.BASE_DIR, normalized)
    if not os.path.exists(full_path):
        return False
    
    # Must be a file
    if not os.path.isfile(full_path):
        return False
    
    return True


def get_allowed_doc_paths() -> List[str]:
    """Get list of all allowed documentation files."""
    docs = []
    base_dir = settings.BASE_DIR
    
    for allowed_path in ALLOWED_DOC_PATHS:
        full_path = os.path.join(base_dir, allowed_path)
        
        if os.path.isfile(full_path):
            docs.append(allowed_path)
        elif os.path.isdir(full_path):
            # Recursively find markdown files
            for root, dirs, files in os.walk(full_path):
                for file in files:
                    if file.endswith(('.md', '.markdown')):
                        rel_path = os.path.relpath(
                            os.path.join(root, file),
                            base_dir
                        )
                        docs.append(rel_path)
    
    return sorted(docs)


def render_markdown(content: str) -> str:
    """
    Render markdown to HTML with syntax highlighting and TOC.
    
    Extensions:
    - toc: Table of contents
    - fenced_code: Code blocks
    - tables: GitHub-style tables
    - codehilite: Syntax highlighting
    """
    md = markdown.Markdown(
        extensions=[
            'toc',
            'fenced_code',
            'tables',
            'codehilite',
            'nl2br',
            'sane_lists',
        ],
        extension_configs={
            'codehilite': {
                'css_class': 'highlight',
                'linenums': False,
            },
            'toc': {
                'permalink': True,
                'permalink_class': 'headerlink',
                'title': 'Table of Contents',
            }
        }
    )
    
    html = md.convert(content)
    toc_html = md.toc if hasattr(md, 'toc') else ''
    
    return html, toc_html


def extract_title(content: str) -> str:
    """Extract title from markdown (first H1 or filename)."""
    # Look for first # heading
    match = re.search(r'^#\s+(.+)$', content, re.MULTILINE)
    if match:
        return match.group(1).strip()
    
    # Fallback to empty
    return "Untitled Document"


def get_breadcrumbs(path: str) -> List[Tuple[str, Optional[str]]]:
    """Generate breadcrumb navigation for path.
    
    Returns:
        List of tuples (name, None) - URLs are generated in template
    """
    parts = path.split('/')
    breadcrumbs = [('Documentation', 'index')]
    
    current_path = ''
    for i, part in enumerate(parts[:-1]):  # Exclude filename
        current_path = os.path.join(current_path, part)
        breadcrumbs.append((part, None))
    
    # Add current file without link
    breadcrumbs.append((parts[-1], None))
    
    return breadcrumbs


def validate_markdown_links(content: str, base_path: str) -> List[str]:
    """
    Check for broken links in markdown content.
    
    Returns list of broken link errors.
    """
    errors = []
    base_dir = settings.BASE_DIR
    
    # Find all markdown links [text](url)
    link_pattern = r'\[([^\]]+)\]\(([^\)]+)\)'
    matches = re.finditer(link_pattern, content)
    
    for match in matches:
        link_text = match.group(1)
        link_url = match.group(2)
        
        # Skip external links
        if link_url.startswith(('http://', 'https://', 'mailto:', '#')):
            continue
        
        # Resolve relative link
        link_path = os.path.normpath(
            os.path.join(os.path.dirname(base_path), link_url)
        )
        full_path = os.path.join(base_dir, link_path)
        
        if not os.path.exists(full_path):
            errors.append(f"Broken link: [{link_text}]({link_url})")
    
    return errors


def sanitize_path(path: str) -> str:
    """Sanitize user input path."""
    # Remove any dangerous characters
    path = path.replace('\0', '')
    
    # Normalize slashes
    path = path.replace('\\', '/')
    
    # Remove multiple slashes
    path = re.sub(r'/+', '/', path)
    
    # Remove leading/trailing slashes
    path = path.strip('/')
    
    return path
