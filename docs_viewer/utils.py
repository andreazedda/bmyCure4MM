"""
Utility functions for documentation viewer with security controls.
"""
import os
import re
import posixpath
from pathlib import Path
from typing import List, Tuple, Optional
from django.conf import settings
import markdown
from markdown.extensions import Extension
from markdown.treeprocessors import Treeprocessor
from django.templatetags.static import static
from urllib.parse import urlsplit, urlunsplit, unquote, quote
import yaml
from pymdownx.superfences import fence_div_format


# Whitelist of allowed documentation directories
ALLOWED_DOC_PATHS = [
    'docs/',
    'tests/',
    'README.md',
    'IMPLEMENTATION_LOG.md',
    'IMPLEMENTATION_SUMMARY.md',
    'deploy/k8s/DEPLOY_K3S.md',
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
    if '..' in normalized or normalized.startswith('/') or Path(normalized).suffix.lower() not in {'.md', '.markdown'}:
        return False
    
    # Check whitelist
    allowed = any(
        (normalized.startswith(allowed_path) if allowed_path.endswith('/') else normalized == allowed_path)
        for allowed_path in ALLOWED_DOC_PATHS
    )
    if not allowed:
        return False
    
    # Resolve symlinks before checking the destination against the same whitelist.
    # Staying within BASE_DIR alone would permit docs/link.md -> private.md.
    full_path = os.path.join(settings.BASE_DIR, normalized)
    try:
        resolved = Path(full_path).resolve().relative_to(Path(settings.BASE_DIR).resolve()).as_posix()
    except (OSError, RuntimeError, ValueError):
        return False
    if Path(resolved).suffix.lower() not in {'.md', '.markdown'} or not any(
        (resolved.startswith(allowed_path) if allowed_path.endswith('/') else resolved == allowed_path)
        for allowed_path in ALLOWED_DOC_PATHS
    ):
        return False
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
            if is_safe_path(allowed_path):
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
                        if is_safe_path(rel_path):
                            docs.append(rel_path)
    
    return sorted(docs)


def document_metadata(content: str) -> dict:
    """Expose authored status/date without interpreting missing metadata as approval."""
    if not content.startswith('---\n'):
        return {}
    closing = content.find('\n---', 4)
    if closing < 0:
        return {}
    try:
        data = yaml.safe_load(content[4:closing])
    except yaml.YAMLError:
        return {}
    if not isinstance(data, dict):
        return {}
    return {key: str(data[key])[:200] for key in ('status', 'last_verified_at', 'audience', 'language', 'lang') if data.get(key)}


def document_language(content: str, path: str) -> str:
    """Prefer declared language; do not guess a single language for mixed documents."""
    metadata = document_metadata(content)
    declared = metadata.get('language') or metadata.get('lang', '')
    if re.fullmatch(r'(?:en|it)(?:-[A-Za-z]{2,8})?', declared):
        return declared
    if path.startswith('docs/it/') or path.endswith('_it.md') or path == 'docs/guides/mathematical_models.md':
        return 'it'
    if path.startswith(('docs/en/', 'docs/development/')) or path.endswith('_en.md'):
        return 'en'
    if path in {'docs/architecture/CURRENT_SYSTEM.md', 'docs/operations/AUTHENTICATION_AND_ACCESS.md',
                'docs/operations/UX_AUDIT_2026_11_06.md', 'deploy/k8s/DEPLOY_K3S.md'}:
        return 'en'
    return ''


def repository_source_url(path: str) -> str:
    """Link known source types to this image's public revision; never serve code/files."""
    roots = ('accounts/', 'simulator/', 'mmportal/', 'docs_viewer/', 'twin_engine/',
             'clinic/', 'chemtools/', 'api/', 'scripts/', 'tests/', 'deploy/', '.github/')
    root_files = {'pyproject.toml', 'uv.lock', 'Dockerfile', 'manage.py', 'mkdocs.yml'}
    candidate = Path(settings.BASE_DIR, path)
    if (path not in root_files and not (path.startswith(roots) and candidate.suffix in {'.py', '.toml', '.yaml', '.yml', '.md'})):
        return ''
    if not candidate.is_file() or not candidate.resolve().is_relative_to(Path(settings.BASE_DIR).resolve()):
        return ''
    revision = os.environ.get('BMY_SOURCE_SHA', '')
    if not re.fullmatch(r'[a-f0-9]{40}', revision):
        revision = 'master'
    return 'https://github.com/andreazedda/bmyCure4MM/blob/' + revision + '/' + quote(path, safe='/')


class DocumentationReferences(Treeprocessor):
    def __init__(self, md, source_path):
        super().__init__(md)
        self.source_path = source_path

    def run(self, root):
        for element in root.iter():
            attribute = 'src' if element.tag == 'img' else 'href' if element.tag == 'a' else None
            if not attribute:
                continue
            original = element.get(attribute, '')
            parts = urlsplit(original)
            if parts.scheme or parts.netloc or not parts.path or parts.path.startswith('/'):
                continue
            resolved = posixpath.normpath(posixpath.join(posixpath.dirname(self.source_path), unquote(parts.path)))
            if attribute == 'href' and is_safe_path(resolved):
                element.set(attribute, urlunsplit(('', '', '/docs/view/' + resolved, parts.query, parts.fragment)))
            elif attribute == 'href' and repository_source_url(resolved):
                element.set(attribute, repository_source_url(resolved) + ('#' + parts.fragment if parts.fragment else ''))
            elif attribute == 'src' and resolved.startswith('docs/assets/images/'):
                image_root = Path(settings.BASE_DIR, 'docs/assets/images').resolve()
                image = Path(settings.BASE_DIR, resolved).resolve()
                if image.is_file() and image.is_relative_to(image_root):
                    element.set(attribute, static('docs_images/' + resolved.removeprefix('docs/assets/images/')))


class DocumentationExtension(Extension):
    def __init__(self, source_path):
        self.source_path = source_path
        super().__init__()

    def extendMarkdown(self, md):
        md.treeprocessors.register(DocumentationReferences(md, self.source_path), 'documentation_references', 5)


def render_markdown(content: str, source_path: str = '') -> tuple[str, str]:
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
            'meta',
            'admonition',
            'attr_list',
            'def_list',
            'pymdownx.superfences',
            'pymdownx.arithmatex',
            'tables',
            'codehilite',
            'nl2br',
            'sane_lists',
            DocumentationExtension(source_path),
        ],
        extension_configs={
            'codehilite': {
                'css_class': 'highlight',
                'linenums': False,
            },
            'pymdownx.arithmatex': {'generic': True},
            'pymdownx.superfences': {
                'custom_fences': [{'name': 'mermaid', 'class': 'mermaid', 'format': fence_div_format}],
            },
            'toc': {
                'permalink': True,
                'permalink_class': 'headerlink',
                'title': '',
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
        parts = urlsplit(link_url)
        if parts.scheme or parts.netloc or not parts.path:
            continue
        
        # Resolve relative link
        link_path = os.path.normpath(
            os.path.join(os.path.dirname(base_path), unquote(parts.path))
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
