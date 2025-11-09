"""
Unit tests for documentation viewer utility functions.
"""
import os
import pytest
from django.conf import settings
from docs_viewer.utils import (
    is_safe_path,
    sanitize_path,
    get_allowed_doc_paths,
    render_markdown,
    extract_title,
    get_breadcrumbs,
    validate_markdown_links,
)


class TestPathSecurity:
    """Test path validation and security controls."""
    
    def test_safe_path_allows_whitelisted_files(self):
        """Whitelist files should be accessible."""
        assert is_safe_path('README.md') is True
        assert is_safe_path('IMPLEMENTATION_LOG.md') is True
    
    def test_safe_path_rejects_path_traversal(self):
        """Path traversal attempts should be blocked."""
        assert is_safe_path('../etc/passwd') is False
        assert is_safe_path('docs/../../etc/passwd') is False
        assert is_safe_path('./docs/../../../etc/passwd') is False
    
    def test_safe_path_rejects_absolute_paths(self):
        """Absolute paths should be blocked."""
        assert is_safe_path('/etc/passwd') is False
        assert is_safe_path('/var/www/html/') is False
    
    def test_safe_path_rejects_non_whitelisted(self):
        """Files outside whitelist should be blocked."""
        assert is_safe_path('manage.py') is False
        assert is_safe_path('mmportal/settings.py') is False
        assert is_safe_path('db.sqlite3') is False
    
    def test_safe_path_rejects_nonexistent_files(self):
        """Non-existent files should be blocked."""
        assert is_safe_path('docs/nonexistent.md') is False
    
    def test_sanitize_path_removes_dangerous_chars(self):
        """Sanitization should remove null bytes and normalize."""
        assert sanitize_path('docs\0/file.md') == 'docs/file.md'
        assert sanitize_path('docs\\file.md') == 'docs/file.md'
        assert sanitize_path('//docs///file.md') == 'docs/file.md'
        assert sanitize_path('/docs/file.md/') == 'docs/file.md'


class TestMarkdownRendering:
    """Test markdown to HTML conversion."""
    
    def test_render_markdown_basic(self):
        """Basic markdown should render correctly."""
        content = "# Hello World\n\nThis is **bold** text."
        html, toc = render_markdown(content)
        
        assert '<h1' in html
        assert 'Hello World' in html
        assert '<strong>bold</strong>' in html
    
    def test_render_markdown_code_blocks(self):
        """Code blocks should render with syntax highlighting."""
        content = """
```python
def hello():
    print("world")
```
"""
        html, toc = render_markdown(content)
        assert '<code' in html or '<pre' in html
    
    def test_render_markdown_tables(self):
        """Tables should render correctly."""
        content = """
| Header 1 | Header 2 |
|----------|----------|
| Cell 1   | Cell 2   |
"""
        html, toc = render_markdown(content)
        assert '<table' in html
        assert '<th' in html
        assert '<td' in html
    
    def test_render_markdown_generates_toc(self):
        """TOC should be generated from headings."""
        content = """
# Title
## Section 1
### Subsection 1.1
## Section 2
"""
        html, toc = render_markdown(content)
        assert 'Section 1' in toc
        assert 'Section 2' in toc


class TestTitleExtraction:
    """Test title extraction from markdown."""
    
    def test_extract_title_from_h1(self):
        """First H1 should be extracted as title."""
        content = "# Main Title\n\nSome content\n\n## Subtitle"
        title = extract_title(content)
        assert title == "Main Title"
    
    def test_extract_title_no_heading(self):
        """Should return default when no heading found."""
        content = "Just plain text without headings."
        title = extract_title(content)
        assert title == "Untitled Document"


class TestBreadcrumbs:
    """Test breadcrumb generation."""
    
    def test_breadcrumbs_root_file(self):
        """Root files should have minimal breadcrumbs."""
        crumbs = get_breadcrumbs('README.md')
        assert len(crumbs) == 2
        assert crumbs[0] == ('Documentation', 'index')
        assert crumbs[1] == ('README.md', None)
    
    def test_breadcrumbs_nested_file(self):
        """Nested files should show full path."""
        crumbs = get_breadcrumbs('docs/en/guide.md')
        assert len(crumbs) == 4
        assert crumbs[0] == ('Documentation', 'index')
        assert crumbs[1][0] == 'docs'
        assert crumbs[2][0] == 'en'
        assert crumbs[3] == ('guide.md', None)


class TestLinkValidation:
    """Test broken link detection."""
    
    def test_validate_links_external_skipped(self):
        """External links should be skipped."""
        content = "[Google](https://google.com) and [Email](mailto:test@example.com)"
        errors = validate_markdown_links(content, 'README.md')
        assert len(errors) == 0
    
    def test_validate_links_anchor_skipped(self):
        """Anchor links should be skipped."""
        content = "[Jump](#section)"
        errors = validate_markdown_links(content, 'README.md')
        assert len(errors) == 0
    
    def test_validate_links_broken_detected(self):
        """Broken relative links should be detected."""
        content = "[Broken](nonexistent.md)"
        errors = validate_markdown_links(content, 'README.md')
        assert len(errors) > 0
        assert 'Broken' in errors[0]


class TestGetAllowedPaths:
    """Test document discovery."""
    
    def test_get_allowed_paths_returns_list(self):
        """Should return list of allowed documents."""
        docs = get_allowed_doc_paths()
        assert isinstance(docs, list)
        assert len(docs) > 0
    
    def test_get_allowed_paths_includes_markdown(self):
        """Should include .md files."""
        docs = get_allowed_doc_paths()
        md_files = [d for d in docs if d.endswith('.md')]
        assert len(md_files) > 0
    
    def test_get_allowed_paths_sorted(self):
        """Results should be sorted."""
        docs = get_allowed_doc_paths()
        assert docs == sorted(docs)
