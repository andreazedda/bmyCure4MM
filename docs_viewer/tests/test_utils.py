"""
Unit tests for documentation viewer utility functions.
"""
import os
import tempfile
from pathlib import Path
from django.test import SimpleTestCase
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


class PathSecurityTests(SimpleTestCase):
    """Test path validation and security controls."""
    
    def test_safe_path_allows_whitelisted_files(self):
        """Whitelist files should be accessible."""
        self.assertTrue(is_safe_path("README.md"))
        self.assertTrue(is_safe_path("docs/governance/INTENDED_USE.md"))
    
    def test_safe_path_rejects_path_traversal(self):
        """Path traversal attempts should be blocked."""
        self.assertFalse(is_safe_path("../etc/passwd"))
        self.assertFalse(is_safe_path("docs/../../etc/passwd"))
        self.assertFalse(is_safe_path("./docs/../../../etc/passwd"))
    
    def test_safe_path_rejects_absolute_paths(self):
        """Absolute paths should be blocked."""
        self.assertFalse(is_safe_path("/etc/passwd"))
        self.assertFalse(is_safe_path("/var/www/html/"))
    
    def test_safe_path_rejects_non_whitelisted(self):
        """Files outside whitelist should be blocked."""
        self.assertFalse(is_safe_path("manage.py"))
        self.assertFalse(is_safe_path("mmportal/settings.py"))
        self.assertFalse(is_safe_path("db.sqlite3"))
    
    def test_safe_path_rejects_nonexistent_files(self):
        """Non-existent files should be blocked."""
        self.assertFalse(is_safe_path("docs/nonexistent.md"))
    
    def test_sanitize_path_removes_dangerous_chars(self):
        """Sanitization should remove null bytes and normalize."""
        self.assertEqual(sanitize_path("docs\0/file.md"), "docs/file.md")
        self.assertEqual(sanitize_path("docs\\file.md"), "docs/file.md")
        self.assertEqual(sanitize_path("//docs///file.md"), "docs/file.md")
        self.assertEqual(sanitize_path("/docs/file.md/"), "docs/file.md")


class MarkdownRenderingTests(SimpleTestCase):
    """Test markdown to HTML conversion."""
    
    def test_render_markdown_basic(self):
        """Basic markdown should render correctly."""
        content = "# Hello World\n\nThis is **bold** text."
        html, toc = render_markdown(content)

        self.assertIn("<h1", html)
        self.assertIn("Hello World", html)
        self.assertIn("<strong>bold</strong>", html)
    
    def test_render_markdown_code_blocks(self):
        """Code blocks should render with syntax highlighting."""
        content = """
```python
def hello():
    print("world")
```
"""
        html, toc = render_markdown(content)
        self.assertTrue("<code" in html or "<pre" in html)
    
    def test_render_markdown_tables(self):
        """Tables should render correctly."""
        content = """
| Header 1 | Header 2 |
|----------|----------|
| Cell 1   | Cell 2   |
"""
        html, toc = render_markdown(content)
        self.assertIn("<table", html)
        self.assertIn("<th", html)
        self.assertIn("<td", html)
    
    def test_render_markdown_generates_toc(self):
        """TOC should be generated from headings."""
        content = """
# Title
## Section 1
### Subsection 1.1
## Section 2
"""
        html, toc = render_markdown(content)
        self.assertIn("Section 1", toc)
        self.assertIn("Section 2", toc)


class TitleExtractionTests(SimpleTestCase):
    """Test title extraction from markdown."""
    
    def test_extract_title_from_h1(self):
        """First H1 should be extracted as title."""
        content = "# Main Title\n\nSome content\n\n## Subtitle"
        title = extract_title(content)
        self.assertEqual(title, "Main Title")
    
    def test_extract_title_no_heading(self):
        """Should return default when no heading found."""
        content = "Just plain text without headings."
        title = extract_title(content)
        self.assertEqual(title, "Untitled Document")


class BreadcrumbsTests(SimpleTestCase):
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


class ReferenceRenderingTests(SimpleTestCase):
    def test_math_figures_links_and_admonitions_use_existing_public_assets(self):
        source = "docs/guides/mathematical_models_en.md"
        content = Path(settings.BASE_DIR, source).read_text()
        html, toc = render_markdown(content, source_path=source)
        self.assertIn('class="arithmatex"', html)
        self.assertIn('class="admonition warning"', html)
        self.assertIn('/static/docs_images/models/logistic_growth.svg', html)
        self.assertIn('class="mermaid"', html)
        self.assertNotIn('toctitle', toc)
        linked, _ = render_markdown('[Model](../models/REGISTRY.md#registry)', source_path=source)
        self.assertIn('/docs/view/docs/models/REGISTRY.md#registry', linked)
        self.assertEqual(validate_markdown_links('[Model](../models/REGISTRY.md#registry)', source), [])

    def test_code_is_highlighted_and_front_matter_is_not_visible(self):
        html, _ = render_markdown('---\nstatus: CURRENT_PARTIAL\n---\n\n# Title\n\n```python\ndef hello():\n    return 1\n```')
        self.assertNotIn('CURRENT_PARTIAL', html)
        self.assertIn('class="k"', html)

    def test_file_whitelist_rejects_assets_and_external_symlinks(self):
        self.assertFalse(is_safe_path('docs/assets/javascripts/mathjax.js'))
        with tempfile.TemporaryDirectory() as root:
            base = Path(root, 'repo')
            Path(base, 'docs').mkdir(parents=True)
            private = Path(root, 'private.md')
            private.write_text('Must not be served')
            Path(base, 'docs/linked.md').symlink_to(private)
            with self.settings(BASE_DIR=base):
                self.assertFalse(is_safe_path('docs/linked.md'))

    def test_symlink_destination_must_be_whitelisted_markdown(self):
        with tempfile.TemporaryDirectory() as root:
            base = Path(root)
            Path(base, 'docs').mkdir()
            Path(base, 'private.md').write_text('Synthetic private sentinel')
            Path(base, 'docs/private.txt').write_text('Synthetic non-document sentinel')
            Path(base, 'docs/public.md').write_text('# Public documentation')
            Path(base, 'docs/internal.md').symlink_to(base / 'private.md')
            Path(base, 'docs/wrong-type.md').symlink_to(base / 'docs/private.txt')
            Path(base, 'docs/allowed.md').symlink_to(base / 'docs/public.md')
            Path(base, 'docs/loop.md').symlink_to(base / 'docs/loop.md')
            with self.settings(BASE_DIR=base):
                self.assertFalse(is_safe_path('docs/internal.md'))
                self.assertFalse(is_safe_path('docs/wrong-type.md'))
                self.assertFalse(is_safe_path('docs/loop.md'))
                self.assertTrue(is_safe_path('docs/allowed.md'))
                self.assertEqual(get_allowed_doc_paths(), ['docs/allowed.md', 'docs/public.md'])


class DocumentationSourceReferenceTests(SimpleTestCase):
    def test_public_source_reference_uses_image_revision_and_does_not_expose_code(self):
        from unittest.mock import patch
        with patch.dict(os.environ, {'BMY_SOURCE_SHA': 'a' * 40}):
            html, _ = render_markdown('[Source](../../accounts/services.py#L12)', source_path='docs/development/DEVELOPMENT.md')
            self.assertIn('https://github.com/andreazedda/bmyCure4MM/blob/' + 'a' * 40 + '/accounts/services.py#L12', html)
            self.assertNotIn('/docs/view/accounts/services.py', html)

    def test_source_reference_rejects_data_and_external_symlinks(self):
        from docs_viewer.utils import repository_source_url
        from django.test import override_settings
        with tempfile.TemporaryDirectory() as root, tempfile.TemporaryDirectory() as external:
            base = Path(root)
            (base / 'accounts').mkdir()
            (base / 'accounts' / 'private.json').write_text('{}')
            (Path(external) / 'secret.py').write_text('# private test only')
            (base / 'accounts' / 'secret.py').symlink_to(Path(external) / 'secret.py')
            with override_settings(BASE_DIR=base):
                self.assertEqual(repository_source_url('accounts/private.json'), '')
                self.assertEqual(repository_source_url('accounts/secret.py'), '')

    def test_actual_operator_runbook_is_readable_via_docs_allowlist(self):
        self.assertTrue(is_safe_path('deploy/k8s/DEPLOY_K3S.md'))
        html, _ = render_markdown('[Operations](../../deploy/k8s/DEPLOY_K3S.md)', source_path='docs/development/DEVELOPMENT.md')
        self.assertIn('/docs/view/deploy/k8s/DEPLOY_K3S.md', html)

    def test_known_and_declared_document_languages_do_not_guess_mixed_docs(self):
        from docs_viewer.utils import document_language
        self.assertEqual(document_language('# Modelli', 'docs/guides/mathematical_models.md'), 'it')
        self.assertEqual(document_language('# Models', 'docs/guides/mathematical_models_en.md'), 'en')
        self.assertEqual(document_language('---\nlanguage: it\n---\n# Guida', 'docs/development/test.md'), 'it')
        self.assertEqual(document_language('# English / Italiano', 'docs/reference/database_objects.md'), '')
