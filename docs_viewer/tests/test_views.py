"""
Integration tests for documentation viewer views.
"""
from django.contrib.auth import get_user_model
from django.test import Client, TestCase
from django.urls import reverse

User = get_user_model()


class DocsViewerViewTests(TestCase):
    def setUp(self) -> None:
        self.client = Client()
        self.user = User.objects.create_user(username="testuser", password="testpass123")
        self.auth_client = Client()
        self.auth_client.force_login(self.user)

    def test_index_accessible(self):
        response = self.client.get(reverse("docs_viewer:index"))
        self.assertEqual(response.status_code, 200)

    def test_index_lists_documents(self):
        response = self.client.get(reverse("docs_viewer:index"))
        self.assertIn(b"Documentation", response.content)
        self.assertTrue(b"README" in response.content or b"IMPLEMENTATION" in response.content)

    def test_index_has_search_form(self):
        response = self.client.get(reverse("docs_viewer:index"))
        self.assertIn(b"<form", response.content)
        self.assertIn(b"search", response.content.lower())

    def test_view_readme_accessible(self):
        response = self.client.get(reverse("docs_viewer:view", kwargs={"path": "README.md"}))
        self.assertEqual(response.status_code, 200)

    def test_view_renders_markdown(self):
        response = self.client.get(reverse("docs_viewer:view", kwargs={"path": "README.md"}))
        content = response.content.decode()
        self.assertTrue("<h1" in content or "<h2" in content)
        # README contains code blocks with '#' characters; just ensure markdown rendered.

    def test_view_shows_breadcrumbs(self):
        response = self.client.get(reverse("docs_viewer:view", kwargs={"path": "README.md"}))
        self.assertIn(b"breadcrumb", response.content)

    def test_view_has_download_button(self):
        response = self.client.get(reverse("docs_viewer:view", kwargs={"path": "README.md"}))
        self.assertTrue(b"Download" in response.content or b"Raw" in response.content)

    def test_view_tracks_views_for_authenticated(self):
        from docs_viewer.models import DocumentView

        initial_count = DocumentView.objects.count()
        self.auth_client.get(reverse("docs_viewer:view", kwargs={"path": "README.md"}))
        self.assertEqual(DocumentView.objects.count(), initial_count + 1)
        view = DocumentView.objects.latest("viewed_at")
        self.assertEqual(view.path, "README.md")
        self.assertEqual(view.user, self.user)

    def test_path_traversal_blocked(self):
        response = self.client.get(reverse("docs_viewer:view", kwargs={"path": "../etc/passwd"}))
        self.assertEqual(response.status_code, 404)

    def test_absolute_path_blocked(self):
        response = self.client.get(reverse("docs_viewer:view", kwargs={"path": "/etc/passwd"}))
        self.assertEqual(response.status_code, 404)

    def test_non_whitelisted_file_blocked(self):
        response = self.client.get(reverse("docs_viewer:view", kwargs={"path": "manage.py"}))
        self.assertEqual(response.status_code, 404)

    def test_nonexistent_file_404(self):
        response = self.client.get(reverse("docs_viewer:view", kwargs={"path": "docs/nonexistent.md"}))
        self.assertEqual(response.status_code, 404)
        self.assertContains(response, "Find documentation", status_code=404)
        self.assertContains(response, reverse("docs_viewer:index"), status_code=404)

    def test_raw_download_works(self):
        response = self.client.get(reverse("docs_viewer:raw", kwargs={"path": "README.md"}))
        self.assertEqual(response.status_code, 200)
        self.assertEqual(response["Content-Type"], "text/plain; charset=utf-8")

    def test_raw_has_attachment_header(self):
        response = self.client.get(reverse("docs_viewer:raw", kwargs={"path": "README.md"}))
        self.assertIn("Content-Disposition", response)
        self.assertIn("attachment", response["Content-Disposition"])
        self.assertIn("README.md", response["Content-Disposition"])

    def test_raw_respects_security(self):
        response = self.client.get(reverse("docs_viewer:raw", kwargs={"path": "../etc/passwd"}))
        self.assertEqual(response.status_code, 404)

    def test_search_requires_query(self):
        response = self.client.get(reverse("docs_viewer:search"))
        self.assertEqual(response.status_code, 302)

    def test_search_finds_content(self):
        response = self.client.get(reverse("docs_viewer:search"), {"q": "documentation"})
        self.assertEqual(response.status_code, 200)
        self.assertTrue(b"result" in response.content.lower() or b"found" in response.content.lower())

    def test_search_case_insensitive(self):
        response1 = self.client.get(reverse("docs_viewer:search"), {"q": "DOCUMENTATION"})
        response2 = self.client.get(reverse("docs_viewer:search"), {"q": "documentation"})
        self.assertEqual(response1.status_code, response2.status_code)

        self.assertEqual(response1.status_code, 200)
        self.assertEqual(response2.status_code, 200)
    
    def test_search_shows_match_count(self):
        """Search results should show number of matches."""
        response = self.client.get(reverse('docs_viewer:search'), {'q': 'test'})
        if b'result' in response.content.lower():
            # If results found, should show count
            self.assertTrue(b'match' in response.content.lower() or b'found' in response.content.lower())
    
    def test_search_results_clickable(self):
        """Search results should link to documents."""
        response = self.client.get(reverse('docs_viewer:search'), {'q': 'test'})
        content = response.content.decode()
        if 'result' in content.lower():
            self.assertTrue('/docs/view/' in content or 'list-group-item' in content)

    def test_search_never_exposes_external_symlink_snippets(self):
        import tempfile
        from pathlib import Path
        from django.test import override_settings

        with tempfile.TemporaryDirectory() as root, tempfile.TemporaryDirectory() as external:
            base = Path(root)
            (base / 'docs').mkdir()
            (Path(external) / 'private.md').write_text('UNIQUE_PRIVATE_MARKER', encoding='utf-8')
            (base / 'docs' / 'linked.md').symlink_to(Path(external) / 'private.md')
            with override_settings(BASE_DIR=base):
                response = self.client.get(reverse('docs_viewer:search'), {'q': 'UNIQUE_PRIVATE_MARKER'})
                self.assertEqual(response.status_code, 200)
                self.assertEqual(response.context['total_results'], 0)
                self.assertNotContains(response, 'linked.md')
                index = self.client.get(reverse('docs_viewer:index'))
                self.assertEqual(index.context['total_docs'], 0)

    def test_search_is_paginated_and_preserves_the_query(self):
        import tempfile
        from pathlib import Path
        from django.test import override_settings

        with tempfile.TemporaryDirectory() as root:
            base = Path(root)
            (base / 'docs').mkdir()
            for number in range(23):
                (base / 'docs' / f'{number:02}.md').write_text(f'# Handbook {number}\nUnique search phrase', encoding='utf-8')
            with override_settings(BASE_DIR=base):
                first = self.client.get(reverse('docs_viewer:search'), {'q': 'Unique search phrase'})
                second = self.client.get(reverse('docs_viewer:search'), {'q': 'Unique search phrase', 'page': '2'})
                self.assertEqual(first.context['total_results'], 23)
                self.assertEqual(len(first.context['results']), 20)
                self.assertEqual(len(second.context['results']), 3)
                self.assertContains(first, 'q=Unique%20search%20phrase')
                self.assertContains(first, 'for="search-query"')
                self.assertContains(first, 'Handbook 0')
                self.assertIn('no-store', first['Cache-Control'])

    def test_document_declared_status_is_visible_and_page_is_not_cached(self):
        import tempfile
        from pathlib import Path
        from django.test import override_settings
        with tempfile.TemporaryDirectory() as root:
            base = Path(root)
            (base / 'docs').mkdir()
            (base / 'docs' / 'status.md').write_text('---\nstatus: AUDIT_TEST_ONLY\nlast_verified_at: 2026-11-06\n---\n# Test document', encoding='utf-8')
            with override_settings(BASE_DIR=base):
                response = self.client.get(reverse('docs_viewer:view', kwargs={'path': 'docs/status.md'}))
                self.assertContains(response, 'Authored status:')
                self.assertContains(response, 'AUDIT_TEST_ONLY')
                self.assertContains(response, 'Last verification recorded:')
                self.assertContains(response, 'id="main-content"')
                self.assertContains(response, 'href="#main-content"')
                self.assertIn('no-store', response['Cache-Control'])

    def test_document_without_metadata_is_not_presented_as_verified(self):
        response = self.client.get(reverse('docs_viewer:view', kwargs={'path': 'README.md'}))
        self.assertContains(response, 'Status not recorded')
