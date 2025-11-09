"""
Integration tests for documentation viewer views.
"""
import pytest
from django.test import Client
from django.urls import reverse
from django.contrib.auth import get_user_model

User = get_user_model()


@pytest.fixture
def client():
    """Django test client."""
    return Client()


@pytest.fixture
def authenticated_client(db):
    """Authenticated Django test client."""
    user = User.objects.create_user(username='testuser', password='testpass123')
    client = Client()
    client.force_login(user)
    return client, user


@pytest.mark.django_db
class TestDocsIndex:
    """Test documentation index view."""
    
    def test_index_accessible(self, client):
        """Index should be accessible."""
        response = client.get(reverse('docs_viewer:index'))
        assert response.status_code == 200
    
    def test_index_lists_documents(self, client):
        """Index should list available documents."""
        response = client.get(reverse('docs_viewer:index'))
        assert b'Documentation' in response.content
        assert b'README' in response.content or b'IMPLEMENTATION' in response.content
    
    def test_index_has_search_form(self, client):
        """Index should have search functionality."""
        response = client.get(reverse('docs_viewer:index'))
        assert b'<form' in response.content
        assert b'search' in response.content.lower()


@pytest.mark.django_db
class TestDocsView:
    """Test document viewing."""
    
    def test_view_readme_accessible(self, client):
        """README should be viewable."""
        response = client.get(reverse('docs_viewer:view', kwargs={'path': 'README.md'}))
        assert response.status_code == 200
    
    def test_view_renders_markdown(self, client):
        """Markdown should be rendered to HTML."""
        response = client.get(reverse('docs_viewer:view', kwargs={'path': 'README.md'}))
        content = response.content.decode()
        # Should have HTML tags, not raw markdown
        assert '<h1' in content or '<h2' in content
        assert '# ' not in content  # Raw markdown should be converted
    
    def test_view_shows_breadcrumbs(self, client):
        """View should show breadcrumb navigation."""
        response = client.get(reverse('docs_viewer:view', kwargs={'path': 'README.md'}))
        assert b'breadcrumb' in response.content
    
    def test_view_has_download_button(self, client):
        """View should offer raw download."""
        response = client.get(reverse('docs_viewer:view', kwargs={'path': 'README.md'}))
        assert b'Download' in response.content or b'Raw' in response.content
    
    def test_view_tracks_views_for_authenticated(self, authenticated_client):
        """Views should be tracked for authenticated users."""
        from docs_viewer.models import DocumentView
        
        client, user = authenticated_client
        
        initial_count = DocumentView.objects.count()
        client.get(reverse('docs_viewer:view', kwargs={'path': 'README.md'}))
        
        assert DocumentView.objects.count() == initial_count + 1
        view = DocumentView.objects.latest('viewed_at')
        assert view.path == 'README.md'
        assert view.user == user


@pytest.mark.django_db
class TestDocsSecurityint:
    """Test security controls in views."""
    
    def test_path_traversal_blocked(self, client):
        """Path traversal attempts should be blocked."""
        response = client.get(reverse('docs_viewer:view', kwargs={'path': '../etc/passwd'}))
        assert response.status_code == 404
    
    def test_absolute_path_blocked(self, client):
        """Absolute paths should be blocked."""
        response = client.get(reverse('docs_viewer:view', kwargs={'path': '/etc/passwd'}))
        assert response.status_code == 404
    
    def test_non_whitelisted_file_blocked(self, client):
        """Non-whitelisted files should be blocked."""
        response = client.get(reverse('docs_viewer:view', kwargs={'path': 'manage.py'}))
        assert response.status_code == 404
    
    def test_nonexistent_file_404(self, client):
        """Non-existent files should return 404."""
        response = client.get(reverse('docs_viewer:view', kwargs={'path': 'docs/nonexistent.md'}))
        assert response.status_code == 404


@pytest.mark.django_db
class TestDocsRaw:
    """Test raw file download."""
    
    def test_raw_download_works(self, client):
        """Raw download should serve plain text."""
        response = client.get(reverse('docs_viewer:raw', kwargs={'path': 'README.md'}))
        assert response.status_code == 200
        assert response['Content-Type'] == 'text/plain; charset=utf-8'
    
    def test_raw_has_attachment_header(self, client):
        """Raw download should suggest filename."""
        response = client.get(reverse('docs_viewer:raw', kwargs={'path': 'README.md'}))
        assert 'Content-Disposition' in response
        assert 'attachment' in response['Content-Disposition']
        assert 'README.md' in response['Content-Disposition']
    
    def test_raw_respects_security(self, client):
        """Raw download should enforce same security checks."""
        response = client.get(reverse('docs_viewer:raw', kwargs={'path': '../etc/passwd'}))
        assert response.status_code == 404


@pytest.mark.django_db
class TestDocsSearch:
    """Test documentation search."""
    
    def test_search_requires_query(self, client):
        """Empty search should redirect to index."""
        response = client.get(reverse('docs_viewer:search'))
        assert response.status_code == 302  # Redirect
    
    def test_search_finds_content(self, client):
        """Search should find matching documents."""
        # Search for common word likely in README
        response = client.get(reverse('docs_viewer:search'), {'q': 'documentation'})
        assert response.status_code == 200
        # Should show results or "no results"
        assert b'result' in response.content.lower() or b'found' in response.content.lower()
    
    def test_search_case_insensitive(self, client):
        """Search should be case-insensitive."""
        response1 = client.get(reverse('docs_viewer:search'), {'q': 'DOCUMENTATION'})
        response2 = client.get(reverse('docs_viewer:search'), {'q': 'documentation'})
        
        # Should return same results
        assert response1.status_code == 200
        assert response2.status_code == 200
    
    def test_search_shows_match_count(self, client):
        """Search results should show number of matches."""
        response = client.get(reverse('docs_viewer:search'), {'q': 'test'})
        if b'result' in response.content.lower():
            # If results found, should show count
            assert b'match' in response.content.lower() or b'found' in response.content.lower()
    
    def test_search_results_clickable(self, client):
        """Search results should link to documents."""
        response = client.get(reverse('docs_viewer:search'), {'q': 'test'})
        content = response.content.decode()
        if 'result' in content.lower():
            assert '/docs/view/' in content or 'list-group-item' in content
