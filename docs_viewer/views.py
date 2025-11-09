"""
Views for documentation viewer with security controls.
"""
import os
from django.shortcuts import render, redirect
from django.http import Http404, HttpResponse
from django.conf import settings
from django.contrib import messages
from django.views.decorators.cache import cache_page
from django.views.decorators.http import require_http_methods

from .utils import (
    is_safe_path,
    get_allowed_doc_paths,
    render_markdown,
    extract_title,
    get_breadcrumbs,
    sanitize_path,
    validate_markdown_links,
)
from .models import DocumentView


@cache_page(60 * 15)  # Cache for 15 minutes
def docs_index(request):
    """List all available documentation files."""
    docs = get_allowed_doc_paths()
    
    # Organize by directory
    organized = {}
    for doc_path in docs:
        if '/' in doc_path:
            directory = doc_path.split('/')[0]
        else:
            directory = 'Root'
        
        if directory not in organized:
            organized[directory] = []
        
        organized[directory].append({
            'path': doc_path,
            'name': os.path.basename(doc_path),
            'url': f'/docs/view/{doc_path}',
        })
    
    context = {
        'organized_docs': organized,
        'total_docs': len(docs),
    }
    
    return render(request, 'docs_viewer/index.html', context)


def docs_view(request, path):
    """
    View a documentation file with security checks.
    
    Security:
    - Path sanitization
    - Whitelist validation
    - Path traversal prevention
    """
    # Sanitize input
    path = sanitize_path(path)
    
    # Security check
    if not is_safe_path(path):
        raise Http404("Documentation file not found or access denied")
    
    # Read file
    full_path = os.path.join(settings.BASE_DIR, path)
    try:
        with open(full_path, 'r', encoding='utf-8') as f:
            content = f.read()
    except Exception as e:
        raise Http404(f"Error reading file: {str(e)}")
    
    # Render markdown
    html_content, toc_html = render_markdown(content)
    
    # Extract metadata
    title = extract_title(content)
    breadcrumbs = get_breadcrumbs(path)
    
    # Validate links (only in DEBUG mode to avoid performance hit)
    broken_links = []
    if settings.DEBUG:
        broken_links = validate_markdown_links(content, path)
        if broken_links:
            messages.warning(
                request,
                f"Found {len(broken_links)} broken links in this document"
            )
    
    # Track view
    if request.user.is_authenticated:
        DocumentView.objects.create(
            path=path,
            user=request.user
        )
    
    context = {
        'title': title,
        'content': html_content,
        'toc': toc_html,
        'breadcrumbs': breadcrumbs,
        'path': path,
        'raw_url': f'/docs/raw/{path}',
        'broken_links': broken_links,
    }
    
    return render(request, 'docs_viewer/view.html', context)


def docs_raw(request, path):
    """
    Serve raw markdown file for download.
    
    Security: Same checks as docs_view
    """
    # Sanitize input
    path = sanitize_path(path)
    
    # Security check
    if not is_safe_path(path):
        raise Http404("Documentation file not found or access denied")
    
    # Read file
    full_path = os.path.join(settings.BASE_DIR, path)
    try:
        with open(full_path, 'r', encoding='utf-8') as f:
            content = f.read()
    except Exception as e:
        raise Http404(f"Error reading file: {str(e)}")
    
    # Serve as plain text
    response = HttpResponse(content, content_type='text/plain; charset=utf-8')
    response['Content-Disposition'] = f'attachment; filename="{os.path.basename(path)}"'
    
    return response


@require_http_methods(["GET"])
def docs_search(request):
    """Search documentation content."""
    query = request.GET.get('q', '').strip()
    
    if not query:
        return redirect('docs_viewer:index')
    
    results = []
    docs = get_allowed_doc_paths()
    
    for doc_path in docs:
        full_path = os.path.join(settings.BASE_DIR, doc_path)
        
        try:
            with open(full_path, 'r', encoding='utf-8') as f:
                content = f.read()
            
            # Case-insensitive search
            if query.lower() in content.lower():
                # Extract context around match
                lines = content.split('\n')
                matching_lines = [
                    (i, line) for i, line in enumerate(lines)
                    if query.lower() in line.lower()
                ]
                
                results.append({
                    'path': doc_path,
                    'name': os.path.basename(doc_path),
                    'url': f'/docs/view/{doc_path}',
                    'matches': len(matching_lines),
                    'preview': matching_lines[0][1][:200] if matching_lines else '',
                })
        except Exception:
            continue
    
    # Sort by number of matches
    results.sort(key=lambda x: x['matches'], reverse=True)
    
    context = {
        'query': query,
        'results': results,
        'total_results': len(results),
    }
    
    return render(request, 'docs_viewer/search.html', context)
