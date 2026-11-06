"""
Views for documentation viewer with security controls.
"""
import os
from django.shortcuts import render, redirect
from django.http import HttpResponse
from django.conf import settings
from django.core.paginator import Paginator
from django.contrib import messages
from django.views.decorators.cache import never_cache
from django.views.decorators.http import require_http_methods

from .utils import (
    is_safe_path,
    get_allowed_doc_paths,
    render_markdown,
    extract_title,
    get_breadcrumbs,
    sanitize_path,
    validate_markdown_links,
    document_metadata,
    document_language,
)
from .models import DocumentView


@never_cache  # Navigation contains the current account and CSRF-protected logout.
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
        
        try:
            with open(os.path.join(settings.BASE_DIR, doc_path), encoding='utf-8') as source:
                label = extract_title(source.read())
        except OSError:
            continue
        organized[directory].append({
            'path': doc_path,
            'name': label if label != 'Untitled Document' else os.path.basename(doc_path),
            'url': f'/docs/view/{doc_path}',
        })
    
    context = {
        'organized_docs': organized,
        'total_docs': len(docs),
    }
    
    return render(request, 'docs_viewer/index.html', context)


@never_cache
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
        return render(request, 'docs_viewer/not_found.html', status=404)
    
    # Read file
    full_path = os.path.join(settings.BASE_DIR, path)
    try:
        with open(full_path, 'r', encoding='utf-8') as f:
            content = f.read()
    except OSError:
        return render(request, 'docs_viewer/not_found.html', status=404)
    
    # Render markdown
    html_content, toc_html = render_markdown(content, source_path=path)
    
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
        'document_metadata': document_metadata(content),
        'document_language': document_language(content, path),
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
        return render(request, 'docs_viewer/not_found.html', status=404)
    
    # Read file
    full_path = os.path.join(settings.BASE_DIR, path)
    try:
        with open(full_path, 'r', encoding='utf-8') as f:
            content = f.read()
    except OSError:
        return render(request, 'docs_viewer/not_found.html', status=404)
    
    # Serve as plain text
    response = HttpResponse(content, content_type='text/plain; charset=utf-8')
    response['Content-Disposition'] = f'attachment; filename="{os.path.basename(path)}"'
    
    return response


@never_cache
@require_http_methods(["GET"])
def docs_search(request):
    """Search documentation content."""
    query = request.GET.get('q', '').strip()[:200]
    
    if not query:
        return redirect('docs_viewer:index')
    
    results = []
    docs = get_allowed_doc_paths()
    
    for doc_path in docs:
        if not is_safe_path(doc_path):
            continue
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
                    'name': extract_title(content) if extract_title(content) != 'Untitled Document' else os.path.basename(doc_path),
                    'url': f'/docs/view/{doc_path}',
                    'matches': len(matching_lines),
                    'preview': matching_lines[0][1][:200] if matching_lines else '',
                })
        except Exception:
            continue
    
    # Sort by number of matches
    results.sort(key=lambda x: x['matches'], reverse=True)
    
    page = Paginator(results, 20).get_page(request.GET.get('page'))
    context = {
        'query': query,
        'results': page.object_list,
        'page': page,
        'total_results': len(results),
    }
    
    return render(request, 'docs_viewer/search.html', context)
