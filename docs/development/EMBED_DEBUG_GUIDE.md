# 🔍 Iframe Embedding & Debug Guide

> Troubleshooting guide for interactive plot rendering in iframe contexts

## Overview

This guide covers debugging and resolving iframe rendering issues, specifically focusing on Plotly visualization embedding within the bmyCure4MM application. It documents the solutions implemented to overcome CDN-based JavaScript loading issues in iframe contexts.

## The Problem: Grey Iframe

### Symptoms
- Iframe appears as grey/blank box on patient detail page
- Same plot opens correctly in new tab/window
- No visible JavaScript errors in console
- X-Frame-Options headers correctly configured

### Root Cause
CDN-based Plotly.js (`https://cdn.plot.ly/plotly-3.1.1.min.js`) fails to load in iframe context due to:
1. Content Security Policy restrictions
2. Cross-origin resource sharing (CORS) limitations
3. Browser security policies blocking external scripts in iframes

## Solution Architecture

### 1. Embed-Safe Endpoint

**Location:** `simulator/views_manage.py`

```python
@login_required
def attempt_plot_embed(request: HttpRequest, pk: int) -> HttpResponse:
    """
    Embed-safe plot endpoint with inline Plotly.js.
    
    Why: CDN-based Plotly fails in iframe. This endpoint injects inline JS.
    """
    attempt = get_object_or_404(SimulationAttempt, pk=pk, scenario__patient=request.user)
    
    # Path validation for security
    plot_path = attempt.plot_file.path if attempt.plot_file else None
    if not plot_path or not os.path.exists(plot_path):
        return HttpResponse("Plot not available", status=404)
    
    # Security: validate path is within MEDIA_ROOT
    if not plot_path.startswith(settings.MEDIA_ROOT):
        return HttpResponse("Invalid path", status=403)
    
    # Read original plot HTML
    with open(plot_path, "r", encoding="utf-8") as f:
        html = f.read()
    
    # Replace CDN script with inline Plotly.js
    plotly_js = plotly.io.get_plotlyjs()  # ~3MB minified
    html = html.replace(
        '<script src="https://cdn.plot.ly/plotly-3.1.1.min.js"></script>',
        f'<script>{plotly_js}</script>'
    )
    
    return HttpResponse(html, content_type="text/html", headers={
        "Cache-Control": "no-store",
        "X-Frame-Options": "SAMEORIGIN"
    })
```

**Key Features:**
- ✅ No external CDN dependencies
- ✅ Path validation prevents directory traversal
- ✅ Inline Plotly.js (~3MB) embedded directly
- ✅ Cache-Control prevents stale data
- ✅ X-Frame-Options allows same-origin embedding

### 2. URL Routing

**Location:** `simulator/urls.py`

```python
urlpatterns = [
    # ... existing patterns
    path(
        "attempts/<int:pk>/plot/embed/",
        views_manage.attempt_plot_embed,
        name="attempt_plot_embed"
    ),
]
```

### 3. Template Integration

**Location:** `clinic/templates/clinic/patient_detail.html`

```html
<!-- OLD: CDN-based plot (fails in iframe) -->
<iframe src="{% url 'simulator:attempt_plot_detail' latest_simulation.pk %}"></iframe>

<!-- NEW: Embed-safe plot (works in iframe) -->
<iframe 
    src="{% url 'simulator:attempt_plot_embed' latest_simulation.pk %}" 
    style="width: 100%; height: 600px; border: 1px solid #dee2e6; border-radius: 0.375rem;"
    loading="lazy">
</iframe>
```

## Debug Logging System

### EmbedDebugMiddleware

**Location:** `mmportal/middleware.py`

```python
import logging
import time

logger = logging.getLogger("embed_debug")

class EmbedDebugMiddleware:
    """Log embed requests for debugging iframe issues."""
    
    def __init__(self, get_response):
        self.get_response = get_response
    
    def __call__(self, request):
        should_log = self._should_log(request)
        
        if should_log:
            self._log("=== EMBED REQUEST START ===")
            self._log(f"Method: {request.method}")
            self._log(f"Path: {request.path}")
            self._log(f"Headers: {dict(request.headers)}")
            start_time = time.time()
        
        response = self.get_response(request)
        
        if should_log:
            duration = time.time() - start_time
            self._log(f"Status: {response.status_code}")
            self._log(f"Response Headers: {dict(response.headers)}")
            self._log(f"Duration: {duration:.3f}s")
            self._log(f"Content-Length: {len(response.content) if hasattr(response, 'content') else 'N/A'}")
            self._log("=== EMBED REQUEST END ===\n")
        
        return response
    
    def _should_log(self, request) -> bool:
        """Only log embed-related requests."""
        embed_paths = ["/sim/attempts/", "/plot/embed/"]
        return any(p in request.path for p in embed_paths)
    
    def _log(self, message: str):
        logger.info(message)
```

**Configuration:** `mmportal/settings.py`

```python
MIDDLEWARE = [
    # ... other middleware
    "mmportal.middleware.EmbedDebugMiddleware",  # Add near end
]

LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'handlers': {
        'embed_debug_file': {
            'level': 'INFO',
            'class': 'logging.FileHandler',
            'filename': BASE_DIR / 'logs' / 'embed_debug.log',
            'encoding': 'utf-8',
        },
    },
    'loggers': {
        'embed_debug': {
            'handlers': ['embed_debug_file'],
            'level': 'INFO',
            'propagate': False,
        },
    },
}
```

### Log File Location

```bash
/Volumes/nvme/Github/bmyCure4MM/logs/embed_debug.log
```

### Reading Logs

```bash
# View recent embed requests
tail -n 50 logs/embed_debug.log

# Watch logs in real-time
tail -f logs/embed_debug.log

# Filter for specific attempt
grep "attempts/42" logs/embed_debug.log

# Check for errors
grep -i "error\|404\|500" logs/embed_debug.log
```

### Sample Log Output

```
=== EMBED REQUEST START ===
Method: GET
Path: /sim/attempts/17/plot/embed/
Headers: {'Host': 'localhost:8000', 'User-Agent': 'Mozilla/5.0...', 'Referer': 'http://localhost:8000/clinic/patients/1/'}
Status: 200
Response Headers: {'Content-Type': 'text/html; charset=utf-8', 'X-Frame-Options': 'SAMEORIGIN', 'Cache-Control': 'no-store'}
Duration: 0.043s
Content-Length: 32648
=== EMBED REQUEST END ===
```

## Troubleshooting Checklist

### Issue: Iframe Still Shows Grey Box

**Check 1: Verify Endpoint**
```bash
# Test embed endpoint directly
curl -I http://localhost:8000/sim/attempts/17/plot/embed/
# Expected: HTTP 200, X-Frame-Options: SAMEORIGIN
```

**Check 2: Verify Plot File Exists**
```python
from simulator.models import SimulationAttempt
attempt = SimulationAttempt.objects.get(pk=17)
print(attempt.plot_file.path)  # Should be valid path
print(os.path.exists(attempt.plot_file.path))  # Should be True
```

**Check 3: Check Logs**
```bash
tail -n 20 logs/embed_debug.log
# Look for: Status 200, Content-Length > 30000 (with inline Plotly)
```

**Check 4: Browser Console**
```javascript
// Open browser console (F12) on patient page
// Check for errors like:
// - "Refused to display in a frame because it set 'X-Frame-Options' to 'deny'"
// - "Failed to load resource"
// - CSP violations
```

### Issue: X-Frame-Options Error

**Symptom:** Console error "Refused to display... X-Frame-Options"

**Solution:** Check settings.py
```python
# Should be SAMEORIGIN (not DENY)
X_FRAME_OPTIONS = 'SAMEORIGIN'
```

**Verify in response:**
```bash
curl -I http://localhost:8000/sim/attempts/17/plot/embed/ | grep X-Frame-Options
# Expected: X-Frame-Options: SAMEORIGIN
```

### Issue: Path Validation Error (403)

**Symptom:** HTTP 403 response, message "Invalid path"

**Cause:** Plot file path outside MEDIA_ROOT (security check)

**Debug:**
```python
from django.conf import settings
from simulator.models import SimulationAttempt

attempt = SimulationAttempt.objects.get(pk=17)
plot_path = attempt.plot_file.path
media_root = settings.MEDIA_ROOT

print(f"Plot path: {plot_path}")
print(f"MEDIA_ROOT: {media_root}")
print(f"Valid: {plot_path.startswith(media_root)}")
```

**Solution:** Ensure plot files stored in MEDIA_ROOT subdirectory

### Issue: Content-Length Too Small

**Symptom:** Log shows Content-Length < 10KB (should be ~32KB with inline Plotly)

**Cause:** Plotly.js not being injected (string replacement failed)

**Debug:**
```python
import plotly.io

# Verify get_plotlyjs() works
plotly_js = plotly.io.get_plotlyjs()
print(f"Plotly.js size: {len(plotly_js)} bytes")  # Should be ~3MB

# Check if CDN script exists in original plot
with open('/path/to/plot.html', 'r') as f:
    html = f.read()
    print('<script src="https://cdn.plot.ly' in html)  # Should be True
```

**Solution:** Verify CDN URL matches exactly in replacement

### Issue: Plot Doesn't Update

**Symptom:** Old plot still showing after new simulation

**Cause:** Browser cache (despite Cache-Control header)

**Solutions:**
1. Hard refresh: Ctrl+Shift+R (Win/Linux) or Cmd+Shift+R (Mac)
2. Clear browser cache
3. Add cache-busting query param:
```html
<iframe src="{% url 'simulator:attempt_plot_embed' latest_simulation.pk %}?t={{ latest_simulation.updated_at|date:'U' }}"></iframe>
```

## Security Considerations

### Path Validation

**Why:** Prevent directory traversal attacks

```python
# Security: validate path is within MEDIA_ROOT
if not plot_path.startswith(settings.MEDIA_ROOT):
    return HttpResponse("Invalid path", status=403)
```

**Attack Example (prevented):**
```
/sim/attempts/17/plot/embed/
  → plot_path = "../../../../etc/passwd"  # Blocked by validation
```

### X-Frame-Options

**Purpose:** Control which domains can embed your content

**Options:**
- `DENY`: No embedding allowed (blocks all iframes)
- `SAMEORIGIN`: Only same domain can embed (✅ recommended)
- `ALLOW-FROM uri`: Specific domain can embed (deprecated)

**Best Practice:** Use `SAMEORIGIN` for internal iframes

### Content Security Policy (CSP)

If you add CSP headers, ensure iframes allowed:

```python
# In middleware or view
response['Content-Security-Policy'] = "frame-ancestors 'self'"
```

### Authentication

Embed endpoint uses `@login_required` and verifies patient ownership:

```python
attempt = get_object_or_404(
    SimulationAttempt, 
    pk=pk, 
    scenario__patient=request.user  # Security: only owner's data
)
```

## Performance Optimization

### Inline Plotly.js Trade-off

**Pros:**
- ✅ Works in iframe (solves CDN issue)
- ✅ No external dependency
- ✅ Offline-capable

**Cons:**
- ❌ ~3MB per response (vs ~400B CDN link)
- ❌ No browser cache across pages
- ❌ Slower initial load

### Mitigation Strategies

**1. Lazy Loading:**
```html
<iframe src="..." loading="lazy"></iframe>
```
Delays iframe load until user scrolls near it.

**2. Conditional Inline:**
```python
# Only inline if CDN fails
if request.GET.get('fallback') == '1':
    # Use inline Plotly.js
else:
    # Try CDN first
```

**3. Service Worker Caching:**
Cache inline Plotly.js in service worker for repeat visits.

**4. Compression:**
Enable gzip compression in web server (reduces ~3MB to ~1MB).

## Testing Procedures

### Manual Testing

1. **Basic Rendering:**
   - Navigate to patient page with simulation
   - Verify iframe shows plot (not grey box)
   - Check plot is interactive (zoom, pan, hover)

2. **New Simulation:**
   - Run new simulation from simulator
   - Return to patient page
   - Verify iframe updates to new plot

3. **Browser Compatibility:**
   - Test in Chrome, Firefox, Safari, Edge
   - Check mobile browsers (iOS Safari, Android Chrome)

4. **Network Conditions:**
   - Test with slow 3G throttling
   - Verify lazy loading works
   - Check load time acceptable

### Automated Testing

```python
from django.test import TestCase, Client
from simulator.models import SimulationAttempt

class EmbedPlotTest(TestCase):
    def test_embed_endpoint_returns_html(self):
        attempt = SimulationAttempt.objects.create(...)
        response = self.client.get(f'/sim/attempts/{attempt.pk}/plot/embed/')
        
        self.assertEqual(response.status_code, 200)
        self.assertIn('text/html', response['Content-Type'])
        self.assertIn('SAMEORIGIN', response['X-Frame-Options'])
        self.assertGreater(len(response.content), 30000)  # Has inline Plotly
    
    def test_path_validation(self):
        # Test that invalid paths are rejected
        # ...
    
    def test_authentication_required(self):
        # Test that unauthenticated requests fail
        # ...
```

## Common Error Messages

| Error | Cause | Solution |
|-------|-------|----------|
| "Refused to display in frame" | X-Frame-Options: DENY | Set to SAMEORIGIN |
| "Plot not available" (404) | Plot file missing | Re-run simulation |
| "Invalid path" (403) | Path outside MEDIA_ROOT | Check file storage configuration |
| Grey box, no errors | CDN Plotly blocked | Use embed endpoint (not CDN endpoint) |
| "Permission denied" | User not simulation owner | Verify authentication |
| Stale plot showing | Browser cache | Hard refresh or add cache-busting param |

## Monitoring

### Log Rotation

Prevent embed_debug.log from growing indefinitely:

```bash
# Add to crontab
0 0 * * 0 gzip /path/to/logs/embed_debug.log && > /path/to/logs/embed_debug.log
```

### Metrics to Track

1. **Embed request latency**: Should be < 100ms
2. **404 rate**: Should be < 1% (indicates missing plots)
3. **403 rate**: Should be 0% (indicates path validation issues)
4. **Content-Length**: Should be ~32KB (indicates inline Plotly working)

## Future Improvements

- [ ] Implement service worker caching for Plotly.js
- [ ] Add progressive enhancement (CDN with inline fallback)
- [ ] Create admin dashboard for iframe metrics
- [ ] Implement automated plot archival/cleanup
- [ ] Add WebSocket for real-time plot updates
- [ ] Support for other visualization libraries (D3.js, Chart.js)

## References

- [Django X-Frame-Options Documentation](https://docs.djangoproject.com/en/4.2/ref/clickjacking/)
- [Plotly.js Documentation](https://plotly.com/javascript/)
- [MDN: CSP frame-ancestors](https://developer.mozilla.org/en-US/docs/Web/HTTP/Headers/Content-Security-Policy/frame-ancestors)
- [OWASP: Clickjacking Defense](https://cheatsheetseries.owasp.org/cheatsheets/Clickjacking_Defense_Cheat_Sheet.html)

---

**Last Updated:** January 2026  
**Maintainer:** Development Team  
**Related:** [Decision Support System](../features/DECISION_SUPPORT_SYSTEM.md)
