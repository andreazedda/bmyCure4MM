"""End-to-end smoke test against the running dev server."""
import requests
from html.parser import HTMLParser

BASE = "http://127.0.0.1:8000"
session = requests.Session()

# 1) Anon should be redirected
routes_anon = ["/", "/patients/", "/patients/new/"]
for r in routes_anon:
    resp = session.get(BASE + r, allow_redirects=False)
    assert resp.status_code in (301, 302), f"FAIL: {r} returned {resp.status_code}"
    print(f"  [PASS] Anon {r} -> redirect {resp.status_code}")

# 2) Login
login_page = session.get(BASE + "/admin/login/")

class CSRFParser(HTMLParser):
    csrf = ""
    def handle_starttag(self, tag, attrs):
        d = dict(attrs)
        if d.get("name") == "csrfmiddlewaretoken":
            self.csrf = d.get("value", "")

p = CSRFParser()
p.feed(login_page.text)
resp = session.post(
    BASE + "/admin/login/",
    data={
        "csrfmiddlewaretoken": p.csrf,
        "username": "admin",
        "password": "admin",
        "next": "/",
    },
    headers={"Referer": BASE + "/admin/login/"},
)
assert resp.status_code in (200, 302), f"Login failed: {resp.status_code}"
print("  [PASS] Admin login OK")

# 3) Authenticated routes
authed_routes = [
    ("/", 200),
    ("/patients/", 200),
    ("/patients/new/", 200),
    ("/regimens/", 200),
    ("/regimens/new/", 200),
    ("/chem/", 200),
    ("/docs/", 200),
]
for path, expected in authed_routes:
    resp = session.get(BASE + path)
    assert resp.status_code == expected, f"FAIL: {path} returned {resp.status_code}"
    print(f"  [PASS] {path} -> {resp.status_code}")

print()
print("All end-to-end checks passed!")
