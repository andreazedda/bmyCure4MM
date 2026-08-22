"""Manual end-to-end smoke check against a running development server."""

from html.parser import HTMLParser

import requests

BASE = "http://127.0.0.1:8000"

class CSRFParser(HTMLParser):
    csrf = ""
    def handle_starttag(self, tag, attrs):
        d = dict(attrs)
        if d.get("name") == "csrfmiddlewaretoken":
            self.csrf = d.get("value", "")


def main() -> None:
    session = requests.Session()

    routes_anon = ["/", "/patients/", "/patients/new/"]
    for route in routes_anon:
        response = session.get(BASE + route, allow_redirects=False)
        assert response.status_code in (301, 302), (
            f"FAIL: {route} returned {response.status_code}"
        )
        print(f"  [PASS] Anon {route} -> redirect {response.status_code}")

    login_page = session.get(BASE + "/admin/login/")
    parser = CSRFParser()
    parser.feed(login_page.text)
    response = session.post(
        BASE + "/admin/login/",
        data={
            "csrfmiddlewaretoken": parser.csrf,
            "username": "admin",
            "password": "admin",
            "next": "/",
        },
        headers={"Referer": BASE + "/admin/login/"},
    )
    assert response.status_code in (200, 302), f"Login failed: {response.status_code}"
    print("  [PASS] Admin login OK")

    authed_routes = [
        ("/", 200),
        ("/patients/", 200),
        ("/patients/new/", 200),
        ("/regimens/", 200),
        ("/regimens/new/", 200),
        ("/chem/", 200),
        ("/docs/", 200),
    ]
    for route, expected in authed_routes:
        response = session.get(BASE + route)
        assert response.status_code == expected, (
            f"FAIL: {route} returned {response.status_code}"
        )
        print(f"  [PASS] {route} -> {response.status_code}")

    print("\nAll end-to-end checks passed!")


if __name__ == "__main__":
    main()
