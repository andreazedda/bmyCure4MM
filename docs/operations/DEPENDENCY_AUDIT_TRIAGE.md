# Dependency Audit Triage

Audit date: 2026-08-13

The broad legacy environment had 161 findings in 25 packages. The canonical
locked environment removes 155 findings. `pip-audit` reports six remaining
Django advisories; all are low or medium severity and none has a Django 4.2
patch. Issue #15 requires one consistent Django 4.2 series, so Django is pinned
to the final 4.2.30 release and each exception is narrowly allowlisted by ID.

| Advisory | Severity | Affected surface | Repository evidence / mitigation |
| --- | --- | --- | --- |
| [GHSA-923m-gv2p-w5qp](https://github.com/advisories/GHSA-923m-gv2p-w5qp) | low | `has_vary_header` whitespace caching | No project call and no Django cache middleware. |
| [GHSA-h7pc-vwp9-298g](https://github.com/advisories/GHSA-h7pc-vwp9-298g) | low | signed-cookie salt collision | Database sessions are used; no signed-cookie session backend. |
| [GHSA-8cjm-8mp7-r2xf](https://github.com/advisories/GHSA-8cjm-8mp7-r2xf) | low | `UpdateCacheMiddleware` header handling | Cache update middleware is absent. |
| [GHSA-3h9f-r86x-qvjx](https://github.com/advisories/GHSA-3h9f-r86x-qvjx) | low | cache middleware private responses | Django cache middleware is absent. |
| [GHSA-crhf-3pfg-w68w](https://github.com/advisories/GHSA-crhf-3pfg-w68w) | medium | `GDALRaster` byte over-read | GeoDjango/GDAL are neither declared nor imported. |
| [GHSA-8qcx-xf44-272x](https://github.com/advisories/GHSA-8qcx-xf44-272x) | medium | `DomainNameValidator` newline acceptance | The validator is not imported; host configuration is parsed and validated separately. |

One additional medium development-only advisory,
[GHSA-6w46-j5rx-g56g](https://github.com/advisories/GHSA-6w46-j5rx-g56g),
affects pytest temporary-directory handling on shared UNIX hosts. The compatible
`pytest-playwright==0.7.1` release requires pytest below 9. Tests run on isolated
CI runners or single-user development machines; untrusted users must not share
that host. Remove the exception as soon as pytest-playwright supports patched
pytest 9.

There are no ignored critical or high findings. `scripts/audit_dependencies.py`
ignores only these exact identifiers, so any new advisory fails CI. A future
Django 5.2 migration must remove the Django exceptions and rerun authorization,
migration, template, numerical, and output-language regression suites.
