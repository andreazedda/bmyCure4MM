---
title: Endpoints and Authorization Map
language: en
status: CURRENT_VERIFIED
owner: Andrea Zedda
audience: developers and maintainers
last_verified_at: 2026-11-06
---

# Endpoints and authorization map

Paths below were checked against the actual URLConfs during the 2026-11-06 audit. `<id>` is a route parameter, not a public authorization capability. Existence of a route does not grant access to its data. The [source atlas](../development/CODE_ATLAS.md) follows the implementation and tests behind these entry points.

## Access and first use / Accesso e primo utilizzo

| Path | Purpose | Implementation / access boundary |
| --- | --- | --- |
| `/` | Public entry | `mmportal.views.home` |
| `/access/request/` | Email/nickname request, no password | `accounts.views.request_access`; request gate and rate budget |
| `/access/request/status/` | Pending request status | Session-limited; knowing an email or UUID is insufficient |
| `/access/verify/` | Email confirmation | GET displays form; purpose token is consumed only by POST |
| `/access/activate/` | Approved member password setup | POST rechecks identity, approval, expiry and revocation |
| `/accounts/login/`, `/accounts/logout/` | Ordinary product authentication | Email backend; logout requires POST |
| `/accounts/security/` | Owner authenticator setup / fresh code | Privileged product session, separate from staff status |
| `/password/reset/`, `/password/reset/confirm/` | Password recovery | Does not reactivate revoked/suspended/expired access |
| `/manage/access-requests/`, `/manage/users/` | Review and account lifecycle | Owner/authorized approver plus fresh MFA |
| `/manage/feedback/` | Private ticket replies | Same privileged product boundary |
| `/feedback/` | Member's own tickets | Active community identity; author-scoped output |
| `/sim/getting-started/` | Guided existing learning entry | Product member access |
| `/learn/`, `/sim/` | Existing scenario catalogue | Same simulator URLConf; permissions still apply |
| `/admin/` | Separate Django administration | Existing staff requirements; public registration never grants staff |

## Existing simulator / Simulatore esistente

| Path | View / behavior |
| --- | --- |
| `/sim/<id>/` | `simulator.views.scenario_detail` |
| `/sim/scenarios/<id>/simulate/` | `views_manage.simulate_scenario`, validated POST, existing real engine |
| `/sim/attempts/<id>/export/` | `views_manage.export_attempt`; current author/access check; `?format=csv` for CSV |
| `/sim/attempts/<id>/plot/embed/` | `views_manage.attempt_plot_embed`; authenticated HTML artifact |
| `/sim/manage/`, `/sim/scenarios/new/` | Existing editor-only scenario management, not community registration |
| `/sim/cohort/`, `/sim/experiments/` | Existing advanced research functions; not granted by approved registration |

`/simulator/` and `/learn/` include the same simulator routes under distinct namespaces. Tests must cover aliases; hiding a link is not authorization. The old reference path `/simulator/scenario/<id>/` was incorrect: the detail route is directly `<id>/` under the chosen prefix.

## Clinical, research and molecular surfaces

These are documented for source navigation, not offered as public onboarding steps.

| Path | Source | Boundary |
| --- | --- | --- |
| `/patients/`, `/patients/new/`, `/patients/<id>/` | `clinic/urls.py`, `clinic/views.py` | Existing clinical permissions; registered community identities denied |
| `/api/patients/`, `/api/assessments/` | `clinic/api.py` | DRF permissions and scoped querysets; account authentication alone is insufficient |
| `/research/` and patient research routes | `twin_engine/urls.py`, `twin_engine/views.py` | Separate governed research access; does not follow from access approval |
| `/chemtools/` | `chemtools/urls.py`, `chemtools/views.py` | Existing molecular utility workflow |
| `/chemtools/job/<id>/` | `chemtools.views.job_detail` | Job selector/authorization; optional dependencies and broker remain separate |
| `/chemtools/job/<id>/status.json` | `chemtools.views.job_status` | Polling endpoint, not proof of worker progress |

`/chem/` is a backward-compatible alias. Both use singular `job/<id>/`; the old `/chem/jobs/<id>/` reference was incorrect. Keep aliases in denial tests and do not expose media to make a job link work.

## Documentation and runtime

| Path | Source / behavior |
| --- | --- |
| `/docs/` | Documentation entry, `docs_viewer.views.docs_index` |
| `/docs/view/docs/development/DEVELOPMENT.md` | Allowed Markdown path; no trailing slash after the filename |
| `/docs/search/?q=...` | Allowed-document search; query should contain no secrets or clinical content |
| `/docs/raw/<path>` | Allowed raw Markdown only |
| `/healthz/` | Minimal web liveness response; not a DB/worker/email health guarantee |
| `/static/...` | Existing static service; check MIME type and actual loaded content |
| `/media/...` | Direct protected media is closed in production; use authorized artifact views |

## Verify routing without loading a clinical record

Use the disposable settings from the [handbook](../development/DEVELOPMENT.md):

```bash
uv run --frozen python manage.py shell -c 'from django.urls import resolve; print(resolve("/sim/1/").view_name); print(resolve("/sim/attempts/1/export/").view_name); print(resolve("/chemtools/job/1/").view_name)'
```

Expected names: `simulator:scenario_detail`, `simulator:attempt_export`, `chemtools:job_detail`. Route resolution alone does not query those IDs or prove access. Run account-denial and artifact-ownership tests for authorization evidence.
