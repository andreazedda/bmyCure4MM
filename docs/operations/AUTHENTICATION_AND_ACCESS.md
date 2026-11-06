---
title: Authentication and Access
status: CURRENT_PARTIAL
owner: Andrea Zedda
audience: users, maintainers and operators
last_verified_at: 2026-11-06
---

# Authentication and access

## Scope and verified production boundary

Registration and private suggestion tickets extend the existing Django application and user model. They use the existing navigation, getting-started guide and learning scenarios. There is no additional demo portal, scenario provisioning, simulation worker or parallel quota system. Scientific equations, parameters, solver, manifests and intended use remain unchanged: E1 research only, never clinical decision support.

The verified production origin is **https://bmycure4mm.clusterlab.uk**, namespace `bmycure4mm`, deployments `web` and `celery`. On 2026-11-06, before this release, both ran source `2c60b63dfb3d42ee573d749c23cf5033aabdbc68` and image digest `sha256:a87b1118813ebf2c1fea827dc9d5038595b4abe2420ba362acc87865f7606626`. Product login was absent. Runtime uses SQLite `/data/db.sqlite3` on the existing NFS PVC, not PostgreSQL. Static files already have a separate service and PVC; adding middleware is unnecessary.

Implementation and local tests do not prove a production rollout or email receipt. See the release evidence below for the current verified state. Broad clinical authorization work (#8), resource governance and the real researcher pilot (#78) remain separate and open.

## Short user guide

1. Open `/access/request/`. Supply email, nickname and a short reason, optional affiliation/role, and acknowledge the privacy/research limits. Never enter patient data, credentials or attachments.
2. Open the verification email and press **Confirm email**. Merely opening the link does not verify or approve anything. This browser can see status at `/access/request/status/`; email or a sequential ID cannot retrieve it.
3. Wait for Andrea's decision. After approval, the activation email lets you set your own password (at least 12 characters). No password is requested before approval or sent by email.
4. Sign in at `/accounts/login/`. The existing `/sim/getting-started/` guide links to learning scenarios. Registration grants no patient, import, advanced research, management or Django admin access. Registered users see only their own saved simulation attempts and exports.
5. Open `/feedback/` to send a suggestion, bug report, usability issue or scientific/interpretation question. The ticket and Andrea's reply persist across logout/login. Internal notes and other users' tickets are hidden.
6. Use the account menu to **Log out**. Password reset at `/password/reset/` cannot reinstate suspended, expired or revoked access.

Access initially expires after the 1–14 days assigned by Andrea. Existing simulator validation still applies; this change does not introduce or claim a new public computation budget. Do not open an unrestricted computational pilot on the strength of registration tests.

## Andrea's product console

Use `/accounts/login/`, then `/accounts/security/` to enroll/confirm a time-based authenticator. Management sessions require a fresh code every 15 minutes. Keep the authenticator key private and retain recovery access through the governed cluster operator. Existing administrators are preserved; new registrants remain `is_staff=false`, `is_superuser=false`.

- `/manage/access-requests/`: verified requests, approve/reject, choose initial expiry, reissue pending activation, renew/restore suspended or expired access, suspend/revoke, and failed email retry. Revoked access cannot be renewed.
- `/manage/users/`: registered users and expiry/status; link back to access controls.
- `/manage/feedback/`: classify, reply and close private tickets. Publishing a reply sends a generic notification without the ticket body.
- Registration controls independently pause requests and activations. Owner access remains available while these gates are closed.

Initially only the configured owner is an approver. Public role/owner/staff/quota/scope fields are rejected. There is no implicit staff promotion or delegation. Django admin remains a separate surface at `/admin/` and owner admin sessions also require authenticator confirmation.

## Lifecycle and identity contract

| Concept | Transitions |
| --- | --- |
| Request | `PENDING_EMAIL_VERIFICATION` → `PENDING_REVIEW` → `APPROVED` / `REJECTED` |
| Access | `PENDING_ACTIVATION` → `ACTIVE` → `SUSPENDED` / `EXPIRED` / `REVOKED` |
| Notification | `queued` → `sent` / `failed`, independently of access |

Email identity is trimmed/lowercased ASCII; dots and plus suffixes are preserved. An additive unique expression index also protects nonempty legacy `auth_user` email identities. Migration refuses duplicates without merging users or printing addresses. Nickname is a display name only.

Lifecycle mutations serialize through a database singleton write lock, including on the deployed SQLite backend; PostgreSQL `select_for_update` is not assumed to work on SQLite. Approval, activation, ticket submissions and outbox creation are transactional/idempotent. Activation validates identity, purpose, approval, expiry, role and the activation gate. Access state and session version are checked on every authenticated request. Expiry/revocation prevent subsequent access to protected pages and exports.

Capabilities use random secrets from maintained standard-library primitives, SHA-256 digests for validation, distinct purposes, one-use POST confirmation, expiry and revocation. Verification lasts 24 hours; activation/reset lasts two hours. Links use only configured `PUBLIC_BASE_URL`; tokens are URL fragments, immediately removed from browser history and never sent in HTTP paths/referrers. Invalid forms never echo tokens. Owner MFA uses `pyotp` and encrypted secrets.

## Runtime configuration and email

Supply `PUBLIC_BASE_URL` with the verified HTTPS origin. Keep `PLATFORM_OWNER_EMAIL`, `PLATFORM_OWNER_USERNAME`, `PLATFORM_OWNER_DISPLAY_NAME`, `ACCOUNT_TOKEN_ENCRYPTION_KEY` and SMTP credentials in a restricted runtime Secret. Templates must use placeholders, never personal addresses. Generate the encryption key with `Fernet.generate_key()` and back it up securely; rotating it without a migration invalidates stored outbox/MFA ciphertext.

SMTP variables: `EMAIL_HOST`, `EMAIL_PORT`, `EMAIL_HOST_USER`, `EMAIL_HOST_PASSWORD`, `DEFAULT_FROM_EMAIL`, `EMAIL_USE_TLS`, `EMAIL_USE_SSL`. Recipient domain does not select the transport. Gmail requires an authorized app password with two-step verification enabled; the account password must not be used. Never use console email in production.

The encrypted durable outbox can be dispatched after commit through the **existing Celery worker**. The target currently has two independent Redis service endpoints; until broker delivery is proven, use the bounded `account-notifications` maintenance CronJob from `deploy/k8s/accounts-maintenance.yaml` as the authoritative fallback, with `ACCOUNT_EMAIL_TASKS_ENABLED=0`. It runs at most once per minute, at most five messages per run, and also enforces expiry/retention. It is suspended in the template until runtime gates pass; it does not run simulations or create another application. `python manage.py deliver_account_notifications` is a bounded five-message recovery command, not another worker deployment. Failed delivery remains visible and retryable. Retries are bounded, back off and share database budgets: default 100 attempts/hour globally, 10 per identity. `sent` means backend acceptance, never verified receipt. SMTP acknowledgement and the database commit cannot be exactly atomic: a crash after SMTP acceptance can redeliver the same message with a stable Message-ID.

Requests/reset allow 5 attempts/hour per identity/IP (100 global); login 20 (400 global), token confirmation 20/IP (400 global), feedback 10 (200 global). Counters are shared in the database, HMAC keyed and fail closed on backend errors. The immediate proxy supplies the observed network address; forwarded IP headers are not trusted for arbitrary rate-limit identity. Owner audit logs contain codes/IDs, not form data, tokens, query strings or email bodies.

Sensitive queued payloads expire within 24 hours; unverified/rejected request reason/affiliation/role are cleared after 30 days. Tickets are retained privately for follow-up; no clinical dumps or screenshots are collected. Expiry and cleanup run during outbox tasks; operations must schedule the bounded housekeeping command when mail traffic is idle.

## Bootstrap, recovery and restore

After migrations, run `python manage.py bootstrap_platform_owner` once in a controlled Job with the runtime Secret. Never run it on every web startup. It is idempotent and serialized, rejects identity conflicts, preserves existing credentials/MFA/revocation and other administrators, and links an existing identity only when username/email match an already privileged active administrator. A new owner is nonstaff/inactive with an unusable password until the explicit owner activation link is confirmed.

Do not consume Andrea's link or choose his password for a smoke test. Until he completes this himself, record `OWNER_FIRST_LOGIN=PENDING_USER_ACTION`. Standard Django password reset is not the initial activation mechanism.

Break-glass recovery requires authenticated cluster operator access and a verified owner identity: close request/activation gates, run `python manage.py recover_platform_owner --revoke-and-reissue`, then deliver the new owner link. This revokes product sessions/capabilities and makes the old password unusable; MFA is preserved. After verifying the owner identity, lost MFA recovery adds `--reset-authenticator` to that explicit recovery command; it records a separate audit event. Never use a bootstrap rerun or public email claim for this.

After restoring a backup, **before exposing traffic**, run `python manage.py quarantine_restored_accounts --post-restore-quarantine`. It closes gates, increments session versions, disables restored credentials, revokes capabilities and removes restored sensitive outbox payloads. Recover the owner explicitly, then review every restored account. Back up encryption/signing keys with restricted access; never publish restored database contents or tokens.

## Release and rollback

1. Verify SSH target, namespace, current web/worker digest, persistence, and restore capability. The historical local kubeconfig is not authoritative. Do not bypass certificate errors or touch other projects' routing/PVCs.
2. Test migrations against an isolated equivalent database. Create a private online SQLite backup on the existing data volume and validate restore with `PRAGMA quick_check`; keep disaster backup/mirror separate. Save current deployment/config metadata privately.
3. Publish only reviewed source through protected PR checks. Use a digest for both web and Celery. The old automatic mutable-`latest` restart has been replaced by an explicit digest release workflow, which requires these operator preflight steps.
4. For initial onboarding keep gates closed. On routine releases preserve the owner’s current gate choices unless a migration needs a documented pause. Run migrations and idempotent owner bootstrap in a controlled one-shot Job; collect static files to the existing static PVC. Set secure cookies and the verified proxy/HTTPS configuration. Do not reset any existing owner's password.
5. Roll out web and Celery by the same digest, verify readiness, source identity, ordinary routes, local JS and admin CSS, email acceptance/receipt separately, permissions and persistence. Open gates only when the owner can govern access.
6. Rollback: close request/activation gates; set both deployments back to the recorded previous digest and check rollout/readiness. Leave additive account tables/index intact. Do not downgrade schema or restore an older clinical database merely to undo application code. Choose the previous clean product image when available. The pre-registration image does not offer product login; disclose that outage if used.

## Requirements and verification evidence

| Requirement | Implementation | Evidence |
| --- | --- | --- |
| Owner/bootstrap, unique identity, explicit single-use activation | `accounts/services.py`, migrations and management commands | `accounts.tests.test_lifecycle` |
| Ordinary nonstaff login, approval required, no reset bypass | backend, lifecycle and middleware | lifecycle tests; browser journey |
| Existing-site guide and private surfaces denied | existing templates, accounts middleware/DRF permission, simulator selectors | `accounts.tests.test_existing_site` |
| Private suggestions, replies, retry deduplication | SupportTicket and existing-site console | existing-site tests; browser journey |
| Real HTTP forms, CSRF, SMTP sink, desktop/mobile and authenticator | `accounts/tests/browser_journey.py` | explicit isolated browser test command below |
| Production transport/receipt, rollout, owner first login | runtime and human actions | Verified independently on 2026-11-06; see bounded runtime evidence below |

Run disposable backend tests with `BMY_TEST_ROOT=/tmp/bmy-registration-check DJANGO_SETTINGS_MODULE=mmportal.settings_test .venv/bin/python manage.py test accounts.tests --noinput`. The optional browser suite uses `uv sync --frozen --group screenshots`, Chromium, then `BMY_TEST_ROOT=/tmp/bmy-registration-browser DJANGO_SETTINGS_MODULE=mmportal.settings_test .venv/bin/python manage.py test accounts.tests.browser_journey --noinput`. It uses real browser forms and SMTP to a loopback mail sink, not simulated authentication or production mail recipients.

## Verified runtime evidence — 2026-11-06

Protected PRs #80 and #81 are merged. Web (one ready replica), existing Celery (two ready replicas) and account notification maintenance use source `2d464f68341378186c156e413cdcafa2dbfd84cd`, index digest `sha256:a37165e803b9d2a8f70d3782b71c5df71293edaefd0bde86ecdf1733f29dcec6`. Runtime source was read back and `/app/.git` is absent. Source lineage now comes from the explicit build SHA. Routine startup migration/static flags are disabled.

The controlled migration/bootstrap/static Job succeeded with `OWNER_BOOTSTRAP=PRESERVED`. Runtime owner state is ACTIVE with confirmed MFA; the owner completed password setup and login personally and opened both request and activation gates. No external researcher was approved automatically. Existing administrator counts were preserved. Gmail TLS/app-password authentication and restricted Secret persistence passed; the actual owner activation message was accepted by SMTP and independently found in the recipient mailbox. This verifies that message’s receipt, not every later notification or spam-folder placement.

HTTPS login and assets return successful responses, HTTP redirects to HTTPS, CSRF/session secure settings and HSTS are active. Production `check --deploy` retains one existing `security.W019` warning because SAMEORIGIN is used for authenticated embedded plots; this is not silently waived or presented as a zero-warning result. Migration drift check returned no changes. A transient 502 occurred during rollout and subsequent smoke returned 200.

The initial release passed 307 local Django tests and 25 account tests against a disposable SQLite database on the actual NFS target (80.46 seconds total; peak RSS 236024 KiB). Later hardening passed its protected CI and real browser forms/CSRF/password/TOTP with two synthetic users and an isolated SMTP sink. The UX follow-up also runs the wider discovered suite and actual browser calculations/exports. Detailed UX findings and bounded verification are in the [audit](UX_AUDIT_2026_11_06.md).

An online SQLite backup was restored into memory on the target and passed `quick_check`; this is not a complete disaster-recovery drill. Release incident: the previous mutable-latest deployment rescheduled during a pre-existing NotReady-node problem and applied additive migrations through old startup flags before the controlled Job. Gates stayed closed at that time and backups existed. Subsequent deployments use explicit digests and controlled Jobs. No clinical restore, schema downgrade or PVC deletion occurred.

| Delivery control | Status | Scope of evidence |
| --- | --- | --- |
| OWNER_BOOTSTRAP | PASS | Idempotent isolated tests; production rerun preserved credentials/MFA |
| OWNER_FIRST_LOGIN | PASS | Owner report and ACTIVE/confirmed-MFA runtime readback |
| OWNER_ACCESS_MANAGEMENT | PASS | Synthetic owner browser console; owner’s runtime gate actions |
| PUBLIC_ACCESS_REQUEST | PASS | Production form/gates; real HTTP submission in isolated browser |
| EMAIL_VERIFICATION | PASS | Purpose-bound POST verification in sink/browser and negative tests |
| OWNER_APPROVAL_REQUIRED | PASS | No activation before verified review; isolated negative tests |
| SINGLE_USE_ACTIVATION | PASS | Real browser password setup; expiry/reuse/purpose tests |
| NON_STAFF_LOGIN | PASS | Real password browser login and persisted role assertions |
| SANDBOX_ISOLATION | PASS | Own existing learning attempts; no separate sandbox provisioning |
| REAL_DEMO_SIMULATION | PASS | Existing real engine on independent synthetic test scenarios |
| PRIVATE_DATA_DENIAL | PASS | Registered-user clinical/research/API/alias denial; not full #8 closure |
| ARTIFACT_AUTHORIZATION | PASS | Own plot/export; cross-user and suspension checks; direct media closed |
| QUOTAS_AND_CANCELLATION | NOT_TESTED | Existing compute behavior retained; no new budget/cancellation claim |
| USER_FEEDBACK | PASS | Two independent authors; private replies persist across login |
| SESSION_REVOCATION_AND_EXPIRY | PASS | Backend/browser access-state and version tests |
| EMAIL_TRANSPORT | PASS | Actual Gmail TLS authentication and SMTP acceptance |
| EMAIL_RECEIPT | PASS | Owner activation message only, independently verified |
| DESKTOP_AND_MOBILE_E2E | PASS | Isolated browser journey; not human researcher acceptance |
| PRODUCTION_DEPLOYMENT | PASS | Digest/source/replica readback for access hardening above |
| PRODUCTION_SMOKE | PASS | Routes, assets, denial, HTTPS and secure headers; bounded surface |
| ROLLBACK_READINESS | PASS | Prior digest and tested online backup; additive schema retained |
| TRACKER_SYNC | PASS | GitHub/Notion/Trello/Asana dated updates; readback recorded in delivery |
| REAL_RESEARCHER_ACCEPTANCE | PENDING_USER_ACTION | Explicit researcher authorization and human #78 acceptance remain separate |

Release evidence is updated after protected CI and runtime verification. The abandoned parallel draft and its tests are not acceptance evidence for this implementation. No researcher is approved automatically; #78 requires Andrea's explicit authorization and real human acceptance.
