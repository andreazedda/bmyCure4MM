# Existing-cluster release: bmyCure4MM

The canonical authentication, email and recovery contract is
[AUTHENTICATION_AND_ACCESS.md](../../docs/operations/AUTHENTICATION_AND_ACCESS.md).
This page contains the concrete operator sequence for the verified cluster.
Daily approvals and tickets use the product console, without these commands.

## Verified target and release

- Origin: `https://bmycure4mm.clusterlab.uk`.
- Namespace: `bmycure4mm`; existing deployments `web` and `celery`.
- Project Caddy proxy in `public-edge`, separate static service/PVC, SQLite
  `/data/db.sqlite3`, existing media/data/static NFS PVCs.
- SSH aliases already configured: `crypto-ms-7c95`, `lenovo3`; do not use a stale
  local kubecontext or bypass certificate verification.
- Access/QR/email hardening source: `2d464f68341378186c156e413cdcafa2dbfd84cd`,
  digest `sha256:a37165e803b9d2a8f70d3782b71c5df71293edaefd0bde86ecdf1733f29dcec6`.
- Reviewed UX/docs/private-result release: PR #82, tested head
  `8d02b441e71b9e2c908dcd79829dc71554a6437d`, merged source
  `55a556b0c5cf4b8dc95ef6be6d0a6bc873cc1a3c`, published digest
  `sha256:56a4beb7e54153f7cba30af4de3e5240feb2f23265e9702b33f0f7e19c43f80d`.


- Previous completed runtime release: PR #84, tested head
  `2a19c73cf149b7afcee2f853e43e38cee5c3d255`, merged source
  `abbb98f8bb5c0f121d9469a19bee36a5e3e2a5a8`, deployed index digest
  `sha256:47bb2811561bd34a385df91e3bdb942576ddb07940951e36b3ea317929465d66`.
  Web and two workers are ready on `lenovo2`; cgroup ceilings match the budgets
  below. Controlled Job `accounts-release-abbb98f` and fresh notification
  maintenance passed; owner/MFA/gates and persistence were preserved.
- Intermediate deployed release: [PR #86](https://github.com/andreazedda/bmyCure4MM/pull/86),
  tested head `49d1a89b12cd0508e8c6aba2340cb193a31eed11`, merged source
  `313e2e4d8a839b9c14a8d853d1f008071a6cb013`, published index digest
  `sha256:0d707ef6c20ee3f6cc4650730610ee045af56955130ab28a5bb98b7bdccfc405`.
  Required/post-merge checks, Pages and GHCR passed; OCI labels and image
  environment match that source. Local verification: 473 Django tests, two
  browser tests in 38.646s, unchanged frozen scientific baseline.
- Target verification `ux-verify-313e2e4`: **PASS**, 83 tests on disposable SQLite
  on the actual NFS PVC (34 accounts, 46 docs, three learning-entry tests),
  38.559s test cases / 96.69s runner, peak RSS 269340 KiB at 1 CPU / 512 MiB.
  Target SMTP/TLS connectivity passed; this does not establish authentication,
  message acceptance or receipt. Controlled Job `accounts-release-313e2e4`
  **passed**: source matched, backup restored with `quick_check=ok`, ACTIVE
  owner/MFA and both open owner-selected gates were preserved, static assets
  were present and the previously empty catalogue has one active synthetic
  scenario. Web and both workers report the expected source/digest at unchanged
  1 CPU / 1 GiB ceilings; owner/MFA/gates/catalogue persisted. All 20 public
  route/asset probes passed.
- **Historical PR #86 static-cache gate FAIL:** returning browsers received stale unversioned
  CSS despite correct collected files; the cache lifetime was 14400 seconds.
  PR #87 introduces source-versioned asset URLs. Accept the corrected release
  only after its protected source checks, published digest, rollout and
  returning-browser cache check pass. The **dated release-evidence comment in
  [GitHub #79](https://github.com/andreazedda/bmyCure4MM/issues/79)** records the
  current deployed source/digest, asset hashes and cache/notification outcomes
  after readback, with the corresponding PR #87 evidence. This source snapshot
  preserves an intermediate finding and the operator procedure; it does not
  predict a future image hash or certify an unobserved rollout.

On 2026-11-06 the control plane became unstable: `lenovo-master` was NotReady
and unreachable over SSH, other API endpoints repeatedly failed etcd readiness
with no leader, and a worker briefly failed GHCR DNS lookup. Image publication
and tests do not establish deployment. The dated tracker update is authoritative
for the latest rollout/smoke outcome. Restore cluster readiness before mutation;
restarting other projects, resetting etcd or deleting PVCs is outside this release.

## 1. Preflight

```bash
ssh -o BatchMode=yes -o StrictHostKeyChecking=yes lenovo3 \
  'sudo -n k3s kubectl --request-timeout=15s get --raw=/readyz'
ssh -o BatchMode=yes -o StrictHostKeyChecking=yes lenovo3 \
  'sudo -n k3s kubectl --request-timeout=15s -n bmycure4mm get deployment web celery -o wide'
```

Require a stable `ok` readiness result. Verify namespace, current digests,
node health/capacity, NFS capacity and recoverable backup metadata. `lenovo2`
was verified Ready without pressure/taints and already serves the project's
static PVC. Recheck it; do not override taints or increase resource budgets.
Use the image digest, never mutable `latest`. The package was verified public;
no broad personal GitHub token needs to be copied into the cluster.

Preserve `bmycure4mm-secrets` and `bmycure4mm-accounts`. **Do not delete or
recreate secrets to release application code.** Replacing the Django signing
key or account encryption key can invalidate sessions, MFA and encrypted mail.
Do not reset an activated owner's password or close gates on routine reruns.

## 2. Render and run the controlled release Job

The Job backs up the actual SQLite database online, tests restoring that backup
into memory, then runs additive migrations, idempotent owner bootstrap,
`ensure_learning_example` and collectstatic on the existing PVC. The catalogue
command uses the existing policy lock and creates exactly one independent
synthetic `Scenario` only when the entire catalogue is empty. An existing
catalogue, including one with only inactive scenarios, is preserved. It does
not create a patient, change scientific defaults or allocate per-user copies.
The Job does not send passwords, run simulations or introduce a new application.
A failed backup check prevents migrations.

Set `BMY_RELEASE_ID` to a unique, reviewed source/attempt identifier and
`BMY_RELEASE_DIGEST` to the corresponding verified digest from the dated release
ledger before running this sequence. The PR #86 digest above records an
intermediate deployment; do not reuse it to claim that the later cache fix has
been deployed. The script rejects unset/placeholder values before rendering.

```bash
set -eu
: "${BMY_RELEASE_ID:?Set the reviewed source/attempt identifier}"
: "${BMY_RELEASE_DIGEST:?Set the verified sha256 image digest}"
export BMY_RELEASE_ID BMY_RELEASE_DIGEST
export BMY_RELEASE_NODE=lenovo2
python3 - <<'PY' > /tmp/bmy-accounts-release.yaml
import os, re
from pathlib import Path
assert re.fullmatch(r'[a-z0-9][a-z0-9-]{0,40}', os.environ['BMY_RELEASE_ID'])
assert re.fullmatch(r'sha256:[a-f0-9]{64}', os.environ['BMY_RELEASE_DIGEST'])
source = Path('deploy/k8s/accounts-release-job.yaml').read_text()
for key in ('ID', 'NODE', 'DIGEST'):
    source = source.replace('RELEASE_' + key, os.environ['BMY_RELEASE_' + key])
assert 'RELEASE_' not in source
print(source)
PY
ssh -o BatchMode=yes -o StrictHostKeyChecking=yes lenovo3 \
  'sudo -n k3s kubectl --request-timeout=20s apply -f -' < /tmp/bmy-accounts-release.yaml
```

Inspect that named Job until completion; do not treat creation as success.
Read only its sanitized completion metadata. Require `backup_restore=ok`, the
expected source, preserved owner/MFA/gates, static assets and the sanitized
`active_scenarios` count. If an intentionally inactive catalogue is preserved,
do not silently activate it to make a smoke check pass. Backups stay in
the private data volume. They are not test fixtures or public artifacts.
A failed Job needs diagnosis and a new attempt ID; never recreate the database.

## 3. Roll out the existing application

Use the existing protected release workflow with the verified digest **after**
the Job passes, or its equivalent scoped commands:

```bash
ssh -o BatchMode=yes -o StrictHostKeyChecking=yes lenovo3 \
  "sudo -n k3s kubectl --request-timeout=20s -n bmycure4mm set image deployment/web web=ghcr.io/andreazedda/bmycure4mm@${BMY_RELEASE_DIGEST}"
ssh -o BatchMode=yes -o StrictHostKeyChecking=yes lenovo3 \
  "sudo -n k3s kubectl --request-timeout=20s -n bmycure4mm set image deployment/celery celery=ghcr.io/andreazedda/bmycure4mm@${BMY_RELEASE_DIGEST}"
```

Keep web/Celery startup `DJANGO_MIGRATE=0` and `DJANGO_COLLECTSTATIC=0`, existing
Secret references, resource limits and storage. If a node fails, schedule only
these project workloads on a verified eligible node; do not force-remove
unreachable workers or terminate jobs belonging to other users.

Render `accounts-maintenance.yaml` with the same digest, retain the verified
eligible node selection, apply it, then explicitly resume the CronJob after
readiness. Its template is suspended intentionally. Confirm a **new successful**
notification run; an old success timestamp is not current delivery evidence.
Keep `ACCOUNT_EMAIL_TASKS_ENABLED=0` until the existing broker topology is proven.

## 4. Bound the existing workloads

Do not apply the historical aggregate `bmycure4mm.yaml`: it does not describe the
current NFS/Caddy topology. The two small strategic patches below touch only
resource requests/limits of the existing containers; they preserve image,
replicas, environment, Secret references, routing and PVCs.

```bash
ssh -o BatchMode=yes -o StrictHostKeyChecking=yes lenovo3 \
  'sudo -n k3s kubectl --request-timeout=20s -n bmycure4mm patch deployment web --type strategic --patch-file /dev/stdin' < deploy/k8s/web-budget.patch.yaml
ssh -o BatchMode=yes -o StrictHostKeyChecking=yes lenovo3 \
  'sudo -n k3s kubectl --request-timeout=20s -n bmycure4mm patch deployment celery --type strategic --patch-file /dev/stdin' < deploy/k8s/celery-budget.patch.yaml
```

| Existing container | CPU request / limit | Memory request / limit |
| --- | --- | --- |
| Web, one replica | 100m / 1 CPU | 384 MiB / 1 GiB |
| Celery, each of two replicas | 100m / 1 CPU | 256 MiB / 1 GiB |
| Notification maintenance | 100m / 500m | 128 MiB / 512 MiB |

Before applying, the actual web/worker resource dictionaries were empty. Idle
measurements on `lenovo2` were 251 MiB web and 131 MiB per worker. The bounded
real-engine/plot/CSV/ownership/suspension integration path on disposable SQLite
on the actual NFS volume passed: 6.401 seconds in the test body, 62.04 seconds
including fresh test-database setup, peak RSS 265408 KiB (about 259 MiB). Its Job
kept the existing 1 CPU / 512 MiB ceiling. The first attempt inherited production
HTTPS redirection into a local HTTP test client and failed before calculation;
only the disposable test environment was corrected, without raising limits or
weakening production HTTPS.

These are conservative container ceilings for the verified account/learning
path, not a claim that all advanced workloads or global queue/admission and
cancellation policies are validated. Recheck rollout, memory and current owner
state; do not raise ceilings just to pass a test. Existing advanced workload
capacity and broker health remain separate governance work.

## 5. Acceptance and rollback

Require current web/worker source and digest readback, HTTPS and redirect,
health/readiness, account JS, admin CSS/JS, docs CSS/figures, owner state and
email outbox status. Check source-versioned asset URLs in a returning browser
with an existing cache, as well as fresh desktop/mobile sessions: correct
collectstatic files and a fresh HTTP 200 do not prove that cached clients have
the new stylesheet. Never consume the
owner's activation token or take a screenshot of a real authenticator seed.

For PR #86 the rollback target is the previous production source
`abbb98f8bb5c0f121d9469a19bee36a5e3e2a5a8`, digest
`sha256:47bb2811561bd34a385df91e3bdb942576ddb07940951e36b3ea317929465d66`.
This restores the earlier application behavior, including the defects fixed
by PR #86; it is not equivalent security acceptance.

1. In `/manage/access-requests/`, close new requests and activations, then
   suspend every externally approved member before reverting. Closing gates
   does not terminate existing members' access. Verify protected pages/exports
   are denied and keep owner governance available. No member is reopened
   automatically after rollback.
2. Repin web, Celery and `account-notifications` to the old digest, preserving
   current Secret references, replicas, eligible node, budgets and PVCs. Check
   rollout/readiness and source identity for web and both workers.
3. Restore the old image's matching static files with a **separate bounded
   one-shot Job** whose sole command is
   `python manage.py collectstatic --noinput`. Use the old digest, existing
   ConfigMap/Secret references and existing static PVC, with the verified node
   and at most the release Job's 1 CPU / 512 MiB and 300-second deadline.
   **Do not use the new `accounts-release-job.yaml` command with the old image:**
   that image has no `ensure_learning_example` command. Routine web startup
   remains `DJANGO_MIGRATE=0` and `DJANGO_COLLECTSTATIC=0`.
4. Check public/login routes, account/admin/docs assets and MIME types,
   owner/MFA access, suspended-member denial and a fresh maintenance outcome.
   Keep the seeded synthetic scenario; do not delete records to imitate an
   earlier database state.

Leave additive account tables/index in place. No schema downgrade, clinical
database restore, secret replacement or PVC deletion is part of an application
rollback. Actual disaster restore requires the separate account/session/token
quarantine procedure before exposing traffic. A loss of cluster/etcd readiness
must be recovered at the cluster level before release or rollback can be verified.

Historical PR #84 rollback target was source `55a556b`, digest
`sha256:56a4beb7e54153f7cba30af4de3e5240feb2f23265e9702b33f0f7e19c43f80d`;
it remains here only for release traceability.
