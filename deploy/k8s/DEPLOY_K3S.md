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


- Current verified application: PR #84, tested head
  `2a19c73cf149b7afcee2f853e43e38cee5c3d255`, merged source
  `abbb98f8bb5c0f121d9469a19bee36a5e3e2a5a8`, deployed index digest
  `sha256:47bb2811561bd34a385df91e3bdb942576ddb07940951e36b3ea317929465d66`.
  Web and two workers are ready on `lenovo2`; cgroup ceilings match the budgets
  below. Controlled Job `accounts-release-abbb98f` and fresh notification
  maintenance passed; owner/MFA/gates and persistence were preserved.

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
into memory, then runs additive migrations, idempotent owner bootstrap and
collectstatic on the existing PVC. It does not send passwords, run simulations
or introduce a new application. A failed backup check prevents migrations.

```bash
export BMY_RELEASE_ID=abbb98f
export BMY_RELEASE_NODE=lenovo2
export BMY_RELEASE_DIGEST=sha256:47bb2811561bd34a385df91e3bdb942576ddb07940951e36b3ea317929465d66
python3 - <<'PY' > /tmp/bmy-accounts-release.yaml
import os
from pathlib import Path
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
expected source, preserved owner/MFA/gates and static assets. Backups stay in
the private data volume. They are not test fixtures or public artifacts.
A failed Job needs diagnosis and a new attempt ID; never recreate the database.

## 3. Roll out the existing application

Use the existing protected release workflow with the verified digest **after**
the Job passes, or its equivalent scoped commands:

```bash
ssh -o BatchMode=yes -o StrictHostKeyChecking=yes lenovo3 \
  'sudo -n k3s kubectl --request-timeout=20s -n bmycure4mm set image deployment/web web=ghcr.io/andreazedda/bmycure4mm@sha256:47bb2811561bd34a385df91e3bdb942576ddb07940951e36b3ea317929465d66'
ssh -o BatchMode=yes -o StrictHostKeyChecking=yes lenovo3 \
  'sudo -n k3s kubectl --request-timeout=20s -n bmycure4mm set image deployment/celery celery=ghcr.io/andreazedda/bmycure4mm@sha256:47bb2811561bd34a385df91e3bdb942576ddb07940951e36b3ea317929465d66'
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
email outbox status. Use desktop/mobile browser checks. Never consume the
owner's activation token or take a screenshot of a real authenticator seed.

For PR #84 the previous clean product image is source `55a556b`, digest
`sha256:56a4beb7e54153f7cba30af4de3e5240feb2f23265e9702b33f0f7e19c43f80d`.
Pause new requests/activations through the owner console if rollback needs an
outage; repin web, Celery and notification maintenance to that digest and verify
readiness. Leave additive account tables/index in place. No schema downgrade,
clinical database restore, secret replacement or PVC deletion is part of an
application rollback. A loss of cluster/etcd readiness must be recovered at the
cluster level before either release or rollback can be verified.
