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
export BMY_RELEASE_ID=55a556b
export BMY_RELEASE_NODE=lenovo2
export BMY_RELEASE_DIGEST=sha256:56a4beb7e54153f7cba30af4de3e5240feb2f23265e9702b33f0f7e19c43f80d
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
  'sudo -n k3s kubectl --request-timeout=20s -n bmycure4mm set image deployment/web web=ghcr.io/andreazedda/bmycure4mm@sha256:56a4beb7e54153f7cba30af4de3e5240feb2f23265e9702b33f0f7e19c43f80d'
ssh -o BatchMode=yes -o StrictHostKeyChecking=yes lenovo3 \
  'sudo -n k3s kubectl --request-timeout=20s -n bmycure4mm set image deployment/celery celery=ghcr.io/andreazedda/bmycure4mm@sha256:56a4beb7e54153f7cba30af4de3e5240feb2f23265e9702b33f0f7e19c43f80d'
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

## 4. Acceptance and rollback

Require current web/worker source and digest readback, HTTPS and redirect,
health/readiness, account JS, admin CSS/JS, docs CSS/figures, owner state and
email outbox status. Use desktop/mobile browser checks. Never consume the
owner's activation token or take a screenshot of a real authenticator seed.

For this release the previous clean product image is the `2d464f6` digest above.
Pause new requests/activations through the owner console if rollback needs an
outage; repin web, Celery and notification maintenance to that digest and verify
readiness. Leave additive account tables/index in place. No schema downgrade,
clinical database restore, secret replacement or PVC deletion is part of an
application rollback. A loss of cluster/etcd readiness must be recovered at the
cluster level before either release or rollback can be verified.
