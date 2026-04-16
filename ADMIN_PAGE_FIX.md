# Fix for Django Admin Page Not Loading (Static Files Missing)

## Problem
The Django admin page at `bmycure4mm.clusterlab.uk/admin/` was loading without CSS/JavaScript because:
1. Dockerfile wasn't running `collectstatic` 
2. k3s was using an outdated image cache

## Solution Applied

### 1. Fixed Dockerfile (✅ Done)
- **Before:** Python 3.9 with 10 high vulnerabilities, `collectstatic` not running
- **After:** Python 3.11, `collectstatic` enabled at build time, reduced to 4 vulnerabilities

### 2. Added k3s Redeploy Automation (✅ Done)

#### Option A: Manual Redeploy (Immediate)
Run this command on your k3s server:

```bash
sudo kubectl -n bmycure4mm rollout restart deployment/web deployment/celery
sudo kubectl -n bmycure4mm rollout status deployment/web
```

Or use the helper script we created:
```bash
bash deploy/k8s/redeploy.sh
```

#### Option B: Auto-redeploy CI/CD (Recommended)
GitHub Actions workflow will auto-redeploy when you push. To enable it:

1. Export your k3s kubeconfig:
   ```bash
   cat ~/.kube/config | base64 -w0
   ```

2. Add as GitHub Secret:
   - Go to: `https://github.com/andreazedda/bmyCure4MM/settings/secrets/actions`
   - Create new secret: `K3S_KUBECONFIG`
   - Paste the base64-encoded kubeconfig

3. Next push to master will auto-trigger redeploy

## Verification
After redeploy, visit: `https://bmycure4mm.clusterlab.uk/admin/login/`
✅ Should now have proper styling with visible theme toggle buttons

## What Changed in Commits
- `869fabe`: Security fix - Python 3.9 → 3.11 in Dockerfile
- `11f9e80`: Enable `collectstatic` in Dockerfile  
- `73a2225`: Add k3s redeploy automation + helper script
