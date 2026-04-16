# INSTANT REDEPLOY - Copy & Paste Solution

The Django admin page will load correctly after this redeploy. All code is fixed and image is built.

## For k3s Server Admin (Run on k3s server):

```bash
# SSH to the k3s server, then run:
sudo kubectl -n bmycure4mm rollout restart deployment/web deployment/celery
sudo kubectl -n bmycure4mm rollout status deployment/web --timeout=5m
sudo kubectl -n bmycure4mm rollout status deployment/celery --timeout=5m

# Verify pods are running:
sudo kubectl -n bmycure4mm get pods -o wide

# Check logs:
sudo kubectl -n bmycure4mm logs -f deployment/web --tail=50
```

## Result Expected:
- Pods will restart (old pods terminate, new pods start)
- New pods will download `ghcr.io/andreazedda/bmycure4mm:latest`
- Entrypoint.sh will run `collectstatic` (populates /app/staticfiles)
- Nginx will serve CSS/JavaScript correctly
- Admin page at `https://bmycure4mm.clusterlab.uk/admin/login/` will load with full styling

## If You Don't Have k3s Server Access:

Contact your infrastructure team to run the kubectl commands above.

Alternatively, apply the manifest directly:
```bash
sudo kubectl apply -f deploy/k8s/bmycure4mm.yaml --force
```

## Verification
After redeploy, visit: https://bmycure4mm.clusterlab.uk/admin/login/

You should see:
- ✅ Proper CSS styling (not black rectangles)
- ✅ Theme toggle buttons working
- ✅ Django admin login form with proper layout
