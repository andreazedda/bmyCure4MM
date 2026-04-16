#!/bin/bash
# Direct k3s redeploy script - run on k3s server or from machine with kubectl access

set -eu

NAMESPACE="bmycure4mm"
REGISTRY="ghcr.io/andreazedda/bmycure4mm"
IMAGE_TAG="latest"

echo "========================================"
echo "BMW Cure 4MM - k3s Redeploy Script"
echo "========================================"
echo ""

# Force pull latest image
echo "[1/4] Pulling latest image from GHCR..."
docker pull "$REGISTRY:$IMAGE_TAG" 2>&1 | tail -5 || echo "Note: Docker pull failed (might be running via containerd/k3s)"

# Restart deployments
echo "[2/4] Restarting web deployment..."
kubectl -n "$NAMESPACE" rollout restart deployment/web

echo "[3/4] Restarting celery deployment..."
kubectl -n "$NAMESPACE" rollout restart deployment/celery

# Wait for rollout
echo "[4/4] Waiting for deployments to be ready (max 5 minutes)..."
if kubectl -n "$NAMESPACE" rollout status deployment/web --timeout=5m; then
    echo "[✓] Web deployment ready!"
else
    echo "[⚠] Web deployment rollout timed out, but deployment was triggered"
fi

if kubectl -n "$NAMESPACE" rollout status deployment/celery --timeout=5m; then
    echo "[✓] Celery deployment ready!"
else
    echo "[⚠] Celery deployment rollout timed out, but deployment was triggered"
fi

echo ""
echo "========================================"
echo "[✓] Redeploy complete!"
echo "========================================"
echo ""
echo "Pod status:"
kubectl -n "$NAMESPACE" get pods -o wide
echo ""
echo "Recent events:"
kubectl -n "$NAMESPACE" get events --sort-by='.lastTimestamp' | tail -10
echo ""
echo "✓ Visit: https://bmycure4mm.clusterlab.uk/admin/login/"
