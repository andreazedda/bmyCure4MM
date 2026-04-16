#!/bin/bash
# Quick redeploy script for k3s cluster
# Forces pulling the latest image and restarting deployments

set -eu

NAMESPACE="bmycure4mm"

echo "[*] Restarting web and celery deployments to pull latest image..."
sudo kubectl -n "$NAMESPACE" rollout restart deployment/web deployment/celery

echo "[*] Waiting for web deployment to be ready..."
sudo kubectl -n "$NAMESPACE" rollout status deployment/web --timeout=5m

echo "[*] Waiting for celery deployment to be ready..."
sudo kubectl -n "$NAMESPACE" rollout status deployment/celery --timeout=5m

echo "[*] Checking pod status..."
sudo kubectl -n "$NAMESPACE" get pods -o wide

echo "[*] Checking web logs..."
sudo kubectl -n "$NAMESPACE" logs -f deployment/web --tail=50 &
LOGS_PID=$!

sleep 3
kill $LOGS_PID 2>/dev/null || true

echo "[✓] Redeploy complete! Latest image is now running."
echo "[*] Verify at: https://bmycure4mm.clusterlab.uk/admin/login/"
