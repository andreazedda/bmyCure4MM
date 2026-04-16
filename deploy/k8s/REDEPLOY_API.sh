#!/bin/bash
# Direct k3s API redeploy - uses Bearer token authentication

set -e

TOKEN="sha256~R77HyEr-Xc6WwY2WpI4QHZZEYdWvNPRZeFBeLujHaSg"
API="https://10.56.245.110:6443"
NAMESPACE="bmycure4mm"

echo "======================================"
echo "Direct k3s API Redeploy"
echo "======================================"
echo ""

# Force pod restart by updating deployment annotation
echo "[1/2] Restarting web deployment via Kubernetes API..."
curl -sSk \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/strategic-merge-patch+json" \
  -d "{\"spec\":{\"template\":{\"metadata\":{\"annotations\":{\"restartedAt\":\"$(date -u +'%Y-%m-%dT%H:%M:%SZ')\"}}}}}" \
  -X PATCH \
  "$API/apis/apps/v1/namespaces/$NAMESPACE/deployments/web"

echo ""
echo "[2/2] Restarting celery deployment via Kubernetes API..."
curl -sSk \
  -H "Authorization: Bearer $TOKEN" \
  -H "Content-Type: application/strategic-merge-patch+json" \
  -d "{\"spec\":{\"template\":{\"metadata\":{\"annotations\":{\"restartedAt\":\"$(date -u +'%Y-%m-%dT%H:%M:%SZ')\"}}}}}" \
  -X PATCH \
  "$API/apis/apps/v1/namespaces/$NAMESPACE/deployments/celery"

echo ""
echo "[✓] Redeploy triggered!"
echo ""
echo "Verify at: https://bmycure4mm.clusterlab.uk/admin/login/"
