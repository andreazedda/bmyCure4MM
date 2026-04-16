#!/bin/bash
# Quick setup script with demo data for learning the platform

set -e  # Exit on error

ADMIN_USERNAME="admin"
ADMIN_EMAIL="admin@example.com"
ADMIN_PASSWORD="$(python3 - <<'PY'
import secrets
import string

alphabet = string.ascii_letters + string.digits
print(''.join(secrets.choice(alphabet) for _ in range(20)))
PY
)"

echo "🚀 bmyCure4MM Quick Start Setup"
echo "================================"
echo ""

# Check if virtual environment is active
if [ -z "$VIRTUAL_ENV" ]; then
    echo "⚠️  Virtual environment not active. Activating..."
    if [ -f "venv/bin/activate" ]; then
        source venv/bin/activate
    else
        echo "❌ Error: venv not found. Run: python3 -m venv venv && source venv/bin/activate"
        exit 1
    fi
fi

echo "✓ Virtual environment active"
echo ""

# Run migrations
echo "📦 Applying database migrations..."
python manage.py migrate --noinput
echo "✓ Migrations applied"
echo ""

# Load demo data
echo "👥 Loading demo patients with assessments..."
python manage.py loaddata clinic/fixtures/demo_patients.json
echo "✓ Demo data loaded"
echo ""

# Create superuser if not exists (non-interactive)
echo "👤 Checking for superuser..."
export ADMIN_USERNAME ADMIN_EMAIL ADMIN_PASSWORD
python manage.py shell << EOF
import os
from django.contrib.auth import get_user_model
User = get_user_model()
username = os.environ['ADMIN_USERNAME']
email = os.environ['ADMIN_EMAIL']
password = os.environ['ADMIN_PASSWORD']
if not User.objects.filter(username=username).exists():
    User.objects.create_superuser(username, email, password)
    print(f'✓ Superuser created: username={username}, password={password}')
else:
    print('✓ Superuser already exists')
EOF
unset ADMIN_USERNAME ADMIN_EMAIL ADMIN_PASSWORD
echo ""

echo "🎉 Setup Complete!"
echo "=================="
echo ""
echo "🌐 Start the server with:"
echo "   python manage.py runserver"
echo ""
echo "Then visit: http://127.0.0.1:8000"
echo ""
echo "📚 Tutorial: Check the dashboard for the '🚀 New to Platform?' card"
echo ""
echo "👥 Demo patients available:"
echo "   - Mario Rossi (R-ISS III, high risk)"
echo "   - Anna Bianchi (R-ISS II, intermediate risk)"
echo "   - Giuseppe Verdi (R-ISS I, low risk)"
echo ""
echo "🔐 Admin login:"
echo "   Username: admin"
echo "   Password: the generated password printed above during setup"
echo "   URL: http://127.0.0.1:8000/admin/"
echo ""
echo "✨ Ready to simulate treatments and see Patient Twin in action!"
