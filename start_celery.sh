#!/bin/bash
# Start Celery worker for MM Portal
# Usage: ./start_celery.sh

# Activate virtual environment
source venv/bin/activate

# Start Celery worker
echo "Starting Celery worker..."
echo "Make sure Redis is running: redis-server"
echo ""

celery -A mmportal worker --loglevel=info

# Press Ctrl+C to stop
