#!/bin/bash
# Start MM Portal development environment
# This script starts Redis, Celery worker, and Django server

set -e

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_DIR"

# Activate virtual environment
if [ -d "venv" ]; then
    source venv/bin/activate
    echo -e "${GREEN}✓ Virtual environment activated${NC}"
else
    echo -e "${RED}✗ Virtual environment not found!${NC}"
    exit 1
fi

# Check if Redis is installed
if ! command -v redis-server &> /dev/null; then
    echo -e "${YELLOW}⚠ Redis not found. Install with: brew install redis${NC}"
    echo -e "${YELLOW}  Running without Celery (synchronous mode)${NC}"
    export CELERY_TASK_ALWAYS_EAGER=True
    python manage.py runserver 8001
    exit 0
fi

# Function to cleanup background processes on exit
cleanup() {
    echo -e "\n${YELLOW}Shutting down services...${NC}"
    if [ ! -z "$REDIS_PID" ]; then
        kill $REDIS_PID 2>/dev/null || true
    fi
    if [ ! -z "$CELERY_PID" ]; then
        kill $CELERY_PID 2>/dev/null || true
    fi
    echo -e "${GREEN}✓ Services stopped${NC}"
    exit 0
}

trap cleanup SIGINT SIGTERM

# Start Redis in background
echo -e "${YELLOW}Starting Redis...${NC}"
redis-server --daemonize yes --port 6379 --dir /tmp --loglevel warning
sleep 1
REDIS_PID=$(pgrep -f "redis-server" | head -1)
echo -e "${GREEN}✓ Redis started (PID: $REDIS_PID)${NC}"

# Start Celery worker in background
echo -e "${YELLOW}Starting Celery worker...${NC}"
celery -A mmportal worker --loglevel=info --logfile=logs/celery.log --detach
sleep 2
CELERY_PID=$(pgrep -f "celery.*worker" | head -1)
echo -e "${GREEN}✓ Celery worker started (PID: $CELERY_PID)${NC}"

# Start Django development server
echo -e "${YELLOW}Starting Django server on http://127.0.0.1:8001${NC}"
echo -e "${GREEN}✓ All services running!${NC}"
echo -e "${YELLOW}Press Ctrl+C to stop all services${NC}"
echo ""

python manage.py runserver 8001

# Cleanup will be called automatically on Ctrl+C
