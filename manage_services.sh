#!/bin/bash
# Manage MM Portal services with supervisord
# Usage: ./manage_services.sh [start|stop|restart|status]

set -e

PROJECT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$PROJECT_DIR"

# Colors
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
RED='\033[0;31m'
NC='\033[0m'

# Check if supervisord is installed
if ! command -v supervisord &> /dev/null; then
    echo -e "${RED}✗ supervisord not found!${NC}"
    echo -e "${YELLOW}Install with: pip install supervisor${NC}"
    exit 1
fi

# Create logs directory
mkdir -p logs

case "${1:-start}" in
    start)
        echo -e "${YELLOW}Starting all services...${NC}"
        supervisord -c supervisord.conf
        sleep 2
        supervisorctl -c supervisord.conf status
        echo -e "${GREEN}✓ All services started!${NC}"
        echo -e "${YELLOW}Access the portal at: http://127.0.0.1:8001${NC}"
        echo -e "${YELLOW}View logs in: logs/${NC}"
        ;;
    stop)
        echo -e "${YELLOW}Stopping all services...${NC}"
        supervisorctl -c supervisord.conf shutdown
        echo -e "${GREEN}✓ All services stopped${NC}"
        ;;
    restart)
        echo -e "${YELLOW}Restarting all services...${NC}"
        supervisorctl -c supervisord.conf restart all
        echo -e "${GREEN}✓ All services restarted${NC}"
        ;;
    status)
        supervisorctl -c supervisord.conf status
        ;;
    logs)
        echo -e "${YELLOW}Showing recent logs (Ctrl+C to stop)...${NC}"
        tail -f logs/*.log
        ;;
    *)
        echo "Usage: $0 {start|stop|restart|status|logs}"
        exit 1
        ;;
esac
