# MM Portal Development Setup

## Quick Start (Recommended)

The easiest way to start the portal is:

```bash
./start_dev.sh
```

This automatically starts all services needed for the portal.

## Prerequisites

### Install Redis (Required for async jobs)

**macOS:**
```bash
brew install redis
```

**Ubuntu/Debian:**
```bash
sudo apt-get update
sudo apt-get install redis-server
```

**Windows:**
Download from https://redis.io/download or use WSL2

### Verify Installation

```bash
redis-server --version
```

## Starting the Portal

### Method 1: All-in-one Script (Easiest)

```bash
./start_dev.sh
```

This starts:
- ✅ Redis server (port 6379)
- ✅ Celery worker (for background jobs)
- ✅ Django server (port 8001)

Press Ctrl+C to stop all services.

### Method 2: Without Redis (Development Mode)

If you don't want to install Redis:

```bash
export CELERY_TASK_ALWAYS_EAGER=True
uv run python manage.py runserver 8001
```

**Note:** Jobs will run synchronously (blocking the request until complete).

### Method 3: Legacy Supervisord scripts

The legacy `manage_services.sh` path depends on Supervisor, which is not part
of the canonical lock and is not a supported reproducible environment. Use the
Docker Compose path below for a production-like local stack.

### Method 4: Docker Compose

```bash
docker-compose up -d          # Start in background
docker-compose logs -f web    # View logs
docker-compose down           # Stop all
```

### Method 5: Tailscale-local Django access

For local tailnet access from another device, use the dedicated helper and setup guide:

```bash
./scripts/run_tailscale_dev.sh
PORT=8000 ./scripts/run_tailscale_dev.sh
```

Guide: [TAILSCALE_LOCALHOST_ACCESS.md](TAILSCALE_LOCALHOST_ACCESS.md)

## Research Cockpit Development

- Research cockpit guide: [../research/RESEARCH_COCKPIT_USER_GUIDE.md](../research/RESEARCH_COCKPIT_USER_GUIDE.md)
- Developer console: [DEVELOPER_CONSOLE.md](DEVELOPER_CONSOLE.md)
- Pre-push safety gate: [PRE_PUSH_SAFETY.md](PRE_PUSH_SAFETY.md)
- Test strategy: [TEST_STRATEGY_RESEARCH_COCKPIT.md](TEST_STRATEGY_RESEARCH_COCKPIT.md)

## Accessing the Portal

Open your browser to: **http://127.0.0.1:8001**

Default credentials:
- Username: `admin`
- Password: (check with your team)

## Troubleshooting

### "Job stuck in Queued status"

**Cause:** Celery worker is not running

**Solution:**
```bash
# Check if Celery is running
ps aux | grep celery

# If not running, start it
celery -A mmportal worker --loglevel=info
```

### "Connection refused to Redis"

**Cause:** Redis server is not running

**Solution:**
```bash
# Check if Redis is running
ps aux | grep redis

# Start Redis manually
redis-server
```

### "Module not found" errors

**Cause:** Virtual environment not activated or dependencies not installed

**Solution:**
```bash
uv sync --frozen --extra chemistry
```

## Development Workflow

1. **Start services once:**
   ```bash
   ./start_dev.sh
   ```

2. **Make code changes** - Django auto-reloads

3. **Run tests:**
   ```bash
   python manage.py test
   ```

   Notes:
   - The repo’s test suite is designed to run via Django’s built-in test runner.
   - Optional/binary dependencies (e.g., RDKit) are best-effort: if they are unavailable or incompatible in your environment, tests should still run.

4. **Stop services:** Press Ctrl+C

## Service Management

### Check Service Status

```bash
# Redis
ps aux | grep redis-server

# Celery
ps aux | grep "celery.*worker"

# Django
ps aux | grep "manage.py runserver"
```

### Stop Individual Services

```bash
# Stop Redis
redis-cli shutdown

# Stop Celery
pkill -f "celery.*worker"
```

## Production Deployment

For production, use:
- Gunicorn/uWSGI instead of `runserver`
- Supervisor or systemd for process management
- Nginx as reverse proxy
- PostgreSQL instead of SQLite

See `docs/DEPLOYMENT.md` for details.
