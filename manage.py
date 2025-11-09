#!/usr/bin/env python3
"""Django's command-line utility for administrative tasks."""
import os
import sys
import subprocess
import signal
import time
import atexit


# Global process references for cleanup
_redis_process = None
_celery_process = None


def _start_background_services():
    """Start Redis and Celery automatically when running development server."""
    global _redis_process, _celery_process
    
    # Only start services for runserver command
    if len(sys.argv) < 2 or sys.argv[1] != 'runserver':
        return
    
    print("\n🚀 Starting background services...\n")
    
    # Check if Redis is already running
    try:
        result = subprocess.run(
            ['redis-cli', 'ping'],
            capture_output=True,
            timeout=1,
        )
        if result.returncode == 0:
            print("✅ Redis already running")
        else:
            raise subprocess.CalledProcessError(1, 'redis-cli')
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError):
        # Try to start Redis
        try:
            print("📦 Starting Redis server...")
            _redis_process = subprocess.Popen(
                ['redis-server', '--port', '6379', '--dir', '/tmp', '--loglevel', 'warning'],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            time.sleep(1)
            
            # Verify Redis started
            result = subprocess.run(['redis-cli', 'ping'], capture_output=True, timeout=1)
            if result.returncode == 0:
                print(f"✅ Redis started (PID: {_redis_process.pid})")
            else:
                print("⚠️  Redis failed to start - using synchronous mode")
                os.environ['CELERY_TASK_ALWAYS_EAGER'] = 'True'
                return
        except FileNotFoundError:
            print("⚠️  Redis not installed - using synchronous mode")
            print("   Install with: brew install redis")
            os.environ['CELERY_TASK_ALWAYS_EAGER'] = 'True'
            return
    
    # Start Celery worker
    try:
        print("📦 Starting Celery worker...")
        from pathlib import Path
        logs_dir = Path('logs')
        logs_dir.mkdir(exist_ok=True)
        
        # Use sys.executable to ensure same Python interpreter
        python_exe = sys.executable
        _celery_process = subprocess.Popen(
            [python_exe, '-m', 'celery', '-A', 'mmportal', 'worker', '--loglevel=info'],
            stdout=open(logs_dir / 'celery.log', 'w'),
            stderr=subprocess.STDOUT,
        )
        time.sleep(2)
        
        if _celery_process.poll() is None:
            print(f"✅ Celery worker started (PID: {_celery_process.pid})")
        else:
            print("⚠️  Celery failed to start")
    except FileNotFoundError:
        print("⚠️  Celery not found - install with: pip install celery")
    except Exception as e:
        print(f"⚠️  Celery error: {e}")
    
    print("\n✨ All services ready!\n")


def _cleanup_background_services():
    """Stop background services on exit."""
    global _redis_process, _celery_process
    
    if _celery_process and _celery_process.poll() is None:
        print("\n📦 Stopping Celery worker...")
        _celery_process.terminate()
        try:
            _celery_process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            _celery_process.kill()
    
    if _redis_process and _redis_process.poll() is None:
        print("📦 Stopping Redis server...")
        _redis_process.terminate()
        try:
            _redis_process.wait(timeout=5)
        except subprocess.TimeoutExpired:
            _redis_process.kill()


def _signal_handler(signum, frame):
    """Handle interrupt signals."""
    _cleanup_background_services()
    sys.exit(0)


def main() -> None:
    """Run administrative tasks."""
    os.environ.setdefault("DJANGO_SETTINGS_MODULE", "mmportal.settings")
    
    try:
        from django.core.management import execute_from_command_line
    except ImportError as exc:  # pragma: no cover - django import guard
        raise ImportError(
            "Couldn't import Django. Are you sure it's installed and "
            "available on your PYTHONPATH environment variable? Did you "
            "forget to activate a virtual environment?"
        ) from exc
    
    # Start background services if running development server
    _start_background_services()
    
    # Register cleanup handlers
    atexit.register(_cleanup_background_services)
    signal.signal(signal.SIGINT, _signal_handler)
    signal.signal(signal.SIGTERM, _signal_handler)
    
    execute_from_command_line(sys.argv)


if __name__ == "__main__":
    main()