"""
Django management command to start all MM Portal services.

Usage:
    python manage.py runservices
"""
import os
import signal
import subprocess
import sys
import time
from pathlib import Path

from django.core.management.base import BaseCommand
from django.core.management import call_command


class Command(BaseCommand):
    help = 'Start all MM Portal services (Redis, Celery, Django)'

    def add_arguments(self, parser):
        parser.add_argument(
            '--port',
            type=int,
            default=8001,
            help='Port for Django development server (default: 8001)',
        )
        parser.add_argument(
            '--no-redis',
            action='store_true',
            help='Run without Redis/Celery (synchronous mode)',
        )
        parser.add_argument(
            '--no-celery',
            action='store_true',
            help='Run without Celery worker',
        )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.redis_process = None
        self.celery_process = None

    def handle(self, *args, **options):
        port = options['port']
        no_redis = options['no_redis']
        no_celery = options['no_celery']

        # Register signal handlers for cleanup
        signal.signal(signal.SIGINT, self._signal_handler)
        signal.signal(signal.SIGTERM, self._signal_handler)

        try:
            if no_redis:
                self.stdout.write(self.style.WARNING(
                    '⚠️  Running in synchronous mode (no Redis/Celery)'
                ))
                os.environ['CELERY_TASK_ALWAYS_EAGER'] = 'True'
                os.environ['CELERY_TASK_EAGER_PROPAGATES'] = 'True'
            else:
                # Start Redis
                if not self._start_redis():
                    self.stdout.write(self.style.ERROR(
                        '✗ Failed to start Redis. Install with: brew install redis'
                    ))
                    self.stdout.write(self.style.WARNING(
                        '  Falling back to synchronous mode...'
                    ))
                    os.environ['CELERY_TASK_ALWAYS_EAGER'] = 'True'
                    no_celery = True

                # Start Celery
                if not no_celery and not self._start_celery():
                    self.stdout.write(self.style.ERROR(
                        '✗ Failed to start Celery worker'
                    ))
                    return

            # Start Django server
            self.stdout.write(self.style.SUCCESS(
                f'✓ Starting Django server on http://127.0.0.1:{port}'
            ))
            self.stdout.write(self.style.WARNING('Press Ctrl+C to stop all services'))
            self.stdout.write('')

            call_command('runserver', f'127.0.0.1:{port}')

        except KeyboardInterrupt:
            self.stdout.write('\n')
            self._cleanup()
        except Exception as e:
            self.stdout.write(self.style.ERROR(f'Error: {e}'))
            self._cleanup()
            sys.exit(1)

    def _start_redis(self):
        """Start Redis server in background."""
        # Check if Redis is already running
        try:
            subprocess.run(
                ['redis-cli', 'ping'],
                capture_output=True,
                timeout=1,
                check=True,
            )
            self.stdout.write(self.style.SUCCESS('✓ Redis already running'))
            return True
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired, FileNotFoundError):
            pass

        # Try to start Redis
        try:
            self.stdout.write('Starting Redis server...')
            self.redis_process = subprocess.Popen(
                ['redis-server', '--port', '6379', '--dir', '/tmp', '--loglevel', 'warning'],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            time.sleep(1)

            # Verify Redis started
            result = subprocess.run(
                ['redis-cli', 'ping'],
                capture_output=True,
                timeout=1,
            )
            if result.returncode == 0:
                self.stdout.write(self.style.SUCCESS(
                    f'✓ Redis started (PID: {self.redis_process.pid})'
                ))
                return True
            else:
                self.stdout.write(self.style.ERROR('✗ Redis failed to start'))
                return False

        except FileNotFoundError:
            return False
        except Exception as e:
            self.stdout.write(self.style.ERROR(f'✗ Redis error: {e}'))
            return False

    def _start_celery(self):
        """Start Celery worker in background."""
        try:
            self.stdout.write('Starting Celery worker...')

            # Create logs directory
            logs_dir = Path('logs')
            logs_dir.mkdir(exist_ok=True)

            self.celery_process = subprocess.Popen(
                ['celery', '-A', 'mmportal', 'worker', '--loglevel=info'],
                stdout=open(logs_dir / 'celery.log', 'w'),
                stderr=subprocess.STDOUT,
            )
            time.sleep(2)

            if self.celery_process.poll() is None:
                self.stdout.write(self.style.SUCCESS(
                    f'✓ Celery worker started (PID: {self.celery_process.pid})'
                ))
                return True
            else:
                self.stdout.write(self.style.ERROR('✗ Celery failed to start'))
                return False

        except FileNotFoundError:
            self.stdout.write(self.style.ERROR('✗ Celery not found. Install with: pip install celery'))
            return False
        except Exception as e:
            self.stdout.write(self.style.ERROR(f'✗ Celery error: {e}'))
            return False

    def _cleanup(self):
        """Stop all background services."""
        self.stdout.write(self.style.WARNING('Shutting down services...'))

        if self.celery_process and self.celery_process.poll() is None:
            self.stdout.write('Stopping Celery worker...')
            self.celery_process.terminate()
            try:
                self.celery_process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                self.celery_process.kill()

        if self.redis_process and self.redis_process.poll() is None:
            self.stdout.write('Stopping Redis server...')
            self.redis_process.terminate()
            try:
                self.redis_process.wait(timeout=5)
            except subprocess.TimeoutExpired:
                self.redis_process.kill()

        self.stdout.write(self.style.SUCCESS('✓ Services stopped'))

    def _signal_handler(self, signum, frame):
        """Handle interrupt signals."""
        raise KeyboardInterrupt
