"""
Django management command to debug and test Celery jobs.
Usage: python manage.py debug_jobs
"""
import logging
from django.core.management.base import BaseCommand
from chemtools import models, tasks

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(name)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)


class Command(BaseCommand):
    help = 'Debug Celery jobs and test task execution'

    def add_arguments(self, parser):
        parser.add_argument(
            '--job-id',
            type=int,
            help='Test specific job by ID',
        )
        parser.add_argument(
            '--test-similarity',
            action='store_true',
            help='Test similarity search with sample SMILES',
        )
        parser.add_argument(
            '--list-queued',
            action='store_true',
            help='List all queued jobs',
        )
        parser.add_argument(
            '--retry-all',
            action='store_true',
            help='Retry all queued jobs',
        )
        parser.add_argument(
            '--mark-failed',
            action='store_true',
            help='Mark stuck jobs as failed',
        )

    def handle(self, *args, **options):
        if options['list_queued']:
            self._list_queued_jobs()
        
        elif options['job_id']:
            self._test_job(options['job_id'])
        
        elif options['test_similarity']:
            self._test_similarity_search()
        
        elif options['retry_all']:
            self._retry_all_queued()
        
        elif options['mark_failed']:
            self._mark_stuck_as_failed()
        
        else:
            self.stdout.write(self.style.WARNING('No action specified. Use --help for options.'))

    def _list_queued_jobs(self):
        """List all jobs in queued state."""
        self.stdout.write(self.style.WARNING('\n=== QUEUED JOBS ===\n'))
        
        queued_jobs = models.ChemJob.objects.filter(
            out_html='',
            out_csv='',
        ).exclude(log__startswith='ERROR')
        
        if not queued_jobs:
            self.stdout.write(self.style.SUCCESS('No queued jobs found.'))
            return
        
        for job in queued_jobs:
            self.stdout.write(f'\nJob ID: {job.pk}')
            self.stdout.write(f'  Type: {job.get_kind_display()}')
            self.stdout.write(f'  Created: {job.created}')
            self.stdout.write(f'  Input: {job.input_a}')
            self.stdout.write(f'  Progress: {job.progress_percent}% - {job.progress_message}')
            self.stdout.write(f'  Log preview: {job.log[:200] if job.log else "(empty)"}...')

    def _test_job(self, job_id):
        """Test a specific job by running it directly."""
        self.stdout.write(self.style.WARNING(f'\n=== TESTING JOB {job_id} ===\n'))
        
        try:
            job = models.ChemJob.objects.get(pk=job_id)
            self.stdout.write(f'Job found: {job.get_kind_display()}')
            self.stdout.write(f'Input A: {job.input_a}')
            self.stdout.write(f'Input B: {job.input_b}')
            
            self.stdout.write(self.style.WARNING('\nRunning task directly (synchronous)...'))
            
            if job.kind == models.ChemJob.SIM:
                tasks.run_similarity_job(job.pk, job.input_a)
            elif job.kind == models.ChemJob.PARAM:
                tasks.run_drug_params_job(job.pk, job.input_a, job.input_b)
            elif job.kind == models.ChemJob.BIND:
                tasks.run_binding_viz_job(job.pk, job.input_a, job.input_b)
            
            job.refresh_from_db()
            self.stdout.write(self.style.SUCCESS('\n✓ Task completed!'))
            self.stdout.write(f'Progress: {job.progress_percent}%')
            self.stdout.write(f'Message: {job.progress_message}')
            self.stdout.write(f'Has output: HTML={bool(job.out_html)}, CSV={bool(job.out_csv)}')
            
            if job.log:
                self.stdout.write(self.style.WARNING('\nLast 500 chars of log:'))
                self.stdout.write(job.log[-500:])
                
        except models.ChemJob.DoesNotExist:
            self.stdout.write(self.style.ERROR(f'Job {job_id} not found!'))
        except Exception as e:
            self.stdout.write(self.style.ERROR(f'\n✗ Error: {e}'))
            import traceback
            self.stdout.write(traceback.format_exc())

    def _test_similarity_search(self):
        """Create and test a new similarity search job."""
        self.stdout.write(self.style.WARNING('\n=== TESTING SIMILARITY SEARCH ===\n'))
        
        # Sample SMILES for Ibuprofen
        test_smiles = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
        
        self.stdout.write(f'Creating test job with SMILES: {test_smiles}')
        
        job = models.ChemJob.objects.create(
            kind=models.ChemJob.SIM,
            input_a=test_smiles,
        )
        
        self.stdout.write(self.style.SUCCESS(f'Created job {job.pk}'))
        
        # Run it directly
        self._test_job(job.pk)

    def _retry_all_queued(self):
        """Retry all queued jobs."""
        self.stdout.write(self.style.WARNING('\n=== RETRYING ALL QUEUED JOBS ===\n'))
        
        queued_jobs = models.ChemJob.objects.filter(
            out_html='',
            out_csv='',
        ).exclude(log__startswith='ERROR')
        
        if not queued_jobs:
            self.stdout.write(self.style.SUCCESS('No queued jobs to retry.'))
            return
        
        for job in queued_jobs:
            self.stdout.write(f'\nRetrying job {job.pk} ({job.get_kind_display()})...')
            try:
                if job.kind == models.ChemJob.SIM:
                    tasks.run_similarity_job.delay(job.pk, job.input_a)
                elif job.kind == models.ChemJob.PARAM:
                    tasks.run_drug_params_job.delay(job.pk, job.input_a, job.input_b)
                elif job.kind == models.ChemJob.BIND:
                    tasks.run_binding_viz_job.delay(job.pk, job.input_a, job.input_b)
                self.stdout.write(self.style.SUCCESS(f'  ✓ Job {job.pk} queued'))
            except Exception as e:
                self.stdout.write(self.style.ERROR(f'  ✗ Failed: {e}'))

    def _mark_stuck_as_failed(self):
        """Mark jobs stuck in processing as failed."""
        self.stdout.write(self.style.WARNING('\n=== MARKING STUCK JOBS AS FAILED ===\n'))
        
        from django.utils import timezone
        from datetime import timedelta
        
        # Jobs older than 10 minutes without output
        cutoff = timezone.now() - timedelta(minutes=10)
        
        stuck_jobs = models.ChemJob.objects.filter(
            created__lt=cutoff,
            out_html='',
            out_csv='',
        ).exclude(log__startswith='ERROR')
        
        if not stuck_jobs:
            self.stdout.write(self.style.SUCCESS('No stuck jobs found.'))
            return
        
        for job in stuck_jobs:
            self.stdout.write(f'Marking job {job.pk} as failed...')
            job.log = f"ERROR: Job stuck in processing. Marked as failed.\n\n{job.log or ''}"
            job.progress_percent = 0
            job.progress_message = "Failed - timeout"
            job.save()
            self.stdout.write(self.style.SUCCESS(f'  ✓ Job {job.pk} marked as failed'))
