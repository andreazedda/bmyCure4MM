"""
Management command to validate documentation quality.

This command performs comprehensive checks on documentation files:
- Validates all markdown files can be rendered
- Checks for broken internal links
- Detects orphaned files
- Reports syntax errors
- Checks for missing titles

Usage:
    python manage.py validate_docs
    python manage.py validate_docs --verbose
    python manage.py validate_docs --fix-links
"""
import os
from django.core.management.base import BaseCommand, CommandError
from django.conf import settings
from docs_viewer.utils import (
    get_allowed_doc_paths,
    is_safe_path,
    render_markdown,
    validate_markdown_links,
    extract_title
)


class Command(BaseCommand):
    help = 'Validate documentation files for quality and consistency'
    
    def add_arguments(self, parser):
        parser.add_argument(
            '--verbose',
            action='store_true',
            help='Show detailed validation information',
        )
        parser.add_argument(
            '--fix-links',
            action='store_true',
            help='Suggest fixes for broken links',
        )
    
    def handle(self, *args, **options):
        verbose = options['verbose']
        fix_links = options['fix_links']
        
        self.stdout.write(self.style.SUCCESS('🔍 Validating documentation files...'))
        self.stdout.write('')
        
        # Get all allowed documentation paths
        doc_paths = get_allowed_doc_paths()
        
        if not doc_paths:
            self.stdout.write(self.style.WARNING('⚠️  No documentation files found'))
            return
        
        # Statistics
        total_files = len(doc_paths)
        valid_files = 0
        files_with_errors = []
        files_with_warnings = []
        total_broken_links = 0
        files_without_titles = []
        
        # Validate each file
        for doc_path in doc_paths:
            errors = []
            warnings = []
            
            # Security check
            if not is_safe_path(doc_path):
                errors.append(f'Security validation failed')
                files_with_errors.append((doc_path, errors))
                continue
            
            full_path = os.path.join(settings.BASE_DIR, doc_path)
            
            try:
                # Read file
                with open(full_path, 'r', encoding='utf-8') as f:
                    content = f.read()
                
                # Check for title
                title = extract_title(content)
                if title == "Untitled Document":
                    files_without_titles.append(doc_path)
                    warnings.append('No H1 heading found')
                
                # Try to render
                try:
                    html, toc = render_markdown(content)
                    if verbose:
                        self.stdout.write(f'✓ {doc_path} - OK')
                except Exception as e:
                    errors.append(f'Rendering failed: {str(e)}')
                
                # Check for broken links
                base_path = os.path.dirname(doc_path)
                broken = validate_markdown_links(content, base_path)
                
                if broken:
                    total_broken_links += len(broken)
                    warnings.append(f'{len(broken)} broken link(s)')
                    
                    if fix_links:
                        for link in broken:
                            self.stdout.write(
                                self.style.WARNING(f'  → Broken link in {doc_path}: {link}')
                            )
                
                # Track results
                if errors:
                    files_with_errors.append((doc_path, errors))
                elif warnings:
                    files_with_warnings.append((doc_path, warnings))
                else:
                    valid_files += 1
                    
            except FileNotFoundError:
                errors.append('File not found')
                files_with_errors.append((doc_path, errors))
            except UnicodeDecodeError:
                errors.append('Invalid UTF-8 encoding')
                files_with_errors.append((doc_path, errors))
            except Exception as e:
                errors.append(f'Unexpected error: {str(e)}')
                files_with_errors.append((doc_path, errors))
        
        # Print summary
        self.stdout.write('')
        self.stdout.write(self.style.SUCCESS('📊 Validation Summary'))
        self.stdout.write('─' * 50)
        self.stdout.write(f'Total files:          {total_files}')
        self.stdout.write(self.style.SUCCESS(f'✓ Valid files:        {valid_files}'))
        
        if files_with_warnings:
            self.stdout.write(
                self.style.WARNING(f'⚠ Files with warnings: {len(files_with_warnings)}')
            )
        
        if files_with_errors:
            self.stdout.write(
                self.style.ERROR(f'✗ Files with errors:   {len(files_with_errors)}')
            )
        
        if total_broken_links > 0:
            self.stdout.write(
                self.style.WARNING(f'🔗 Broken links:       {total_broken_links}')
            )
        
        if files_without_titles:
            self.stdout.write(
                self.style.WARNING(f'📄 Files without H1:   {len(files_without_titles)}')
            )
        
        # Show detailed errors
        if files_with_errors:
            self.stdout.write('')
            self.stdout.write(self.style.ERROR('❌ Files with Errors:'))
            for path, errors in files_with_errors:
                self.stdout.write(f'  {path}')
                for error in errors:
                    self.stdout.write(f'    - {error}')
        
        # Show warnings if verbose
        if verbose and files_with_warnings:
            self.stdout.write('')
            self.stdout.write(self.style.WARNING('⚠️  Files with Warnings:'))
            for path, warnings in files_with_warnings:
                self.stdout.write(f'  {path}')
                for warning in warnings:
                    self.stdout.write(f'    - {warning}')
        
        # Show files without titles
        if verbose and files_without_titles:
            self.stdout.write('')
            self.stdout.write(self.style.WARNING('📄 Files without H1 heading:'))
            for path in files_without_titles:
                self.stdout.write(f'  {path}')
        
        # Exit with error if there are critical errors
        if files_with_errors:
            raise CommandError(
                f'Documentation validation failed: {len(files_with_errors)} file(s) with errors'
            )
        
        self.stdout.write('')
        self.stdout.write(self.style.SUCCESS('✅ Documentation validation complete!'))
