"""
Tests for Django settings security requirements.

Validates that the application properly enforces SECRET_KEY security.
"""

import os
import sys
from pathlib import Path
from unittest import TestCase
from unittest.mock import patch

import django
from django.core.exceptions import ImproperlyConfigured


class SecretKeySecurityTests(TestCase):
    """Test SECRET_KEY security validation."""

    def setUp(self):
        """Set up test environment."""
        # Store original environment
        self.original_env = os.environ.copy()
        
    def tearDown(self):
        """Restore original environment."""
        os.environ.clear()
        os.environ.update(self.original_env)
        
        # Clear Django settings module cache
        if 'mmportal.settings' in sys.modules:
            del sys.modules['mmportal.settings']

    def test_missing_secret_key_raises_error(self):
        """Test that missing DJANGO_SECRET_KEY raises ImproperlyConfigured."""
        # Remove SECRET_KEY from environment
        os.environ.pop('DJANGO_SECRET_KEY', None)
        
        with self.assertRaises(ImproperlyConfigured) as context:
            import mmportal.settings
            
        self.assertIn('DJANGO_SECRET_KEY', str(context.exception))
        self.assertIn('required', str(context.exception).lower())

    def test_insecure_secret_key_raises_error(self):
        """Test that SECRET_KEY containing 'insecure' raises ImproperlyConfigured."""
        os.environ['DJANGO_SECRET_KEY'] = 'django-insecure-test-key'
        
        with self.assertRaises(ImproperlyConfigured) as context:
            import mmportal.settings
            
        self.assertIn('insecure', str(context.exception).lower())

    def test_default_secret_key_raises_error(self):
        """Test that default SECRET_KEY values raise ImproperlyConfigured."""
        default_keys = [
            'django-insecure-please-change-me',
            'your-secret-key-here-change-me-in-production',
            'change-me',
            'changeme',
        ]
        
        for key in default_keys:
            with self.subTest(key=key):
                os.environ['DJANGO_SECRET_KEY'] = key
                
                # Clear module cache
                if 'mmportal.settings' in sys.modules:
                    del sys.modules['mmportal.settings']
                
                with self.assertRaises(ImproperlyConfigured) as context:
                    import mmportal.settings
                    
                self.assertIn('insecure', str(context.exception).lower())

    def test_valid_secret_key_succeeds(self):
        """Test that a valid SECRET_KEY is accepted."""
        # Generate a valid-looking key
        os.environ['DJANGO_SECRET_KEY'] = 'django-valid-xh8$p9k#2m@5n!qwerty-random-key-12345'
        
        try:
            import mmportal.settings
            # If we get here, settings loaded successfully
            self.assertTrue(True)
        except ImproperlyConfigured:
            self.fail("Valid SECRET_KEY should not raise ImproperlyConfigured")

    def test_error_message_includes_generation_command(self):
        """Test that error message includes SECRET_KEY generation instructions."""
        os.environ.pop('DJANGO_SECRET_KEY', None)
        
        with self.assertRaises(ImproperlyConfigured) as context:
            import mmportal.settings
            
        error_msg = str(context.exception)
        self.assertIn('get_random_secret_key', error_msg)
        self.assertIn('python -c', error_msg)


if __name__ == '__main__':
    import unittest
    unittest.main()
