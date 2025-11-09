"""
Models for documentation viewer.
"""
from django.db import models


class DocumentView(models.Model):
    """Track documentation views for analytics."""
    
    path = models.CharField(max_length=255, db_index=True)
    viewed_at = models.DateTimeField(auto_now_add=True, db_index=True)
    user = models.ForeignKey(
        'auth.User',
        on_delete=models.SET_NULL,
        null=True,
        blank=True
    )
    
    class Meta:
        verbose_name = "Document View"
        verbose_name_plural = "Document Views"
        ordering = ['-viewed_at']
    
    def __str__(self):
        return f"{self.path} at {self.viewed_at}"


class DocumentFeedback(models.Model):
    """User feedback on documentation quality."""
    
    RATING_CHOICES = [
        (1, '⭐'),
        (2, '⭐⭐'),
        (3, '⭐⭐⭐'),
        (4, '⭐⭐⭐⭐'),
        (5, '⭐⭐⭐⭐⭐'),
    ]
    
    path = models.CharField(max_length=255, db_index=True)
    rating = models.IntegerField(choices=RATING_CHOICES)
    comment = models.TextField(blank=True)
    user = models.ForeignKey(
        'auth.User',
        on_delete=models.SET_NULL,
        null=True,
        blank=True
    )
    created_at = models.DateTimeField(auto_now_add=True)
    
    class Meta:
        verbose_name = "Document Feedback"
        verbose_name_plural = "Document Feedback"
        ordering = ['-created_at']
        unique_together = [['path', 'user']]
    
    def __str__(self):
        return f"{self.path} - {self.get_rating_display()}"
