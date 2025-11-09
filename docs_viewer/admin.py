"""
Admin interface for documentation viewer.
"""
from django.contrib import admin
from .models import DocumentView, DocumentFeedback


@admin.register(DocumentView)
class DocumentViewAdmin(admin.ModelAdmin):
    """Admin for document views tracking."""
    
    list_display = ['path', 'user', 'viewed_at']
    list_filter = ['viewed_at', 'path']
    search_fields = ['path', 'user__username']
    date_hierarchy = 'viewed_at'
    readonly_fields = ['path', 'user', 'viewed_at']
    
    def has_add_permission(self, request):
        """Disable manual addition."""
        return False


@admin.register(DocumentFeedback)
class DocumentFeedbackAdmin(admin.ModelAdmin):
    """Admin for document feedback."""
    
    list_display = ['path', 'rating', 'user', 'created_at']
    list_filter = ['rating', 'created_at']
    search_fields = ['path', 'comment', 'user__username']
    date_hierarchy = 'created_at'
    readonly_fields = ['path', 'user', 'created_at']

