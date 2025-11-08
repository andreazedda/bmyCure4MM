"""Tests for help_item API caching with ETag and 304 responses."""

from __future__ import annotations

import re

from django.urls import reverse

from simulator.models_help import HelpArticle


def test_help_item_returns_etag_header(client, db):
    """Verify that help_item returns an ETag header."""
    article = HelpArticle.objects.create(
        slug="test_article",
        title_en="Test Title",
        body_en="Test body content.",
        title_it="Titolo Test",
        body_it="Contenuto del corpo test.",
    )
    url = reverse("api_help_item", kwargs={"slug": article.slug}) + "?lang=en"
    response = client.get(url)
    assert response.status_code == 200
    assert "ETag" in response.headers
    etag = response.headers["ETag"]
    # ETag should be a 40-character hex string (salted_hmac output)
    assert re.fullmatch(r"[0-9a-f]{40}", etag)


def test_help_item_returns_304_with_matching_etag(client, db):
    """Verify that help_item returns 304 when If-None-Match matches ETag."""
    article = HelpArticle.objects.create(
        slug="test_article",
        title_en="Test Title",
        body_en="Test body content.",
        title_it="Titolo Test",
        body_it="Contenuto del corpo test.",
    )
    url = reverse("api_help_item", kwargs={"slug": article.slug}) + "?lang=en"

    # First request to get ETag
    response1 = client.get(url)
    assert response1.status_code == 200
    etag = response1.headers["ETag"]

    # Second request with If-None-Match header
    response2 = client.get(url, HTTP_IF_NONE_MATCH=etag)
    assert response2.status_code == 304
    assert response2.content == b""


def test_help_item_returns_cache_control_header(client, db):
    """Verify that help_item returns Cache-Control header."""
    article = HelpArticle.objects.create(
        slug="test_article",
        title_en="Test Title",
        body_en="Test body content.",
        title_it="Titolo Test",
        body_it="Contenuto del corpo test.",
    )
    url = reverse("api_help_item", kwargs={"slug": article.slug}) + "?lang=en"
    response = client.get(url)
    assert response.status_code == 200
    assert "Cache-Control" in response.headers
    assert "max-age=600" in response.headers["Cache-Control"]
    assert "public" in response.headers["Cache-Control"]


def test_help_item_returns_last_modified_header(client, db):
    """Verify that help_item returns Last-Modified header."""
    article = HelpArticle.objects.create(
        slug="test_article",
        title_en="Test Title",
        body_en="Test body content.",
        title_it="Titolo Test",
        body_it="Contenuto del corpo test.",
    )
    url = reverse("api_help_item", kwargs={"slug": article.slug}) + "?lang=en"
    response = client.get(url)
    assert response.status_code == 200
    assert "Last-Modified" in response.headers


def test_help_item_etag_changes_after_update(client, db):
    """Verify that ETag changes when article is updated."""
    article = HelpArticle.objects.create(
        slug="test_article",
        title_en="Test Title",
        body_en="Test body content.",
        title_it="Titolo Test",
        body_it="Contenuto del corpo test.",
    )
    url = reverse("api_help_item", kwargs={"slug": article.slug}) + "?lang=en"

    # First request
    response1 = client.get(url)
    etag1 = response1.headers["ETag"]

    # Update article
    article.body_en = "Updated body content."
    article.save()

    # Second request after update
    response2 = client.get(url)
    etag2 = response2.headers["ETag"]

    # ETags should be different
    assert etag1 != etag2


def test_help_item_etag_differs_by_language(client, db):
    """Verify that ETag differs for different languages."""
    article = HelpArticle.objects.create(
        slug="test_article",
        title_en="Test Title",
        body_en="Test body content.",
        title_it="Titolo Test",
        body_it="Contenuto del corpo test.",
    )

    url_en = reverse("api_help_item", kwargs={"slug": article.slug}) + "?lang=en"
    url_it = reverse("api_help_item", kwargs={"slug": article.slug}) + "?lang=it"

    response_en = client.get(url_en)
    response_it = client.get(url_it)

    etag_en = response_en.headers["ETag"]
    etag_it = response_it.headers["ETag"]

    # ETags should differ for different languages
    assert etag_en != etag_it
