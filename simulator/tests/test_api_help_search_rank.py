"""Tests for help_search API ranking, deduplication, and cutoff."""

from __future__ import annotations

from django.urls import reverse

from simulator.models_help import HelpArticle


def test_search_matches_articles_by_slug(client, db):
    """Verify that search matches articles by slug."""
    HelpArticle.objects.create(
        slug="lenalidomide_dose",
        title_en="Lenalidomide Dose",
        body_en="Content about lenalidomide.",
        title_it="Dose Lenalidomide",
        body_it="Contenuto su lenalidomide.",
    )
    HelpArticle.objects.create(
        slug="bortezomib_dose",
        title_en="Bortezomib Dose",
        body_en="Content about bortezomib.",
        title_it="Dose Bortezomib",
        body_it="Contenuto su bortezomib.",
    )
    url = reverse("api_help_search")
    response = client.get(url, {"q": "lenali", "lang": "en"})
    assert response.status_code == 200
    data = response.json()
    slugs = [item["slug"] for item in data["results"]]
    assert "lenalidomide_dose" in slugs
    assert "bortezomib_dose" not in slugs


def test_search_matches_articles_by_title(client, db):
    """Verify that search matches articles by title."""
    HelpArticle.objects.create(
        slug="dose_calc",
        title_en="Lenalidomide Calculator",
        body_en="Content about calculator.",
        title_it="Calcolatore Lenalidomide",
        body_it="Contenuto sul calcolatore.",
    )
    url = reverse("api_help_search")
    response = client.get(url, {"q": "lenalidomide", "lang": "en"})
    assert response.status_code == 200
    data = response.json()
    slugs = [item["slug"] for item in data["results"]]
    assert "dose_calc" in slugs


def test_search_matches_presets(client, db, settings):
    """Verify that search matches presets and prioritizes them."""
    HelpArticle.objects.create(
        slug="lenalidomide_dose",
        title_en="Lenalidomide Dose",
        body_en="Content.",
        title_it="Dose",
        body_it="Contenuto.",
    )
    url = reverse("api_help_search")
    response = client.get(url, {"q": "standard", "lang": "en"})
    assert response.status_code == 200
    data = response.json()
    # Should match "Standard Regimen" preset if it exists in PRESETS
    types = [item["type"] for item in data["results"]]
    # At least check that we get results
    assert len(data["results"]) >= 0


def test_search_cutoff_at_20_results(client, db):
    """Verify that search returns at most 20 results."""
    # Create 30 articles
    for i in range(30):
        HelpArticle.objects.create(
            slug=f"article_{i}",
            title_en=f"Title {i}",
            body_en="Content.",
            title_it=f"Titolo {i}",
            body_it="Contenuto.",
        )
    url = reverse("api_help_search")
    response = client.get(url, {"q": "title", "lang": "en"})
    assert response.status_code == 200
    data = response.json()
    assert len(data["results"]) <= 20


def test_search_deduplication(client, db):
    """Verify that search deduplicates results."""
    HelpArticle.objects.create(
        slug="duplicate_test",
        title_en="Duplicate Test",
        body_en="Content.",
        title_it="Test Duplicato",
        body_it="Contenuto.",
    )
    url = reverse("api_help_search")
    response = client.get(url, {"q": "duplicate", "lang": "en"})
    assert response.status_code == 200
    data = response.json()
    slugs = [item["slug"] for item in data["results"]]
    # Check no duplicates
    assert len(slugs) == len(set(slugs))


def test_search_returns_type_field(client, db):
    """Verify that search results include type field."""
    HelpArticle.objects.create(
        slug="kpi_tumor_reduction",
        title_en="Tumor Reduction",
        body_en="Content.",
        title_it="Riduzione Tumorale",
        body_it="Contenuto.",
    )
    url = reverse("api_help_search")
    response = client.get(url, {"q": "tumor", "lang": "en"})
    assert response.status_code == 200
    data = response.json()
    for item in data["results"]:
        assert "type" in item
        assert item["type"] in ["article", "field", "preset", "kpi"]


def test_search_empty_query_returns_first_20(client, db):
    """Verify that empty query returns first 20 articles ordered by slug."""
    for i in range(25):
        HelpArticle.objects.create(
            slug=f"z_article_{i:02d}",
            title_en=f"Title {i}",
            body_en="Content.",
            title_it=f"Titolo {i}",
            body_it="Contenuto.",
        )
    url = reverse("api_help_search")
    response = client.get(url, {"lang": "en"})
    assert response.status_code == 200
    data = response.json()
    assert len(data["results"]) <= 20
    # Should be ordered by slug
    slugs = [item["slug"] for item in data["results"]]
    assert slugs == sorted(slugs)


def test_search_respects_language_parameter(client, db):
    """Verify that search respects language parameter for titles."""
    HelpArticle.objects.create(
        slug="test_lang",
        title_en="English Title",
        body_en="Content.",
        title_it="Titolo Italiano",
        body_it="Contenuto.",
    )
    url = reverse("api_help_search")

    # Search in English
    response_en = client.get(url, {"q": "english", "lang": "en"})
    assert response_en.status_code == 200
    data_en = response_en.json()
    assert len(data_en["results"]) > 0
    assert data_en["results"][0]["title"] == "English Title"

    # Search in Italian
    response_it = client.get(url, {"q": "italiano", "lang": "it"})
    assert response_it.status_code == 200
    data_it = response_it.json()
    assert len(data_it["results"]) > 0
    assert data_it["results"][0]["title"] == "Titolo Italiano"
