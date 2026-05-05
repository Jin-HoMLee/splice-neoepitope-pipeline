"""Tests for zotero_add.crossref_to_zotero — preprint vs journal-article paths."""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from zotero_add import crossref_to_zotero, _first_or_empty, _is_preprint


_COLLECTION = "Z38GTJNW"
_TAGS = ["morning-reading"]


def _journal_crossref():
    return {
        "type": "journal-article",
        "title": ["Splice neoantigens in cancer"],
        "author": [{"given": "Jin-Ho", "family": "Lee"}],
        "container-title": ["Bioinformatics Advances"],
        "volume": "4",
        "issue": "1",
        "page": "vbae080",
        "published-print": {"date-parts": [[2024, 6, 15]]},
        "DOI": "10.1093/bioadv/vbae080",
        "URL": "https://doi.org/10.1093/bioadv/vbae080",
        "ISSN": ["2635-0041"],
        "abstract": "<jats:p>Cancer abstract.</jats:p>",
    }


def _biorxiv_crossref():
    return {
        "type": "posted-content",
        "subtype": "preprint",
        "title": ["Benchmarking foundation models for splice site annotation"],
        "author": [{"given": "Some", "family": "Author"}],
        "container-title": [],
        "institution": [{"name": "bioRxiv"}],
        "publisher": "openRxiv",
        "published": {"date-parts": [[2026, 2, 22]]},
        "DOI": "10.64898/2026.02.22.707219",
        "URL": "https://www.biorxiv.org/content/10.64898/2026.02.22.707219v1",
    }


def test_journal_article_path():
    item = crossref_to_zotero(_journal_crossref(), _COLLECTION, _TAGS)
    assert item["itemType"] == "journalArticle"
    assert item["publicationTitle"] == "Bioinformatics Advances"
    assert item["volume"] == "4"
    assert item["issue"] == "1"
    assert item["pages"] == "vbae080"
    assert item["ISSN"] == "2635-0041"
    assert item["date"] == "2024-6-15"
    assert item["abstractNote"] == "Cancer abstract."


def test_preprint_path_no_crash_on_empty_container_title():
    item = crossref_to_zotero(_biorxiv_crossref(), _COLLECTION, _TAGS)
    assert item["itemType"] == "preprint"
    assert item["repository"] == "bioRxiv"
    assert item["DOI"] == "10.64898/2026.02.22.707219"
    assert item["archiveID"] == "10.64898/2026.02.22.707219"
    assert "publicationTitle" not in item
    assert "ISSN" not in item


def test_preprint_with_abstract_populates_abstract_note():
    data = _biorxiv_crossref()
    data["abstract"] = "<jats:p>Preprint abstract.</jats:p>"
    item = crossref_to_zotero(data, _COLLECTION, _TAGS)
    assert item["itemType"] == "preprint"
    assert item["abstractNote"] == "Preprint abstract."


def test_preprint_with_missing_institution_falls_back_to_empty():
    data = _biorxiv_crossref()
    data.pop("institution")
    item = crossref_to_zotero(data, _COLLECTION, _TAGS)
    assert item["itemType"] == "preprint"
    assert item["repository"] == ""


def test_first_or_empty_handles_empty_list_and_none():
    assert _first_or_empty([]) == ""
    assert _first_or_empty(None) == ""
    assert _first_or_empty(["x"]) == "x"
    assert _first_or_empty(["x", "y"]) == "x"


def test_is_preprint_detects_both_signals():
    assert _is_preprint({"type": "posted-content", "subtype": "preprint"})
    assert _is_preprint({"type": "posted-content"})
    assert _is_preprint({"subtype": "preprint"})
    assert not _is_preprint({"type": "journal-article"})
    assert not _is_preprint({})


def test_journal_with_empty_container_title_does_not_crash():
    """Defensive: a malformed journal-article record with empty container-title
    should not crash even though preprint detection won't fire."""
    data = _journal_crossref()
    data["container-title"] = []
    item = crossref_to_zotero(data, _COLLECTION, _TAGS)
    assert item["itemType"] == "journalArticle"
    assert item["publicationTitle"] == ""
