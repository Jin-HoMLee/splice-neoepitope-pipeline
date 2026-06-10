"""Tests for zotero_add — crossref/datacite mappers + CrossRef-404 DataCite fallback."""

import io
import sys
import urllib.error
from contextlib import redirect_stdout
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

import zotero_add
from zotero_add import (
    crossref_to_zotero,
    datacite_to_zotero,
    _datacite_is_preprint,
    _first_or_empty,
    _is_preprint,
)


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


def test_biorxiv_preprint_matches_native_ris_export():
    """bioRxiv's own .ris export uses TY=JOUR + JF=bioRxiv. Verify our output
    mirrors that convention: itemType=journalArticle, publicationTitle=<institution>."""
    item = crossref_to_zotero(_biorxiv_crossref(), _COLLECTION, _TAGS)
    assert item["itemType"] == "journalArticle"
    assert item["publicationTitle"] == "bioRxiv"
    assert item["DOI"] == "10.64898/2026.02.22.707219"
    # Legacy preprint-itemType fields removed
    assert "repository" not in item
    assert "archiveID" not in item
    # CrossRef preprints don't carry ISSN
    assert "ISSN" not in item


def test_preprint_with_abstract_populates_abstract_note():
    data = _biorxiv_crossref()
    data["abstract"] = "<jats:p>Preprint abstract.</jats:p>"
    item = crossref_to_zotero(data, _COLLECTION, _TAGS)
    assert item["itemType"] == "journalArticle"
    assert item["abstractNote"] == "Preprint abstract."


def test_preprint_with_missing_institution_falls_back_to_empty():
    data = _biorxiv_crossref()
    data.pop("institution")
    item = crossref_to_zotero(data, _COLLECTION, _TAGS)
    assert item["itemType"] == "journalArticle"
    assert item["publicationTitle"] == ""


def test_medrxiv_preprint_uses_medrxiv_as_publication_title():
    """Same code path as bioRxiv — verify the institution label propagates verbatim."""
    data = _biorxiv_crossref()
    data["institution"] = [{"name": "medRxiv"}]
    item = crossref_to_zotero(data, _COLLECTION, _TAGS)
    assert item["itemType"] == "journalArticle"
    assert item["publicationTitle"] == "medRxiv"


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


# --- DataCite fallback (arXiv / DataCite-minted DOIs) ----------------------
#
# NOTE: built from the documented DataCite JSON:API v2 schema (`data.attributes`)
# + the field list in Issue #641. The single permitted live probe of
# api.datacite.org was blocked by the sandbox network allowlist ("Host not in
# allowlist"), so the live-DOI dry-run AC is deferred to human verification.


def _arxiv_datacite():
    """Real-shaped DataCite `data.attributes` for an arXiv DOI.

    publisher='arXiv', personal creators carry givenName/familyName,
    titles[].title, descriptions[] with an Abstract entry, dates[]
    (Submitted/Issued), publicationYear, types.resourceTypeGeneral='Preprint'.
    """
    return {
        "doi": "10.48550/arxiv.2512.06592",
        "identifiers": [{"identifier": "2512.06592", "identifierType": "arXiv"}],
        "creators": [
            {"name": "King, Alex", "nameType": "Personal",
             "givenName": "Alex", "familyName": "King", "affiliation": []},
            {"name": "Doe, Jane", "nameType": "Personal",
             "givenName": "Jane", "familyName": "Doe", "affiliation": []},
        ],
        "titles": [{"title": "On fine-tuning Boltz-2 for protein-protein affinity prediction"}],
        "publisher": "arXiv",
        "container": {},
        "publicationYear": 2025,
        # Verified live shape: Issued is year-only; Submitted carries the
        # full v1 submission timestamp (the real, more-precise date).
        "dates": [
            {"date": "2025-12-06T23:07:10Z", "dateType": "Submitted", "dateInformation": "v1"},
            {"date": "2025", "dateType": "Issued"},
        ],
        "types": {
            "ris": "GEN", "bibtex": "misc", "citeproc": "article-journal",
            "schemaOrg": "CreativeWork", "resourceType": "Article",
            "resourceTypeGeneral": "Preprint",
        },
        "url": "https://arxiv.org/abs/2512.06592",
        "descriptions": [
            {"description": "We fine-tune Boltz-2 for affinity.", "descriptionType": "Abstract"},
        ],
    }


def test_datacite_arxiv_maps_to_preprint_zotero_item():
    item = datacite_to_zotero(_arxiv_datacite(), _COLLECTION, _TAGS)
    assert item["itemType"] == "journalArticle"
    assert item["title"] == "On fine-tuning Boltz-2 for protein-protein affinity prediction"
    assert item["publicationTitle"] == "arXiv"
    assert item["DOI"] == "10.48550/arxiv.2512.06592"
    assert item["url"] == "https://arxiv.org/abs/2512.06592"
    # full-precision Submitted (time-stripped) wins over the year-only Issued
    assert item["date"] == "2025-12-06"
    # creators mirror the crossref mapper's firstName/lastName shape
    assert item["creators"][0] == {
        "creatorType": "author", "firstName": "Alex", "lastName": "King",
    }
    assert len(item["creators"]) == 2
    assert item["abstractNote"] == "We fine-tune Boltz-2 for affinity."
    assert {"tag": "morning-reading"} in item["tags"]
    assert item["collections"] == [_COLLECTION]
    # preprint path: no journal-only fields leak in
    assert "ISSN" not in item


def test_datacite_is_preprint_for_arxiv_and_preprint_resource_type():
    assert _datacite_is_preprint(_arxiv_datacite())
    assert _datacite_is_preprint({"types": {"resourceTypeGeneral": "Preprint"}})
    assert not _datacite_is_preprint({"publisher": "Zenodo", "types": {"resourceTypeGeneral": "Dataset"}})


def test_datacite_publisher_as_object_is_handled():
    """Newer DataCite records may nest publisher as {'name': ...}."""
    data = _arxiv_datacite()
    data["publisher"] = {"name": "arXiv"}
    item = datacite_to_zotero(data, _COLLECTION, _TAGS)
    assert item["publicationTitle"] == "arXiv"


def test_datacite_date_falls_back_to_publication_year():
    data = _arxiv_datacite()
    data["dates"] = []
    item = datacite_to_zotero(data, _COLLECTION, _TAGS)
    assert item["date"] == "2025"


def test_datacite_prefers_full_date_over_year_issued():
    """Real arXiv shape (verified live on 10.48550/arXiv.2512.06592): a year-only
    Issued alongside a full-precision Submitted timestamp. Pick the most-precise
    date (Submitted, time-stripped to YYYY-MM-DD) regardless of list order — not
    the year-only Issued that happens to come first."""
    data = _arxiv_datacite()
    data["dates"] = [
        {"date": "2025", "dateType": "Issued"},
        {"date": "2025-12-06T23:07:10Z", "dateType": "Submitted", "dateInformation": "v1"},
    ]
    item = datacite_to_zotero(data, _COLLECTION, _TAGS)
    assert item["date"] == "2025-12-06"


def test_datacite_name_only_creator_splits_family_given():
    """Organizational / name-only creators ('Family, Given') still parse."""
    data = _arxiv_datacite()
    data["creators"] = [{"name": "Smith, John", "nameType": "Personal"}]
    item = datacite_to_zotero(data, _COLLECTION, _TAGS)
    assert item["creators"][0] == {
        "creatorType": "author", "firstName": "John", "lastName": "Smith",
    }


def test_datacite_non_preprint_dataset_uses_journal_path():
    data = {
        "doi": "10.5281/zenodo.123",
        "creators": [{"givenName": "A", "familyName": "B", "nameType": "Personal"}],
        "titles": [{"title": "A dataset"}],
        "publisher": "Zenodo",
        "publicationYear": 2024,
        "dates": [{"date": "2024-01-01", "dateType": "Issued"}],
        "types": {"resourceTypeGeneral": "Dataset"},
        "url": "https://zenodo.org/records/123",
        "descriptions": [],
    }
    item = datacite_to_zotero(data, _COLLECTION, _TAGS)
    assert item["itemType"] == "journalArticle"
    assert item["publicationTitle"] == "Zenodo"
    assert item["date"] == "2024-01-01"


def _make_http_404(doi):
    raise urllib.error.HTTPError("https://api.crossref.org/works/x", 404, "Not Found", {}, None)


def test_crossref_404_falls_back_to_datacite(monkeypatch):
    captured = {}

    def _fake_datacite(doi):
        captured["doi"] = doi
        return _arxiv_datacite()

    monkeypatch.setattr(zotero_add, "fetch_crossref", _make_http_404)
    monkeypatch.setattr(zotero_add, "fetch_datacite", _fake_datacite)
    monkeypatch.setattr(zotero_add, "fetch_pubmed", lambda doi: (None, None, None))
    monkeypatch.setenv("ZOTERO_USER_ID", "0000")
    monkeypatch.setenv("ZOTERO_API_KEY", "key")
    monkeypatch.setattr(sys, "argv",
                        ["zotero_add.py", "10.48550/arXiv.2512.06592", "--dry-run"])

    buf = io.StringIO()
    with redirect_stdout(buf):
        zotero_add.main()
    out = buf.getvalue()

    assert captured["doi"] == "10.48550/arXiv.2512.06592"
    assert '"publicationTitle": "arXiv"' in out
    assert "On fine-tuning Boltz-2" in out


def test_unresolvable_doi_exits_when_both_registries_404(monkeypatch):
    def _datacite_404(doi):
        raise urllib.error.HTTPError("https://api.datacite.org/dois/x", 404, "Not Found", {}, None)

    monkeypatch.setattr(zotero_add, "fetch_crossref", _make_http_404)
    monkeypatch.setattr(zotero_add, "fetch_datacite", _datacite_404)
    monkeypatch.setattr(zotero_add, "fetch_pubmed", lambda doi: (None, None, None))
    monkeypatch.setenv("ZOTERO_USER_ID", "0000")
    monkeypatch.setenv("ZOTERO_API_KEY", "key")
    monkeypatch.setattr(sys, "argv", ["zotero_add.py", "10.9999/nonexistent", "--dry-run"])

    import pytest
    with pytest.raises(SystemExit):
        zotero_add.main()


def test_crossref_non_404_http_error_still_exits_without_datacite(monkeypatch):
    """A 500 from CrossRef is not a 'not found' — must not silently fall back."""
    def _crossref_500(doi):
        raise urllib.error.HTTPError("https://api.crossref.org/works/x", 500, "Server Error", {}, None)

    def _datacite_should_not_be_called(doi):
        raise AssertionError("DataCite must not be queried on a non-404 CrossRef error")

    monkeypatch.setattr(zotero_add, "fetch_crossref", _crossref_500)
    monkeypatch.setattr(zotero_add, "fetch_datacite", _datacite_should_not_be_called)
    monkeypatch.setattr(zotero_add, "fetch_pubmed", lambda doi: (None, None, None))
    monkeypatch.setenv("ZOTERO_USER_ID", "0000")
    monkeypatch.setenv("ZOTERO_API_KEY", "key")
    monkeypatch.setattr(sys, "argv", ["zotero_add.py", "10.1093/x", "--dry-run"])

    import pytest
    with pytest.raises(SystemExit):
        zotero_add.main()
