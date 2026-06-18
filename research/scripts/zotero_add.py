#!/usr/bin/env python3
"""Add a paper to the Splice Neoepitope Pipeline Zotero collection by DOI.

Usage:
    python scripts/zotero_add.py <DOI> [--tags tag1 tag2 ...]

Reads ZOTERO_USER_ID and ZOTERO_API_KEY from the environment (or .env file).
Fetches metadata from CrossRef, converts to Zotero JSON, and posts to the
'Splice Neoepitope Pipeline' collection (Z38GTJNW).
"""

import argparse
import json
import os
import sys
import urllib.request
import urllib.parse
import urllib.error

ZOTERO_COLLECTION = "Z38GTJNW"
DEFAULT_TAGS = ["morning-reading", "splice-neoepitope-pipeline"]


def load_env(path=".env"):
    if not os.path.exists(path):
        return
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#") and "=" in line:
                key, _, value = line.partition("=")
                os.environ.setdefault(key.strip(), value.strip())


def fetch_pubmed(doi):
    """Return (full_date_str, abstract, pmid) from PubMed, or (None, None, None) if not found."""
    import xml.etree.ElementTree as ET

    encoded = urllib.parse.quote(doi, safe="")
    search_url = (
        f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        f"?db=pubmed&term={encoded}[doi]&retmode=json"
    )
    try:
        with urllib.request.urlopen(search_url) as resp:
            ids = json.loads(resp.read()).get("esearchresult", {}).get("idlist", [])
        if not ids:
            return None, None, None
        pmid = ids[0]

        fetch_url = (
            f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
            f"?db=pubmed&id={pmid}&rettype=abstract&retmode=xml"
        )
        with urllib.request.urlopen(fetch_url) as resp:
            root = ET.fromstring(resp.read())

        article = root.find(".//PubmedArticle/MedlineCitation/Article")
        if article is None:
            return None, None, pmid

        pub_date = article.find("Journal/JournalIssue/PubDate")
        date_str = None
        if pub_date is not None:
            year = pub_date.findtext("Year", "")
            month = pub_date.findtext("Month", "")
            day = pub_date.findtext("Day", "")
            month_map = {
                "Jan": "1", "Feb": "2", "Mar": "3", "Apr": "4",
                "May": "5", "Jun": "6", "Jul": "7", "Aug": "8",
                "Sep": "9", "Oct": "10", "Nov": "11", "Dec": "12",
            }
            month = month_map.get(month, month)
            date_str = "-".join(p for p in [year, month, day] if p) or None

        abstract_texts = article.findall(".//AbstractText")
        abstract = " ".join(
            (f"{el.get('Label')}: " if el.get("Label") else "") + (el.text or "")
            for el in abstract_texts
        ).strip() or None

        return date_str, abstract, pmid

    except Exception:
        return None, None, None


def fetch_crossref(doi):
    encoded = urllib.parse.quote(doi, safe="")
    url = f"https://api.crossref.org/works/{encoded}"
    req = urllib.request.Request(url, headers={"User-Agent": "splice-neoepitope-pipeline/1.0"})
    with urllib.request.urlopen(req) as resp:
        return json.loads(resp.read())["message"]


def _first_or_empty(value):
    """Return value[0] if value is a non-empty list, else ''. Tolerates None."""
    if isinstance(value, list) and value:
        return value[0]
    return ""


def _is_preprint(data):
    """CrossRef marks preprints with type='posted-content' and/or subtype='preprint'."""
    return data.get("subtype") == "preprint" or data.get("type") == "posted-content"


def crossref_to_zotero(data, collection, tags, pubmed_date=None, pubmed_abstract=None, pmid=None):
    import re

    authors = [
        {
            "creatorType": "author",
            "firstName": a.get("given", ""),
            "lastName": a.get("family", ""),
        }
        for a in data.get("author", [])
    ]

    # Prefer published-print (journal issue date) over published-online
    date_parts = (
        data.get("published-print", data.get("published", {}))
        .get("date-parts", [[""]])
    )
    crossref_date = "-".join(str(p) for p in date_parts[0]) if date_parts else ""
    # Use PubMed date if it has more precision (more components)
    date = pubmed_date if pubmed_date and pubmed_date.count("-") > crossref_date.count("-") else crossref_date

    # Abstract: CrossRef first (strips JATS tags), fall back to PubMed
    raw_abstract = data.get("abstract", "")
    abstract = re.sub(r"<[^>]+>", "", raw_abstract).strip() or pubmed_abstract or ""

    title = _first_or_empty(data.get("title"))

    if _is_preprint(data):
        # Mirror bioRxiv's own .ris export: TY=JOUR, JF=bioRxiv (publicationTitle).
        # journalArticle renders correctly across CSL citation styles; the host
        # (bioRxiv, medRxiv, …) is taken from CrossRef's institution[].name.
        institutions = data.get("institution") or []
        publication_title = institutions[0].get("name", "") if institutions else ""
        item = {
            "itemType": "journalArticle",
            "title": title,
            "creators": authors,
            "publicationTitle": publication_title,
            "date": date,
            "DOI": data.get("DOI", ""),
            "url": data.get("URL", ""),
            "PMID": pmid or "",
            "collections": [collection],
            "tags": [{"tag": t} for t in tags],
        }
    else:
        item = {
            "itemType": "journalArticle",
            "title": title,
            "creators": authors,
            "publicationTitle": _first_or_empty(data.get("container-title")),
            "volume": data.get("volume", ""),
            "issue": data.get("issue", ""),
            "pages": data.get("page", ""),
            "date": date,
            "DOI": data.get("DOI", ""),
            "url": data.get("URL", ""),
            "ISSN": _first_or_empty(data.get("ISSN")),
            "PMID": pmid or "",
            "collections": [collection],
            "tags": [{"tag": t} for t in tags],
        }

    if abstract:
        item["abstractNote"] = abstract
    return item


def fetch_datacite(doi):
    """DataCite fallback for DOIs CrossRef doesn't know (arXiv 10.48550, dataset DOIs).

    Returns the JSON:API `data.attributes` object — mirror of fetch_crossref()
    returning CrossRef's `message`. arXiv mints through DataCite, so a CrossRef
    404 on an arXiv DOI is resolved here instead.
    """
    encoded = urllib.parse.quote(doi, safe="")
    url = f"https://api.datacite.org/dois/{encoded}"
    req = urllib.request.Request(url, headers={"User-Agent": "splice-neoepitope-pipeline/1.0"})
    with urllib.request.urlopen(req) as resp:
        return json.loads(resp.read())["data"]["attributes"]


def _datacite_publisher(data):
    """DataCite `publisher` may be a plain string or a {'name': ...} object."""
    publisher = data.get("publisher", "")
    if isinstance(publisher, dict):
        return publisher.get("name", "")
    return publisher or ""


def _datacite_is_preprint(data):
    """arXiv (publisher) or an explicit Preprint resourceTypeGeneral → preprint path.

    Mirrors _is_preprint()'s role for CrossRef; routes arXiv through the same
    journalArticle/publicationTitle-as-host treatment bioRxiv gets.
    """
    if _datacite_publisher(data).lower() == "arxiv":
        return True
    resource_type = (data.get("types") or {}).get("resourceTypeGeneral", "")
    return resource_type.lower() == "preprint"


def _datacite_creators(creators):
    """Map DataCite creators[] → Zotero author dicts (firstName/lastName).

    Personal creators carry givenName/familyName; name-only personal creators
    carry just 'name', often as 'Family, Given' — split on the comma. An
    Organizational creator's name is taken whole (no comma split), since a comma
    in an org name ('Chen, Wang & Associates') is not a Family,Given separator.
    nameType is optional in DataCite, so only an explicit 'Organizational'
    suppresses the split — absent/Personal still splits.
    """
    authors = []
    for c in creators or []:
        given = c.get("givenName", "")
        family = c.get("familyName", "")
        if not given and not family:
            name = c.get("name", "")
            if c.get("nameType") == "Organizational":
                family = name
            elif "," in name:
                family_part, _, given_part = name.partition(",")
                family = family_part.strip()
                given = given_part.strip()
            else:
                family = name
        authors.append({
            "creatorType": "author",
            "firstName": given,
            "lastName": family,
        })
    return authors


def _datacite_date(data):
    """Return the most-precise of the Issued/Submitted dates, then publicationYear.

    DataCite dates[] carry a dateType; arXiv records expose BOTH a year-only
    Issued ("2025") and a full-precision Submitted timestamp ("2025-12-06T...",
    the v1 submission date — verified live on 10.48550/arXiv.2512.06592). A fixed
    Issued-before-Submitted order would return the year-only Issued and drop the
    real date, so pick whichever entry carries the most precision instead.
    Timestamps are trimmed to the date portion (YYYY-MM-DD), matching the
    CrossRef path's `-`-joined date format. Falls back to the year-only
    publicationYear when no Issued/Submitted entry exists.
    """
    best = ""
    best_precision = -1
    for d in data.get("dates") or []:
        if d.get("dateType") in ("Issued", "Submitted") and d.get("date"):
            candidate = str(d["date"]).split("T")[0]
            precision = candidate.count("-")
            if precision > best_precision:
                best = candidate
                best_precision = precision
    if best:
        return best
    year = data.get("publicationYear")
    return str(year) if year else ""


def _datacite_abstract(data):
    for d in data.get("descriptions") or []:
        if d.get("descriptionType") == "Abstract" and d.get("description"):
            return d["description"].strip()
    return ""


def datacite_to_zotero(data, collection, tags, pubmed_date=None, pubmed_abstract=None, pmid=None):
    """Map a DataCite `data.attributes` record → the same Zotero item shape
    crossref_to_zotero() produces. arXiv routes through the preprint path
    (itemType=journalArticle, publicationTitle=<publisher>), consistent with
    how the CrossRef path treats bioRxiv."""
    authors = _datacite_creators(data.get("creators"))
    title = _first_or_empty([t.get("title", "") for t in data.get("titles") or []])
    publisher = _datacite_publisher(data)

    # Prefer DataCite date; let a more-precise PubMed date win if one exists.
    datacite_date = _datacite_date(data)
    date = pubmed_date if pubmed_date and pubmed_date.count("-") > datacite_date.count("-") else datacite_date

    abstract = _datacite_abstract(data) or pubmed_abstract or ""

    if _datacite_is_preprint(data):
        # Host (arXiv, …) goes in publicationTitle, mirroring the bioRxiv path.
        item = {
            "itemType": "journalArticle",
            "title": title,
            "creators": authors,
            "publicationTitle": publisher,
            "date": date,
            "DOI": data.get("doi", ""),
            "url": data.get("url", ""),
            "PMID": pmid or "",
            "collections": [collection],
            "tags": [{"tag": t} for t in tags],
        }
    else:
        # Non-preprint DataCite (datasets, repository DOIs): publisher (or the
        # container title, when present) acts as the publication title.
        container = data.get("container") or {}
        publication_title = container.get("title") or publisher
        item = {
            "itemType": "journalArticle",
            "title": title,
            "creators": authors,
            "publicationTitle": publication_title,
            "date": date,
            "DOI": data.get("doi", ""),
            "url": data.get("url", ""),
            "PMID": pmid or "",
            "collections": [collection],
            "tags": [{"tag": t} for t in tags],
        }

    if abstract:
        item["abstractNote"] = abstract
    return item


def post_to_zotero(item, user_id, api_key):
    url = f"https://api.zotero.org/users/{user_id}/items"
    payload = json.dumps([item]).encode()
    req = urllib.request.Request(
        url,
        data=payload,
        headers={
            "Zotero-API-Key": api_key,
            "Content-Type": "application/json",
        },
        method="POST",
    )
    with urllib.request.urlopen(req) as resp:
        return json.loads(resp.read())


def update_note(item_key, note_content, user_id, api_key):
    """Update (or create) the first note child of item_key."""
    base = f"https://api.zotero.org/users/{user_id}/items"
    headers = {"Zotero-API-Key": api_key, "Content-Type": "application/json"}

    # Fetch existing note children
    req = urllib.request.Request(f"{base}/{item_key}/children", headers=headers)
    with urllib.request.urlopen(req) as resp:
        children = json.loads(resp.read())

    notes = [c for c in children if c["data"]["itemType"] == "note"]

    if notes:
        note = notes[0]
        note_key = note["key"]
        version = note["version"]
        payload = json.dumps({
            "key": note_key,
            "version": version,
            "parentItem": item_key,
            "itemType": "note",
            "note": note_content,
            "tags": note["data"].get("tags", []),
            "relations": {},
        }).encode()
        req = urllib.request.Request(
            f"{base}/{note_key}",
            data=payload,
            headers={**headers, "If-Unmodified-Since-Version": str(version)},
            method="PUT",
        )
        with urllib.request.urlopen(req) as resp:
            resp.read()
        print(f"Note updated on item {item_key} (note key: {note_key}).")
    else:
        # No existing note — create one
        result = post_to_zotero(
            {"itemType": "note", "parentItem": item_key, "note": note_content, "tags": [], "relations": {}},
            user_id, api_key,
        )
        if result.get("successful"):
            note_key = list(result["successful"].values())[0]["key"]
            print(f"Note created on item {item_key} (note key: {note_key}).")
        else:
            sys.exit(f"Failed to create note: {result.get('failed')}")


def main():
    parser = argparse.ArgumentParser(description="Add a paper to Zotero by DOI.")
    parser.add_argument("doi", nargs="?", help="DOI of the paper to add")
    parser.add_argument("--tags", nargs="*", default=[], help="Additional tags")
    parser.add_argument("--note", help="Relevance note to attach to the Zotero entry")
    parser.add_argument("--update-note", metavar="ITEM_KEY", help="Update the note on an existing Zotero item (skips DOI add)")
    parser.add_argument("--dry-run", action="store_true", help="Print item without posting to Zotero")
    args = parser.parse_args()

    if not args.update_note and not args.doi:
        sys.exit("Error: doi is required when not using --update-note.")
    if args.update_note and not args.note:
        sys.exit("Error: --update-note requires --note.")

    load_env()

    user_id = os.environ.get("ZOTERO_USER_ID")
    api_key = os.environ.get("ZOTERO_API_KEY")
    if not user_id or not api_key:
        sys.exit("Error: ZOTERO_USER_ID and ZOTERO_API_KEY must be set in .env or environment.")

    if args.update_note:
        update_note(args.update_note, args.note, user_id, api_key)
        return

    print(f"Fetching metadata for DOI: {args.doi}")
    source = "crossref"
    try:
        data = fetch_crossref(args.doi)
    except urllib.error.HTTPError as e:
        if e.code == 404:
            # arXiv / DataCite-minted DOIs 404 on CrossRef — try DataCite.
            print("Not found on CrossRef; trying DataCite (arXiv / DataCite-minted DOIs)...")
            try:
                data = fetch_datacite(args.doi)
                source = "datacite"
            except urllib.error.HTTPError as e2:
                sys.exit(f"DOI not found on CrossRef or DataCite: {e2.code} {e2.reason}")
            except urllib.error.URLError as e2:
                # DNS / connection / TLS failure reaching DataCite — clean message,
                # not a raw traceback (Issue #702). HTTPError is caught above first.
                sys.exit(f"Network error reaching DataCite: {e2.reason}")
        else:
            sys.exit(f"CrossRef lookup failed: {e.code} {e.reason}")
    except urllib.error.URLError as e:
        # A bare URLError (DNS / connection / TLS) is not a 404 — do not fall back;
        # exit with a clear message instead of propagating a traceback (Issue #702).
        # Ordered after `except HTTPError` since HTTPError subclasses URLError.
        sys.exit(f"Network error reaching CrossRef: {e.reason}")

    print("Fetching supplementary data from PubMed...")
    pubmed_date, pubmed_abstract, pmid = fetch_pubmed(args.doi)

    tags = DEFAULT_TAGS + args.tags
    if source == "datacite":
        item = datacite_to_zotero(data, ZOTERO_COLLECTION, tags, pubmed_date, pubmed_abstract, pmid)
        is_preprint = _datacite_is_preprint(data)
    else:
        item = crossref_to_zotero(data, ZOTERO_COLLECTION, tags, pubmed_date, pubmed_abstract, pmid)
        is_preprint = _is_preprint(data)

    print(f"Title: {item['title']}")
    print(f"Authors: {len(item['creators'])} authors")
    if is_preprint:
        # itemType is journalArticle for journals and preprints alike; is_preprint
        # (set per-source above) drives the display label, not the item shape.
        print(f"Preprint ({item['publicationTitle'] or '—'}): {item['date']}")
    else:
        # .get() defaults keep this safe for the DataCite non-preprint path,
        # which omits volume/issue/pages/ISSN (only CrossRef journals carry them).
        print(f"Journal: {item['publicationTitle']} {item.get('volume', '')}({item.get('issue', '')}): {item.get('pages', '')}, {item['date']}")
        print(f"ISSN: {item.get('ISSN', '—')}")
    print(f"Abstract: {'yes' if item.get('abstractNote') else 'not available'}")
    print(f"Tags: {[t['tag'] for t in item['tags']]}")

    if args.dry_run:
        print("\n[dry-run] Item JSON:")
        print(json.dumps(item, indent=2))
        if args.note:
            print("\n[dry-run] Note:")
            print(args.note)
        return

    print("Posting to Zotero...")
    result = post_to_zotero(item, user_id, api_key)
    successful = result.get("successful", {})
    if not successful:
        failed = result.get("failed", {})
        sys.exit(f"Failed to add item: {failed}")

    key = list(successful.values())[0]["key"]
    print(f"Added successfully. Zotero key: {key}")

    if args.note:
        note_item = {
            "itemType": "note",
            "parentItem": key,
            "note": f"<p>{args.note}</p>",
            "tags": [],
            "collections": [],
            "relations": {},
        }
        note_result = post_to_zotero(note_item, user_id, api_key)
        if note_result.get("successful"):
            print("Note attached.")
        else:
            print(f"Warning: note failed to attach: {note_result.get('failed')}")


if __name__ == "__main__":
    main()
