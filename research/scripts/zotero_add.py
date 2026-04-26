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

    # ISSN: prefer print ISSN (index 0)
    issn_list = data.get("ISSN", [])
    issn = issn_list[0] if issn_list else ""

    # Abstract: CrossRef first (strips JATS tags), fall back to PubMed
    raw_abstract = data.get("abstract", "")
    abstract = re.sub(r"<[^>]+>", "", raw_abstract).strip() or pubmed_abstract or ""

    item = {
        "itemType": "journalArticle",
        "title": data.get("title", [""])[0],
        "creators": authors,
        "publicationTitle": data.get("container-title", [""])[0],
        "volume": data.get("volume", ""),
        "issue": data.get("issue", ""),
        "pages": data.get("page", ""),
        "date": date,
        "DOI": data.get("DOI", ""),
        "url": data.get("URL", ""),
        "ISSN": issn,
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
    try:
        data = fetch_crossref(args.doi)
    except urllib.error.HTTPError as e:
        sys.exit(f"CrossRef lookup failed: {e.code} {e.reason}")

    print("Fetching supplementary data from PubMed...")
    pubmed_date, pubmed_abstract, pmid = fetch_pubmed(args.doi)

    tags = DEFAULT_TAGS + args.tags
    item = crossref_to_zotero(data, ZOTERO_COLLECTION, tags, pubmed_date, pubmed_abstract, pmid)

    print(f"Title: {item['title']}")
    print(f"Authors: {len(item['creators'])} authors")
    print(f"Journal: {item['publicationTitle']} {item['volume']}({item['issue']}): {item['pages']}, {item['date']}")
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
