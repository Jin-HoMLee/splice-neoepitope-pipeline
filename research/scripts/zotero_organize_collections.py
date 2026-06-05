#!/usr/bin/env python3
"""Reorganize the 'Splice Neoepitope Pipeline' Zotero collection (Z38GTJNW)
into 8 pipeline-stage sub-collections, ensure every item carries the project
tag, and report duplicate titles for manual merge.

Issue #675. Idempotent and dry-run-first: running with no flag mutates nothing
and prints the full plan. Pass --apply to execute creates + assignments + tags.

  python research/scripts/zotero_organize_collections.py            # dry-run (read-only)
  python research/scripts/zotero_organize_collections.py --apply    # execute

Design notes
------------
* The spine is the *pipeline stage*; manuscript-section and Issue# stay as tags.
* Items are ADDED to sub-collections; their existing parent (Z38GTJNW)
  membership is preserved, so the reorg is non-destructive and reversible.
* Genuinely cross-stage papers get membership in two sub-collections.
* Duplicate detection is REPORT-ONLY -- the Zotero web API has no faithful
  "merge items" primitive (it lives in the desktop client's Duplicate Items
  pane), and deleting a copy would drop its notes/attachments. The script names
  the duplicate + its child counts so the merge can be done safely by hand.
"""

import argparse
import json
import os
import sys
import time
import urllib.error
import urllib.request

ZOTERO_COLLECTION = "Z38GTJNW"
PROJECT_TAG = "splice-neoepitope-pipeline"
API_BASE = "https://api.zotero.org"

# Pipeline-stage spine. Numeric prefix forces pipeline order in the Zotero
# left pane (collections otherwise sort alphabetically).
FOLDERS = {
    1: "1. Splice detection & prediction",
    2: "2. Splice-neoantigen biology & discovery",
    3: "3. MHC presentation & immunogenicity",
    4: "4. TCR-pMHC structure & binding",
    5: "5. Immunopeptidomics & MS validation",
    6: "6. Reference & normal data",
    7: "7. Clinical translation & vaccines",
    8: "8. Methods & infrastructure",
}

# (distinctive lowercase title substring, [folder numbers]).  A paper in two
# folders is an intentional cross-stage spanner.  Substrings avoid non-ASCII
# (e.g. C/EBPbeta, en-dashes) so matching is robust.
MATCHERS = [
    # 1 · detection & prediction
    ("mmsplice: modular modeling", [1]),
    ("predicting splicing from primary sequence", [1]),                 # SpliceAI
    ("advancing regulatory variant effect prediction with alphagenome", [1]),
    ("benchmarking foundation models for splice site", [1]),
    ("prediction of tumor-specific splicing from somatic mutations", [1, 2]),  # splice2neo (spanner)
    # 2 · splice-neoantigen biology & discovery
    ("asneo: identification of personalized", [2]),
    ("neoantigens in cancer immunotherapy: focusing on alternative", [2]),
    ("alternative splicing: from tumorigenesis", [2]),
    ("splicemutr", [2, 1]),                                              # tool+biology (spanner)
    ("splicing neoantigen discovery with snaf", [2, 1]),                # tool+biology (spanner)
    ("alternative splicing of rcan1", [2]),                             # C/EBPb RCAN1 GBM
    ("decoding neoantigen-encoding tumor-specific transcripts", [2, 5]),  # neoTST HCC (spanner)
    ("efficient and effective identification of cancer neoantigens from tumor only", [2]),  # ENEO
    ("harnessing alternative splicing for off-the-shelf", [2, 7]),      # HCC mRNA (spanner)
    ("large-scale transcript variants dictate neoepitopes", [2]),
    ("mis-splicing-derived neoantigens and cognate", [2]),
    ("single-cell mapping of alternative splicing linked to checkpoint", [2]),
    ("targeting srsf1 improves cancer immunotherapy", [2, 7]),          # splice-modulator (spanner)
    ("transcript-targeted antigen mapping reveals the potential of postn", [2, 5]),  # POSTN (spanner)
    ("tumour-wide rna splicing aberrations generate actionable public", [2]),
    ("corest complex inhibition alters rna splicing", [2]),
    ("harnessing tumor-specific transcript diversity uncovers", [2, 5]),  # FOXA2 PDAC (spanner)
    ("identification of immunogenic neoantigens from intron retention", [2, 5]),  # CRC IR (spanner)
    # 3 · MHC presentation & immunogenicity
    ("a neoantigen fitness model predicts tumour response", [3]),
    ("mhcflurry 2.0", [3]),
    ("improve: a feature model to predict neoepitope immunogenicity", [3]),
    ("immunostruct", [3]),
    ("leveraging attention-based deep multiple instance", [3]),
    ("neoguider", [3]),
    ("cnneopp", [3]),
    ("neoprecis: enhancing immunotherapy response prediction", [3]),
    # 4 · TCR-pMHC structure & binding
    ("accurate structure prediction of biomolecular interactions with alphafold 3", [4]),
    ("chai-1", [4]),
    ("structural characterization and alphafold modeling of human t cell receptor recognition of nras", [4, 2]),  # NRAS (spanner)
    ("benchmarking tcr-pmhc structure prediction: a unified evaluation", [4]),
    ("boltz-2: towards accurate and efficient binding affinity", [4]),
    ("comparative analysis of tcr and tcr-pmhc complex structure prediction", [4]),
    ("nettcr-struc", [4]),
    ("on fine-tuning boltz-2", [4]),
    ("t cell receptor specificity landscape revealed through de novo peptide design", [4]),  # HERMES
    ("ai predicted tcr-pmhc structures differentiate immune interactions", [4]),  # (duplicate title)
    ("immset", [4]),
    ("language modeling materializes a world model of protein biology", [4]),  # ESMFold2
    ("predicting tcr-pmhc binding by reinforcement learning", [4]),
    ("binding prediction and generalization to unseen peptides", [4]),  # Structure-Based TCR-pMHC
    ("t2pmhc: a structure-informed graph neural network", [4]),
    ("tcrlens", [4]),
    # 5 · immunopeptidomics & MS validation
    ("from prediction to precision: how immunopeptidomics advances", [5]),
    ("mhc1-tip enables single-tube", [5]),
    ("non-canonical transcription and splicing shape the colorectal", [5, 2]),  # CRC immunopeptidome (spanner)
    ("sensitive detection of cancer antigens enabled by user-defined peptide libraries", [5]),
    # 6 · reference & normal data
    ("snaptron: querying splicing patterns", [6]),
    ("the gtex consortium atlas", [6]),
    ("transcriptomic and genomic testing to guide individualized treatment in chemoresistant gastric", [6]),
    ("recount3", [6]),
    # 7 · clinical translation & vaccines
    ("base editing screens map mutations affecting interferon", [7]),   # immune-evasion mechanism (loose fit)
    ("personalized rna neoantigen vaccines stimulate t cells in pancreatic", [7]),
    ("individualised neoantigen therapy mrna-4157", [7]),
    ("rna neoantigen vaccines prime long-lived cd8", [7]),
    ("adjuvant personalized multivalent neoantigen dna vaccination", [7]),
    ("bridging clinical gaps in personalized cancer neoantigen vaccines", [7]),
    ("individualized mrna vaccines evoke durable t cell immunity in adjuvant tnbc", [7]),
    ("intismeran autogene plus pembrolizumab", [7]),
    ("personalized cancer vaccines in the clinical trial pipeline", [7]),
    ("pipe dream to pipeline", [7]),
    ("prime-target neoantigen vaccination unleashes", [7]),
    ("safety-first tcr-t", [7]),
    # 8 · methods & infrastructure
    ("splicing neoepitope prediction is sensitive to methodological", [8]),
    ("snakemake hackathon 2026", [8]),
]

# Manuscript-section facet -> saved searches (virtual, auto-updating views over
# the existing manuscript-* tags). Each search OR-matches its tag plus known
# case-variants/synonyms (Zotero defaults to AND, so a joinMode=any condition is
# required). Verified API format 2026-06-05.
SAVED_SEARCHES = [
    ("Cited in INTRODUCTION", ["manuscript-INTRODUCTION", "introduction-citation"]),
    ("Cited in METHODS", ["manuscript-METHODS", "methods-citation"]),
    ("Cited in DISCUSSION",
     ["manuscript-DISCUSSION", "manuscript-discussion", "manuscript-citation-candidate"]),
]


def load_env(path=".env"):
    """Load KEY=VALUE pairs from a .env file, walking up from cwd to find it."""
    here = os.path.abspath(os.getcwd())
    while True:
        candidate = os.path.join(here, path)
        if os.path.exists(candidate):
            with open(candidate) as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith("#") and "=" in line:
                        k, _, v = line.partition("=")
                        os.environ.setdefault(k.strip(), v.strip())
            return candidate
        parent = os.path.dirname(here)
        if parent == here:
            return None
        here = parent


def api_request(method, url, key, body=None, extra_headers=None):
    headers = {"Zotero-API-Key": key, "Zotero-API-Version": "3"}
    data = None
    if body is not None:
        headers["Content-Type"] = "application/json"
        data = json.dumps(body).encode()
    if extra_headers:
        headers.update(extra_headers)
    req = urllib.request.Request(url, data=data, headers=headers, method=method)
    with urllib.request.urlopen(req) as resp:
        raw = resp.read()
        return resp.status, (json.loads(raw) if raw else None), dict(resp.headers)


def get_paginated(path, key):
    """GET all items across pages from a users/<id> sub-path."""
    out, start = [], 0
    while True:
        url = f"{API_BASE}/users/{os.environ['ZOTERO_USER_ID']}/{path}"
        url += ("&" if "?" in url else "?") + f"limit=100&start={start}"
        status, data, hdr = api_request("GET", url, key)
        out.extend(data)
        total = int(hdr.get("Total-Results", len(out)))
        start += len(data)
        if not data or start >= total:
            return out


def norm(title):
    return " ".join((title or "").lower().split())


def folders_for(title):
    t = norm(title)
    nums = []
    for sub, folder_nums in MATCHERS:
        if sub in t:
            nums.extend(folder_nums)
    return sorted(set(nums))


def plan_saved_searches(uid, key, apply):
    """Print (and, if apply, create) the manuscript-section saved searches.
    Idempotent: a search whose name already exists is left untouched."""
    existing = get_paginated("searches", key)
    have = {s["data"]["name"] for s in existing}
    print("\nSaved searches (manuscript-section facet):")
    to_make = []
    for name, tags in SAVED_SEARCHES:
        if name in have:
            print(f"  EXISTS  {name}")
        else:
            print(f"  {'CREATE' if apply else 'WILL CREATE'}  {name}  <- {' || '.join(tags)}")
            to_make.append((name, tags))
    if apply and to_make:
        body = []
        for name, tags in to_make:
            conds = [{"condition": "joinMode", "operator": "any", "value": ""}]
            conds += [{"condition": "tag", "operator": "is", "value": t} for t in tags]
            body.append({"name": name, "conditions": conds})
        status, data, _ = api_request("POST", f"{API_BASE}/users/{uid}/searches", key, body=body)
        if data.get("failed"):
            sys.exit(f"ERROR creating saved searches: {data['failed']}")
        made = data.get("successful") or data.get("success") or {}
        print(f"  -> created {len(made)} saved search(es)")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--apply", action="store_true", help="execute mutations (default: dry-run)")
    ap.add_argument("--saved-searches", action="store_true",
                    help="also create the manuscript-section saved searches (with --apply)")
    args = ap.parse_args()

    load_env()
    key = os.environ.get("ZOTERO_API_KEY")
    uid = os.environ.get("ZOTERO_USER_ID")
    if not key or not uid:
        sys.exit("ERROR: ZOTERO_API_KEY / ZOTERO_USER_ID not found (.env or env).")

    mode = "APPLY" if args.apply else "DRY-RUN"
    print(f"=== Zotero reorg [{mode}] — collection {ZOTERO_COLLECTION} ===\n")

    # 1. Existing sub-collections (idempotent: reuse by name).
    existing = get_paginated(f"collections/{ZOTERO_COLLECTION}/collections", key)
    name_to_key = {c["data"]["name"]: c["key"] for c in existing}

    print("Sub-collections:")
    to_create = []
    for n in sorted(FOLDERS):
        nm = FOLDERS[n]
        if nm in name_to_key:
            print(f"  EXISTS  {nm}  ({name_to_key[nm]})")
        else:
            print(f"  CREATE  {nm}")
            to_create.append(nm)

    if args.apply and to_create:
        url = f"{API_BASE}/users/{uid}/collections"
        body = [{"name": nm, "parentCollection": ZOTERO_COLLECTION} for nm in to_create]
        status, data, _ = api_request("POST", url, key, body=body)
        succ = data.get("success") or {}
        for idx, ck in succ.items():
            name_to_key[to_create[int(idx)]] = ck
        if data.get("failed"):
            sys.exit(f"ERROR creating collections: {data['failed']}")
        print(f"  -> created {len(succ)} collection(s)")

    folder_key = {n: name_to_key.get(FOLDERS[n]) for n in FOLDERS}

    # 2. Items.
    items = get_paginated(f"collections/{ZOTERO_COLLECTION}/items/top", key)
    print(f"\nItems: {len(items)} top-level\n")

    per_folder = {n: [] for n in FOLDERS}
    tag_plan, unmatched, by_title, spanners = [], [], {}, []
    for it in items:
        d = it["data"]
        title = d.get("title", "(no title)")
        by_title.setdefault(norm(title), []).append(it)
        cur_cols = set(d.get("collections", []))
        cur_tags = {t["tag"] for t in d.get("tags", [])}

        nums = folders_for(title)
        if not nums:
            unmatched.append(title)
        for n in nums:
            already = bool(folder_key[n]) and folder_key[n] in cur_cols  # True only on re-runs
            per_folder[n].append((title, already))
        if len(nums) > 1:
            spanners.append((title, nums))
        if PROJECT_TAG not in cur_tags:
            tag_plan.append(it)

    print("Planned assignment by sub-collection  (+ = will add, . = already filed):")
    to_add = 0
    for n in sorted(FOLDERS):
        rows = sorted(per_folder[n])
        print(f"\n  {FOLDERS[n]}  ({len(rows)})")
        for title, already in rows:
            print(f"    {'.' if already else '+'} {title[:74]}")
            to_add += 0 if already else 1
    print(f"\n  -> {to_add} membership(s) to add across {len(items)} items")

    if spanners:
        print("\nCross-stage papers (filed in two sub-collections):")
        for title, nums in spanners:
            print(f"  [{', '.join(map(str, nums))}]  {title[:66]}")

    print(f"\nTag fixes (+{PROJECT_TAG}):")
    for it in tag_plan:
        print(f"  +  {it['data'].get('title','')[:72]}")
    if not tag_plan:
        print("  (none — every item already carries the project tag)")

    print("\nDuplicate titles (REPORT-ONLY — merge by hand in the desktop client):")
    dup_found = False
    for t, group in by_title.items():
        if len(group) > 1:
            dup_found = True
            print(f"  '{group[0]['data'].get('title','')[:60]}' x{len(group)}")
            for it in group:
                k = it["key"]
                children = get_paginated(f"items/{k}/children", key)
                nnotes = sum(1 for c in children if c["data"].get("itemType") == "note")
                natt = sum(1 for c in children if c["data"].get("itemType") == "attachment")
                ntags = len(it["data"].get("tags", []))
                print(f"      {k}  v{it['version']}  tags={ntags}  notes={nnotes}  attachments={natt}")
    if not dup_found:
        print("  (none)")

    if unmatched:
        print("\n!! UNMATCHED (no folder assigned — fix MATCHERS before --apply):")
        for t in unmatched:
            print(f"  ?  {t[:72]}")

    if not args.apply:
        plan_saved_searches(uid, key, apply=False)

    # 3. Apply mutations.
    if args.apply:
        if unmatched:
            sys.exit("\nABORT: unmatched items present; refusing to apply a partial mapping.")
        print("\nApplying...")
        changed = 0
        # rebuild plans now that folder_key is populated post-create
        for it in items:
            d = it["data"]
            cur_cols = list(d.get("collections", []))
            cur_tags = list(d.get("tags", []))
            nums = folders_for(d.get("title", ""))
            target_keys = [folder_key[n] for n in nums if folder_key[n]]
            new_cols = cur_cols + [k for k in target_keys if k not in cur_cols]
            new_tags = cur_tags + (
                [{"tag": PROJECT_TAG}] if PROJECT_TAG not in {t["tag"] for t in cur_tags} else []
            )
            patch = {}
            if set(new_cols) != set(cur_cols):
                patch["collections"] = new_cols
            if len(new_tags) != len(cur_tags):
                patch["tags"] = new_tags
            if not patch:
                continue
            url = f"{API_BASE}/users/{uid}/items/{it['key']}"
            for attempt in range(3):
                try:
                    api_request("PATCH", url, key, body=patch,
                                extra_headers={"If-Unmodified-Since-Version": str(it["version"])})
                    changed += 1
                    break
                except urllib.error.HTTPError as e:
                    if e.code == 412 and attempt < 2:  # version conflict: refetch + retry
                        _, fresh, _ = api_request(
                            "GET", f"{API_BASE}/users/{uid}/items/{it['key']}", key)
                        it["version"] = fresh["version"]
                        time.sleep(1)
                        continue
                    raise
        print(f"  -> patched {changed} item(s)")
        if args.saved_searches:
            plan_saved_searches(uid, key, apply=True)
        print("\nDONE. Duplicates (if any) still need a manual merge in the desktop client.")
    else:
        print("\n(dry-run — nothing was modified. Re-run with --apply to execute.)")


if __name__ == "__main__":
    main()
