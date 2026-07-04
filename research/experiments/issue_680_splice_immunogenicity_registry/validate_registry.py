#!/usr/bin/env python3
"""Validate registry.tsv against the documented labeling scheme (Issue #735).

Run: research/.venv/bin/python validate_registry.py
Exit 0 = clean; exit 1 = violations printed to stderr.
"""
import json
import os
import sys
import urllib.error
import urllib.request
import pandas as pd
from pathlib import Path

from labeling_constants import (GRADES, STRENGTHS, ASSAY_CONTEXTS, VENUE_TYPES, EFFECTOR,
                                DETECTION, IVS_MARKER, VENUE_BY_SOURCE_SUBSTR,
                                ZOTERO_COLLECTION, ZOTERO_ITEMTYPE_TO_VENUE)

HERE = Path(__file__).resolve().parent
REGISTRY = HERE / "registry.tsv"

REQUIRED_NEW_COLS = ["evidence_strength", "label_rationale", "junction_id", "junction_mapping_grade",
                     "assay_context", "venue_type"]


def violations(df: pd.DataFrame) -> list[str]:
    out = []
    for col in REQUIRED_NEW_COLS:
        if col not in df.columns:
            out.append(f"missing required column: {col}")
    if out:
        return out  # schema not yet present; stop here

    for i, r in df.iterrows():
        rid = f"row {i} ({r['peptide']})"
        if r["junction_mapping_grade"] not in GRADES:
            out.append(f"{rid}: bad junction_mapping_grade {r['junction_mapping_grade']!r}")
        if r["evidence_strength"] not in STRENGTHS:
            out.append(f"{rid}: bad evidence_strength {r['evidence_strength']!r}")
        if not str(r["label_rationale"]).strip():
            out.append(f"{rid}: empty label_rationale")
        # consistency: strong positives must name an effector readout
        ro = str(r["readout"]).lower()
        if r["label"] == "positive" and r["evidence_strength"] == "strong":
            if not any(k in ro for k in EFFECTOR):
                out.append(f"{rid}: strong positive but readout has no effector term")
        if r["label"] == "positive" and r["evidence_strength"] == "weak":
            if any(k in ro for k in EFFECTOR):
                out.append(f"{rid}: weak positive but readout names an effector term")
            if not any(k in ro for k in DETECTION):
                out.append(f"{rid}: weak positive but readout has no detection term")
        # consistency: negative subtypes cross-check the IVS-context rule
        if r["label"] == "negative" and r["evidence_strength"] == "soft":
            if IVS_MARKER not in ro:
                out.append(f"{rid}: soft negative but readout has no IVS context")
        if r["label"] == "negative" and r["evidence_strength"] == "hard":
            if IVS_MARKER in ro:
                out.append(f"{rid}: hard negative but readout indicates IVS context (should be soft)")
        # negatives must carry a negative subtype, never a positive strength
        if r["label"] == "negative" and r["evidence_strength"] in {"strong", "weak"}:
            out.append(f"{rid}: negative row with positive-only evidence_strength {r['evidence_strength']!r}")
        # grade 'none' must record a reason in notes
        if r["junction_mapping_grade"] == "none" and not str(r["notes"]).strip():
            out.append(f"{rid}: grade 'none' with no reason in notes")
        # positive rows must not resolve to na evidence_strength
        if r["label"] == "positive" and r["evidence_strength"] == "na":
            out.append(f"{rid}: positive row resolved to na evidence_strength (needs manual review)")
        # grade <-> junction_id consistency
        grade = r["junction_mapping_grade"]
        jid = str(r["junction_id"]).strip()
        if grade in {"coords", "event-id"} and not jid:
            out.append(f"{rid}: {grade} grade with empty junction_id")
        if grade in {"gene-mechanism", "none"} and jid:
            out.append(f"{rid}: {grade} grade with non-empty junction_id")
        # assay_context (#823): controlled vocabulary + cross-checks
        ac = r["assay_context"]
        if ac not in ASSAY_CONTEXTS:
            out.append(f"{rid}: bad assay_context {ac!r}")
        # healthy_donor_ivs <-> IVS marker in readout (ties to the #735 IVS rule, both directions)
        if (IVS_MARKER in ro) != (ac == "healthy_donor_ivs"):
            out.append(f"{rid}: assay_context {ac!r} inconsistent with IVS marker in readout")
        # prevalence_only is exactly the presentation-prevalence tier
        if (ac == "prevalence_only") != (r["tier"] == "presentation-prevalence"):
            out.append(f"{rid}: prevalence_only must match the presentation-prevalence tier")
        # na is exactly the negative-control-not-splice tier (mirror of the above)
        if (ac == "na") != (r["tier"] == "negative-control-not-splice"):
            out.append(f"{rid}: na assay_context must match the negative-control-not-splice tier")
        # venue_type (#1001): controlled vocabulary. The out-of-vocab sentinel
        # 'unclassified' fails here, so a source absent from the venue map cannot
        # be folded venue-unmarked (esp. a preprint).
        if r["venue_type"] not in VENUE_TYPES:
            out.append(f"{rid}: bad venue_type {r['venue_type']!r} "
                       f"(source not classified in derive_venue_type.py)")
    return out


def decoy_violations(path: Path) -> list[str]:
    """Validate the materialized #681 presented-decoy seed (a committed deliverable).

    A missing file is a failure, not a skip: the seed is part of the #735 deliverable
    set, so its absence means the registry is incomplete rather than simply un-annotated.
    """
    if not path.exists():
        return [f"missing decoy seed file: {path.name}"]
    dd = pd.read_csv(path, sep="\t", dtype=str).fillna("")
    out = []
    for col in ("peptide", "decoy_tier"):
        if col not in dd.columns:
            out.append(f"decoy seed missing column: {col}")
    if out:
        return out
    if len(dd) != 13:
        out.append(f"expected 13 presented decoys, got {len(dd)}")
    if not (dd["decoy_tier"] == "presented-decoy").all():
        out.append("decoy_tier not all 'presented-decoy'")
    if not dd["peptide"].is_unique:
        out.append("duplicate decoy peptide")
    if not dd["peptide"].str.len().between(8, 11).all():
        out.append("decoy peptide length out of [8, 11]")
    return out


def _zotero_creds() -> tuple[str, str] | None:
    """Return (user_id, api_key) from the environment or the project-root .env,
    or None if unavailable. The cross-check is offline-optional: no creds -> skip."""
    uid, key = os.environ.get("ZOTERO_USER_ID"), os.environ.get("ZOTERO_API_KEY")
    if uid and key:
        return uid, key
    for parent in [HERE, *HERE.parents]:
        env = parent / ".env"
        if env.exists():
            vals = {}
            for raw in env.read_text().splitlines():
                line = raw.strip()
                if line.startswith("export "):
                    line = line[len("export "):]
                if line.startswith("#") or "=" not in line:
                    continue
                k, _, v = line.partition("=")
                vals[k.strip()] = v.strip().strip('"').strip("'")
            uid, key = vals.get("ZOTERO_USER_ID"), vals.get("ZOTERO_API_KEY")
            if uid and key:
                return uid, key
            return None
    return None


def _fetch_zotero_doi_itemtypes(uid: str, key: str) -> dict[str, str]:
    """Map normalized DOI -> itemType across all top-level items in the collection.
    Paginated (the collection exceeds one page). Raises on any HTTP/JSON error."""
    out, start, limit = {}, 0, 100
    while True:
        url = (f"https://api.zotero.org/users/{uid}/collections/{ZOTERO_COLLECTION}"
               f"/items/top?format=json&limit={limit}&start={start}")
        req = urllib.request.Request(url, headers={"Zotero-API-Key": key})
        with urllib.request.urlopen(req, timeout=20) as resp:
            page = json.loads(resp.read())
        if not page:
            break
        for it in page:
            d = it.get("data", {})
            doi = str(d.get("DOI", "")).strip().lower()
            if doi:
                out[doi] = d.get("itemType", "")
        if len(page) < limit:
            break
        start += limit
    return out


def zotero_venue_crosscheck(df: pd.DataFrame) -> tuple[list[str], bool]:
    """Offline-optional cross-check: assert each source's mapped venue_type agrees
    with the itemType Zotero records for that source's DOI. Returns (violations,
    ran). ran=False means the check was skipped (no creds / network error) - never
    a failure, so hermetic CI/offline runs stay green; a genuine venue<->itemType
    disagreement (creds present) is a real violation. DOI-keyed, so it sidesteps
    the parent-vs-attachment-key trap (only top-level items carry DOIs)."""
    creds = _zotero_creds()
    if creds is None:
        return [], False
    try:
        doi_types = _fetch_zotero_doi_itemtypes(*creds)
    except (urllib.error.URLError, TimeoutError, ValueError, json.JSONDecodeError) as e:
        print(f"NOTE: Zotero venue cross-check skipped (fetch failed: {e}).", file=sys.stderr)
        return [], False

    present = {s.lower() for s in df["source"].unique()}
    out = []
    for substr, spec in VENUE_BY_SOURCE_SUBSTR.items():
        if not any(substr in s for s in present):
            continue  # source not currently in the registry
        doi = spec["doi"].lower()
        itemtype = doi_types.get(doi)
        if itemtype is None:
            out.append(f"source {substr!r}: DOI {doi} not found in Zotero collection "
                       f"{ZOTERO_COLLECTION} (add the item, or fix the DOI in labeling_constants.py)")
            continue
        expected = ZOTERO_ITEMTYPE_TO_VENUE.get(itemtype)
        if expected is None:
            out.append(f"source {substr!r}: Zotero itemType {itemtype!r} has no venue_type mapping "
                       f"(extend ZOTERO_ITEMTYPE_TO_VENUE)")
        elif expected != spec["venue"]:
            out.append(f"source {substr!r}: venue_type {spec['venue']!r} disagrees with Zotero "
                       f"itemType {itemtype!r} (implies {expected!r})")
    return out, True


def main() -> int:
    df = pd.read_csv(REGISTRY, sep="\t", dtype=str).fillna("")
    v = violations(df)
    if v:
        print(f"FAIL: {len(v)} violation(s):", file=sys.stderr)
        for line in v:
            print("  -", line, file=sys.stderr)
        return 1
    print(f"PASS: {len(df)} rows valid against the labeling scheme.")

    # venue_type (#1001) advisory summary: surface the preprint-row count so a
    # non-peer-reviewed source in the ground-truth set is visible at a glance.
    # A preprint row is legitimate (it just must be marked), so this reports, it
    # does not reject.
    venue_counts = df["venue_type"].value_counts()
    n_preprint = int(venue_counts.get("preprint", 0))
    print("venue_type: " + ", ".join(f"{v}={int(n)}" for v, n in venue_counts.items()))
    if n_preprint:
        print(f"ADVISORY: {n_preprint} preprint-sourced row(s) in the set "
              f"(non-peer-reviewed; consider down-weighting in sensitivity analysis).")

    # venue_type <-> Zotero itemType cross-check (#1001). Offline-optional: runs only
    # when ZOTERO_* creds are present (skipped cleanly in hermetic CI / offline), and
    # fails on a genuine disagreement between the stored venue and Zotero's itemType.
    vx, ran = zotero_venue_crosscheck(df)
    if vx:
        print(f"FAIL: {len(vx)} venue_type/Zotero cross-check violation(s):", file=sys.stderr)
        for line in vx:
            print("  -", line, file=sys.stderr)
        return 1
    print("PASS: venue_type agrees with Zotero itemType." if ran
          else "SKIP: Zotero venue cross-check (no ZOTERO_* creds; hermetic run).")

    DECOY = HERE / "decoy_negatives" / "presented_decoys_681.tsv"
    dv = decoy_violations(DECOY)
    if dv:
        print(f"FAIL: {len(dv)} decoy violation(s):", file=sys.stderr)
        for line in dv:
            print("  -", line, file=sys.stderr)
        return 1
    print("PASS: 13 presented decoys valid.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
