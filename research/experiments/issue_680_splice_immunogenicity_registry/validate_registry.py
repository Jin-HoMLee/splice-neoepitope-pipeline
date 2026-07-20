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
                                ZOTERO_COLLECTION, ZOTERO_ITEMTYPE_TO_VENUE, PREPRINT_MARKERS,
                                PEPTIDE_STATUSES, PEPTIDE_STATUS_PRESENT, PEPTIDE_NULL_STATUSES,
                                JUNCTION_EVIDENCE_GRADES, SCORABLE_TIER, TIER_ALLOWED_LABELS,
                                IN_VIVO_MODELS, IN_VIVO_NONE, IN_VIVO_MARKERS)
from registry_dedup import duplicate_keys, row_identity

HERE = Path(__file__).resolve().parent
REGISTRY = HERE / "registry.tsv"

REQUIRED_NEW_COLS = ["evidence_strength", "label_rationale", "junction_id", "junction_mapping_grade",
                     "assay_context", "venue_type", "peptide_status", "in_vivo_model"]


def violations(df: pd.DataFrame) -> list[str]:
    out = []
    for col in REQUIRED_NEW_COLS:
        if col not in df.columns:
            out.append(f"missing required column: {col}")
    if out:
        return out  # schema not yet present; stop here

    for i, r in df.iterrows():
        pep = str(r["peptide"]).strip()
        jid = str(r["junction_id"]).strip()
        # name the row by its coalesced identity: a peptide-null row is nameable by
        # its junction, and a row with neither is named as such rather than blank.
        try:
            rid = f"row {i} ({row_identity(r)})"
        except ValueError:
            rid = f"row {i} (<no identity>)"
            out.append(f"{rid}: at-least-one-non-null violated - a row needs a "
                       f"junction_id, a peptide, or both")
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
        # tier<->label firewall (#1178, #1233; made TWO-directional in #1237): every row's
        # `label` must be in TIER_ALLOWED_LABELS for its `tier`. This subsumes the old
        # positive-only guard AND closes its mirror gap - an `untested`-only tier
        # (presentation-prevalence) or a control tier wearing `label=negative` was caught by
        # neither the positive firewall (fired on positive) nor the negative-subtype checks
        # (which fire on soft/hard/strong/weak and are blind to the `na` strength those tiers
        # carry). A label-only consumer could then pull an untested peptide into a scored set, or
        # a never-tested row into the negative set - silent corruption of the registry's whole
        # product, in either direction. A tier absent from the table has an empty allowed set, so
        # an unknown tier is rejected here too. Turns the section-4 rule from documented into
        # enforced, both directions. (Not `!= SCORABLE_TIER`: functional-nonscorable and candidate
        # positives are legitimate without a scorable sequence; the *scored* set is
        # scorable_positive_mask.)
        allowed = TIER_ALLOWED_LABELS.get(r["tier"], set())
        if r["label"] not in allowed:
            out.append(f"{rid}: label {r['label']!r} on tier {r['tier']!r} - tier {r['tier']!r} "
                       f"may carry only {sorted(allowed)} (LABELING_SCHEME.md section 4)")
        # grade -> junction_id consistency (#1086). A `coords`/`event-id` grade asserts
        # the source published a junction identifier, so the column must carry it. The
        # converse block is deliberately gone: it used to forbid a `gene-mechanism`/
        # `none` row from EVER carrying a junction_id, which barred 64% of rows from
        # gaining one from a later authoritative recovery. The grade records what the
        # *source* published; the column records what we hold. The no-inference rule
        # (LABELING_SCHEME.md section 6) still bars deriving a junction from a peptide.
        grade = r["junction_mapping_grade"]
        if grade in JUNCTION_EVIDENCE_GRADES and not jid:
            out.append(f"{rid}: {grade} grade with empty junction_id")
        # peptide_status (#1086): a typed null, harmonized with length / tier / grade
        ps = r["peptide_status"]
        if ps not in PEPTIDE_STATUSES:
            out.append(f"{rid}: bad peptide_status {ps!r}")
        length = str(r["length"]).strip()
        if pep:
            if ps != PEPTIDE_STATUS_PRESENT:
                out.append(f"{rid}: peptide present but peptide_status {ps!r} "
                           f"(must be {PEPTIDE_STATUS_PRESENT!r})")
            if length != str(len(pep)):
                out.append(f"{rid}: peptide present but length {length!r} disagrees with "
                           f"len(peptide)={len(pep)}")
        else:
            if ps not in PEPTIDE_NULL_STATUSES:
                out.append(f"{rid}: peptide absent but peptide_status {ps!r} (must be one "
                           f"of {sorted(PEPTIDE_NULL_STATUSES)})")
            # a null peptide nulls its derivatives
            if length:
                out.append(f"{rid}: peptide absent but non-empty length {length!r}")
            if r["tier"] == SCORABLE_TIER:
                out.append(f"{rid}: peptide absent but tier {SCORABLE_TIER!r} - the "
                           f"benchmark keys on sequence, so this row cannot be scored")
            # with no sequence the junction IS the row, so it needs real junction evidence
            if grade not in JUNCTION_EVIDENCE_GRADES:
                out.append(f"{rid}: peptide-null row needs a coords/event-id grade, "
                           f"got {grade!r}")
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
        # in_vivo_model (#1120): the SETTING axis, orthogonal to assay_context's
        # T-cell-source axis. Controlled vocabulary.
        ivm = r["in_vivo_model"]
        if ivm not in IN_VIVO_MODELS:
            out.append(f"{rid}: bad in_vivo_model {ivm!r}")
        # Cross-check, both directions: the column must agree with the in-vivo marker
        # in `readout`. The derivation satisfies this by construction, so this can
        # only fire on a hand-edited registry - which is precisely what it guards
        # (same shape as the healthy_donor_ivs <-> IVS cross-check below).
        has_marker = any(m in str(r["readout"]).lower() for m in IN_VIVO_MARKERS)
        if has_marker and ivm == IN_VIVO_NONE:
            out.append(f"{rid}: readout names an in-vivo animal readout but "
                       f"in_vivo_model is {IN_VIVO_NONE!r} (derive it, don't hand-edit)")
        if not has_marker and ivm != IN_VIVO_NONE:
            out.append(f"{rid}: in_vivo_model {ivm!r} but readout names no in-vivo readout")
        # The coupling that keeps the two axes honest: `animal_syngeneic` is the ONLY
        # assay_context in which the animal is the T-cell source, and it is coherent
        # only in a syngeneic (immunocompetent) host. A xenograft host is
        # immunodeficient - it has no T cells to be the source of - so this pairing is
        # the one that catches a future attempt to re-merge the setting into the
        # source axis.
        if ac == "animal_syngeneic" and ivm != "syngeneic":
            out.append(f"{rid}: assay_context 'animal_syngeneic' requires "
                       f"in_vivo_model 'syngeneic', got {ivm!r} - an immunodeficient "
                       f"(xenograft) host cannot be the source of the responding T cells")
        # venue_type (#1001): controlled vocabulary. The out-of-vocab sentinel
        # 'unclassified' fails here, so a source absent from the venue map cannot
        # be folded venue-unmarked (esp. a preprint).
        if r["venue_type"] not in VENUE_TYPES:
            out.append(f"{rid}: bad venue_type {r['venue_type']!r} "
                       f"(source not classified in derive_venue_type.py)")
        # preprint-marker guard (#1001 review finding 1): the source-keyed derivation
        # can't distinguish a preprint from the journal version of an already-mapped
        # study, so catch it here - a source naming a preprint server must be `preprint`.
        src_l = str(r["source"]).lower()
        if r["venue_type"] != "preprint" and any(m in src_l for m in PREPRINT_MARKERS):
            out.append(f"{rid}: source names a preprint server but venue_type is "
                       f"{r['venue_type']!r} (map it as 'preprint' in labeling_constants.py)")

    # dedup on the full identity triple (#1086). Not on the coalesced identity alone:
    # one junction legitimately carries several distinct peptides.
    for key in duplicate_keys(df):
        out.append(f"duplicate identity key (junction_id, peptide, hla) = {key}")

    out.extend(source_key_violations(df))
    return out


def source_key_violations(df: pd.DataFrame) -> list[str]:
    """`source` must map 1:1 onto studies (#1106).

    `source` is documented as the stable curation key, but every consumer resolves
    it through a *lowercased substring* map (VENUE_BY_SOURCE_SUBSTR here,
    derive_assay_context.py and derive_venue_type.py likewise). So two spellings of
    one study ("IRIS" and "IRIS (Pan/Xing, PNAS)") both match the same key and derive
    identically - the split is invisible to every other check in this file, and it
    was: the registry carried 15 strings for 12 studies with nothing firing.

    It is not harmless. `source` is also the group-by key for study-level analysis,
    so a split study is silently double-counted (the #737 sparsity notebook was
    counting Xiong twice), and a *future* source whose name is a substring of another
    would mis-derive its venue and DOI outright.

    So: fail when two distinct `source` strings collide on one substring key. This
    is a registry-level invariant, not a row-level one - it cannot be expressed by
    looking at any single row, which is why it lived undetected in a file otherwise
    dense with per-row guards.
    """
    out = []
    by_key: dict[str, set[str]] = {}
    for src in sorted({str(s) for s in df["source"]}):
        src_l = src.lower()
        for substr in VENUE_BY_SOURCE_SUBSTR:
            if substr in src_l:
                by_key.setdefault(substr, set()).add(src)
                break  # first match wins - mirrors derive_venue_type.py exactly
        else:
            # Unmapped source: already fails loudly via the venue_type check (it
            # derives to the `unclassified` sentinel), so don't double-report here.
            continue

    for substr, sources in sorted(by_key.items()):
        if len(sources) > 1:
            spellings = " | ".join(repr(s) for s in sorted(sources))
            out.append(
                f"source key collision: {len(sources)} distinct `source` spellings "
                f"resolve to the same key {substr!r} -> {spellings}. "
                f"`source` is the curation key and the study group-by key; pick one "
                f"spelling per study and record the rewrite in PROVENANCE.md."
            )
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
            # this .env lacks creds; keep walking up (e.g. research/.env -> root .env)
            continue
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
    except urllib.error.HTTPError as e:
        # creds ARE present but the request was rejected (403 bad key / 404 wrong
        # collection): the operator intended to verify, so surface this loudly
        # rather than a quiet skip - but still don't hard-fail (could be transient).
        print(f"WARNING: Zotero venue cross-check could NOT run - HTTP {e.code} "
              f"(check ZOTERO_API_KEY / collection {ZOTERO_COLLECTION}). Not verified.",
              file=sys.stderr)
        return [], False
    except (urllib.error.URLError, TimeoutError, ValueError, json.JSONDecodeError) as e:
        print(f"NOTE: Zotero venue cross-check skipped (offline / fetch failed: {e}).", file=sys.stderr)
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
    print("venue_type: " + ", ".join(f"{vt}={int(ct)}" for vt, ct in venue_counts.items()))
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
