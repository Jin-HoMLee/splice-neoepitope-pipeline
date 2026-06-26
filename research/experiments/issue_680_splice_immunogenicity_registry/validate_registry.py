#!/usr/bin/env python3
"""Validate registry.tsv against the documented labeling scheme (Issue #735).

Run: research/.venv/bin/python validate_registry.py
Exit 0 = clean; exit 1 = violations printed to stderr.
"""
import sys
import pandas as pd
from pathlib import Path

from labeling_constants import GRADES, STRENGTHS, ASSAY_CONTEXTS, EFFECTOR, DETECTION, IVS_MARKER

HERE = Path(__file__).resolve().parent
REGISTRY = HERE / "registry.tsv"

REQUIRED_NEW_COLS = ["evidence_strength", "label_rationale", "junction_id", "junction_mapping_grade",
                     "assay_context"]


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


def main() -> int:
    df = pd.read_csv(REGISTRY, sep="\t", dtype=str).fillna("")
    v = violations(df)
    if v:
        print(f"FAIL: {len(v)} violation(s):", file=sys.stderr)
        for line in v:
            print("  -", line, file=sys.stderr)
        return 1
    print(f"PASS: {len(df)} rows valid against the labeling scheme.")

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
