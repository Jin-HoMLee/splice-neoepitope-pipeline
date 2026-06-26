#!/usr/bin/env python3
"""Validate registry.tsv against the documented labeling scheme (Issue #735).

Run: research/.venv/bin/python validate_registry.py
Exit 0 = clean; exit 1 = violations printed to stderr.
"""
import sys
import pandas as pd
from pathlib import Path

HERE = Path(__file__).resolve().parent
REGISTRY = HERE / "registry.tsv"

GRADES = {"coords", "event-id", "gene-mechanism", "none"}
STRENGTHS = {"strong", "weak", "hard", "soft", "na"}
REQUIRED_NEW_COLS = ["evidence_strength", "label_rationale", "junction_id", "junction_mapping_grade"]

# Effector vs detection-only keyword sets used to cross-check positives.
EFFECTOR = ("ifn", "elispot", "cytotox", "granzyme", "gzmb", "cd107", "cd137",
            "degranul", "killing", "ldh", "caspase", "incucyte", "tcr", "tnf",
            "stabiliz", "activation", "in vivo")
DETECTION = ("tetramer", "dextramer", "multimer")


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
        # grade 'none' must record a reason in notes
        if r["junction_mapping_grade"] == "none" and not str(r["notes"]).strip():
            out.append(f"{rid}: grade 'none' with no reason in notes")
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

    # decoy seed check
    DECOY = HERE / "decoy_negatives" / "presented_decoys_681.tsv"
    if DECOY.exists():
        dd = pd.read_csv(DECOY, sep="\t", dtype=str).fillna("")
        assert len(dd) == 13, f"expected 13 presented decoys, got {len(dd)}"
        assert (dd["decoy_tier"] == "presented-decoy").all(), "bad decoy_tier value"
        assert dd["peptide"].is_unique, "duplicate decoy peptide"
        assert dd["peptide"].str.len().between(8, 11).all(), "decoy peptide length out of range"
        print("PASS: 13 presented decoys valid.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
