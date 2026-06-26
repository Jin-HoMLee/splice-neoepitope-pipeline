#!/usr/bin/env python3
"""Derive evidence_strength + seed label_rationale (Issue #735, Task 3)."""
import pandas as pd
from pathlib import Path

HERE = Path(__file__).resolve().parent
REG = HERE / "registry.tsv"
EFFECTOR = ("ifn", "elispot", "cytotox", "granzyme", "gzmb", "cd107", "cd137",
            "degranul", "killing", "ldh", "caspase", "incucyte", "tcr", "tnf",
            "activation", "in vivo")
DETECTION = ("tetramer", "dextramer", "multimer")


def strength(r):
    ro = str(r["readout"]).lower()
    tier = r["tier"]
    if r["label"] == "untested" or tier in ("presentation-prevalence", "negative-control-not-splice"):
        return "na"
    if r["label"] == "negative":
        return "soft" if "healthy-donor ivs" in ro or "ivs" in ro else "hard"
    # positive
    if any(k in ro for k in EFFECTOR):
        return "strong"
    if any(k in ro for k in DETECTION):
        return "weak"
    return "na"  # validator will catch; investigate by hand


RATIONALE = {
    "candidate-negative": "tested-negative in healthy-donor IVS; failed-to-prime != intrinsically non-immunogenic (soft negative).",
    "hard-negative-true-splice": "splice-derived, MS-presented, functional assay performed, no T-cell response (hard negative).",
    "presentation-prevalence": "presented + prevalence/survival evidence but no functional T-cell assay performed; untested.",
    "negative-control-not-splice": "constitutive same-locus control; fails the splice gate; control only.",
    "functional-nonscorable": "functionally validated but exact sequence unpublished; cannot be keyed.",
    "candidate": "positive whose 2nd adversarial-verify pass did not confirm; held below scorable pending re-verify.",
}


def rationale(r):
    if r["tier"] in RATIONALE:
        return RATIONALE[r["tier"]]
    es = r["evidence_strength"]
    if es == "strong":
        return f"effector readout ({r['readout']}) -> strong positive."
    if es == "weak":
        return f"antigen-specific detection only ({r['readout']}) -> weak positive."
    return "REVIEW: rule did not assign a rationale."


df = pd.read_csv(REG, sep="\t", dtype=str).fillna("")
df["evidence_strength"] = df.apply(strength, axis=1)
df["label_rationale"] = df.apply(rationale, axis=1)
df.to_csv(REG, sep="\t", index=False)
print(df["evidence_strength"].value_counts().to_string())
print("rows needing manual review:",
      int((df["label_rationale"].str.startswith("REVIEW")).sum()))
