#!/usr/bin/env python3
"""Derive the assay_context column (Issue #823).

assay_context records *which immunological system* produced each row's
functional readout, so a scoring run can weight rows by assay realism without
re-parsing free-text notes. The assignment is **source-keyed**: the `source`
column is the stable curation key, and the per-source rationale is recorded
first-hand in PROVENANCE.md. The only sub-source split is Merlotti's
ex-vivo-vs-TIL distinction, read from the per-context detail in `notes`.

`unspecified` is assigned where the held provenance records the assay but not
its T-cell source (SNAF / Kim / Xiong / Fisher / POSTN / Kwok) - the
verify-against-source rule forbids guessing the patient-vs-donor system, and a
scoring run should apply no context weight to those rows.

This script OWNS the assay_context column: re-running overwrites every value
from source. So when an `unspecified` row is later pinned (the future
read-the-methods pass), make the pin by extending the source-keyed logic
below - never by editing registry.tsv directly, or the next run silently
reverts it.

Run: research/.venv/bin/python derive_assay_context.py
"""
import pandas as pd
from pathlib import Path

from labeling_constants import IVS_MARKER

HERE = Path(__file__).resolve().parent
REG = HERE / "registry.tsv"


def assay_context(r):
    src = str(r["source"]).lower()
    tier = r["tier"]
    notes = str(r["notes"]).lower()

    # tier overrides (independent of source): a prevalence/presentation-only row
    # or a constitutive non-splice control never carries a T-cell-assay context.
    if tier == "presentation-prevalence":
        return "prevalence_only"
    if tier == "negative-control-not-splice":
        return "na"

    # source-keyed assignment (rationale per source in PROVENANCE.md):
    if "manoharan" in src:
        return "healthy_donor_ivs"   # healthy-donor moDC-primed IVS (positives + soft negatives)
    if "bigot" in src:
        return "patient_exvivo"      # A2+ patient PBMC/TIL ex vivo (S2 panel)
    if "merlotti" in src:
        # ex-vivo blood detection where present; otherwise Day-20 TIL / draining-LN.
        # accept the spaced variant too so a future "ex vivo" doesn't silently fall to TIL.
        return "patient_exvivo" if ("ex-vivo" in notes or "ex vivo" in notes) else "patient_til"
    if "iris" in src:
        return "cloned_tcr"          # JPTCR engineered-TCR IFN-g / cytotoxicity readout
    if "col6a3" in src:
        return "cloned_tcr"          # Kim/Immatics 2022: affinity-enhanced + natural COL6A3 TCRs, no primary patient/donor detection

    # SNAF / Kim / Xiong / Fisher / POSTN / Kwok: assay reported, T-cell source
    # not determinable from held provenance -> no guess.
    return "unspecified"


df = pd.read_csv(REG, sep="\t", dtype=str).fillna("")
df["assay_context"] = df.apply(assay_context, axis=1)
df.to_csv(REG, sep="\t", index=False)

print(df["assay_context"].value_counts().to_string())
# sanity: healthy_donor_ivs <-> IVS marker in readout (ties to the #735 IVS rule)
ivs_rows = df["readout"].str.lower().str.contains(IVS_MARKER, regex=False)
mismatch = (ivs_rows != (df["assay_context"] == "healthy_donor_ivs")).sum()
print("IVS-marker/assay_context mismatches:", int(mismatch))
