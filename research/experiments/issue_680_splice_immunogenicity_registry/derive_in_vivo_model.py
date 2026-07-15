#!/usr/bin/env python3
"""Derive the in_vivo_model column ([Issue #1120](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1120)).

## Why this column exists, and why it is NOT an assay_context value

`assay_context` is a **T-cell-source** axis: every one of its values names *where the
responding T cells came from* (patient blood, patient TIL, healthy-donor IVS, an
engineered clone). An animal model is not a T-cell source - it is a **setting**, i.e.
where the experiment was run. The two are independent, and one single-valued column
cannot carry both.

#1120 was originally filed asking for an "animal model" value inside `assay_context`.
That is a category error, and grounding it against the registry is what showed why:

- **`FLLDGSANV`** (COL6A3, Kim/Immatics 2022) - in-vivo tumor control in **NSG**
  (NOD/SCID/gamma-chain-/-) mice with **adoptively transferred human T cells**
  (verified first-hand, PMC10130759).
- **`IFSESETRAKF`** (RCAN1-4, Xiong 2025) - **orthotopic xenograft** in
  **immunodeficient NSG** mice with **TCR-T cell transfer** (verified first-hand in
  the Zotero PDF: *"Immunodeficient NSG mice were intracranially injected..."*,
  *"...xenograft murine model after TCR-T cell transfer"*).

**An NSG mouse cannot produce T cells at all.** In both rows the responding T cells
are human; the mouse supplies the tumor bed and nothing immunological. Recording
either row's *T-cell source* as "animal" would therefore be flatly false - and would
overwrite the `cloned_tcr` fact that is the whole point of the column.

There IS a genuine animal-T-cell-source case - a **syngeneic, immunocompetent** host
responding with its own T cells (the Burbage-2023 mouse exon-TE shape, #699). That
case is `assay_context = animal_syngeneic`, and we hold no such row yet.

## The rule

Keyed on the in-vivo marker in `readout` first, so a row can only carry a model type
if it actually reports an in-vivo readout; the model *type* is then source-keyed from
the first-hand rationale in PROVENANCE.md. A source with an in-vivo readout but no
determinable model type derives to `unspecified` - the same no-guess convention as
`assay_context`, never an inference from a keyword.

This script OWNS the in_vivo_model column: re-running overwrites every value from
source. Pin a new row by extending the source-keyed logic below - never by editing
registry.tsv directly, or the next run silently reverts it.

Run: research/.venv/bin/python derive_in_vivo_model.py
"""
import sys
from pathlib import Path

import pandas as pd

from labeling_constants import IN_VIVO_MARKERS, IN_VIVO_MODELS, IN_VIVO_NONE

HERE = Path(__file__).resolve().parent
REG = HERE / "registry.tsv"

# Source-keyed model type, for sources that report an in-vivo animal readout.
# Both entries verified first-hand (see the module docstring + PROVENANCE.md);
# neither is inferred from the readout string.
MODEL_BY_SOURCE_SUBSTR = {
    "xiong":  "xenograft",   # orthotopic GBM xenograft, immunodeficient NSG, TCR-T transfer
    "col6a3": "xenograft",   # NSG + adoptively transferred human T cells (Kim/Immatics)
}


def has_in_vivo_readout(readout: str) -> bool:
    """Does this row report an in-vivo animal readout at all?

    Gate the whole column on this. Keying the model type on `source` alone would
    wrongly stamp a model onto every row of an in-vivo source - including Xiong's
    `VFVDGLCRAKF`, a constitutive non-splice control that was never put in a mouse.
    """
    r = str(readout).lower()
    return any(m in r for m in IN_VIVO_MARKERS)


def in_vivo_model(r) -> str:
    if not has_in_vivo_readout(r["readout"]):
        return IN_VIVO_NONE

    src = str(r["source"]).lower()
    for substr, model in MODEL_BY_SOURCE_SUBSTR.items():
        if substr in src:
            return model

    # In-vivo readout reported, but the model type is not established first-hand.
    # No guess: an in-vivo mouse experiment is xenograft-vs-syngeneic on evidence,
    # not on the plausibility of the source.
    return "unspecified"


def main() -> int:
    df = pd.read_csv(REG, sep="\t", dtype=str).fillna("")

    # Capture the ON-DISK column BEFORE deriving. This ordering is the whole point:
    # the first draft of this script compared the *just-derived* column against the
    # readout marker and called it a falsifier - but `df.apply` had already overwritten
    # any hand-edit, so the comparison was the derivation checked against itself and
    # was structurally green no matter what was in the file. A check that cannot fail
    # is not a check. (Caught in the PR #1186 bot review; the real cross-check has
    # always been the one in validate_registry.py, which reads the file fresh.)
    on_disk = df["in_vivo_model"] if "in_vivo_model" in df.columns else None

    df["in_vivo_model"] = df.apply(in_vivo_model, axis=1)

    bad = set(df["in_vivo_model"]) - IN_VIVO_MODELS
    if bad:
        print(f"FAIL: derived out-of-vocabulary in_vivo_model {bad}", file=sys.stderr)
        return 1

    # Real falsifier: a hand-edit is a row whose on-disk value disagrees with what the
    # derivation produces. This script OWNS the column, so such an edit is about to be
    # silently reverted - say so loudly instead of reverting it in silence.
    reverted = []
    if on_disk is not None:
        drift = on_disk != df["in_vivo_model"]
        for _, r in df[drift].iterrows():
            reverted.append(f"  {r['peptide'] or r['junction_id']}: "
                            f"{on_disk[r.name]!r} (on disk) -> {r['in_vivo_model']!r} (derived)")

    df.to_csv(REG, sep="\t", index=False)
    print(df["in_vivo_model"].value_counts().to_string())

    if reverted:
        print(f"\nWARNING: {len(reverted)} hand-edited in_vivo_model value(s) REVERTED "
              f"by the derivation (this script owns the column):", file=sys.stderr)
        for line in reverted:
            print(line, file=sys.stderr)

    print("\nin-vivo rows:")
    for _, r in df[df["in_vivo_model"] != IN_VIVO_NONE].iterrows():
        print(f"  {r['peptide'] or r['junction_id']:<14} {r['source']:<32} "
              f"assay_context={r['assay_context']:<12} in_vivo_model={r['in_vivo_model']}")
    return 1 if reverted else 0


if __name__ == "__main__":
    sys.exit(main())
