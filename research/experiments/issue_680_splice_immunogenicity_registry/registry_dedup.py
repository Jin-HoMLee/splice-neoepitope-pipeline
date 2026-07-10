#!/usr/bin/env python3
"""Identity, dedup, and the junction-resolution view for the two-resolution registry (#1086).

Row identity is the *coalesced* key `COALESCE(junction_id, peptide)`: a coordinate-first
source (junction, no sequence) and a sequence-first source (sequence, no junction) each
need their own key column nullable, so neither column alone can serve as the identity.

Dedup keys on the full `(junction_id, peptide, hla)` triple rather than on the coalesced
identity, because one junction legitimately yields several distinct peptides - the Kim
2025 constitutive-intron event `ci@16:719606:720123:+|16:719606:719607:+` carries three
A*02:01 peptides in the live registry. Immunogenicity is a property of the peptide-HLA
pair and is irreducible, so those are three rows, not one row to merge. Junction-level
consumers get that grouping from `junction_view` instead.
"""
from collections import Counter, defaultdict

import pandas as pd


def _s(value) -> str:
    """Normalize a registry cell to a stripped string (NaN and None read as empty)."""
    if value is None or (isinstance(value, float) and pd.isna(value)):
        return ""
    return str(value).strip()


def row_identity(row) -> str:
    """The coalesced junction-or-peptide identity of a row.

    Raises ValueError when both identity columns are empty - the at-least-one-non-null
    invariant, surfaced here so callers cannot silently key a row on a blank.
    """
    junction, peptide = _s(row.get("junction_id")), _s(row.get("peptide"))
    if not junction and not peptide:
        raise ValueError(
            "at-least-one-non-null violated: a row needs a junction_id, a peptide, or both"
        )
    return junction or peptide


def dedup_key(row) -> tuple[str, str, str]:
    """The full identity triple. Either identity column may be empty."""
    return _s(row.get("junction_id")), _s(row.get("peptide")), _s(row.get("hla"))


def duplicate_keys(df: pd.DataFrame) -> list[tuple[str, str, str]]:
    """Every (junction_id, peptide, hla) triple occurring more than once, sorted."""
    counts = Counter(dedup_key(r) for r in df.to_dict("records"))
    return sorted(k for k, n in counts.items() if n > 1)


def junction_view(df: pd.DataFrame) -> dict[str, list[str]]:
    """Group the registry at junction resolution: junction_id -> its sorted peptides.

    Rows with no junction_id are omitted (they have no junction-level identity). A
    junction whose rows are all peptide-null maps to an empty list, which is a real
    junction-level positive - not a missing value.
    """
    view: dict[str, list[str]] = defaultdict(list)
    for row in df.to_dict("records"):
        junction, peptide = _s(row.get("junction_id")), _s(row.get("peptide"))
        if not junction:
            continue
        peptides = view[junction]  # touch, so a wholly peptide-null junction still appears
        if peptide:
            peptides.append(peptide)
    return {j: sorted(p) for j, p in view.items()}
