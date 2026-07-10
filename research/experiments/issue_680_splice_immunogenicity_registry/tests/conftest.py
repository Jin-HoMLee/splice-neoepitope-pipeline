import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))


# A real, currently-valid registry row (SNAF IPDSQGNDI), used as the template every
# test mutates. Keeping it verbatim means a test failure is about the invariant under
# test, not about an unrelated field drifting out of the controlled vocabularies.
_TEMPLATE = {
    "peptide": "IPDSQGNDI",
    "gene": "SLC45A2",
    "hla": "HLA-C*04",
    "hla_resolution": "2-digit",
    "length": "9",
    "splice_mechanism": "alt-5'",
    "splice_mechanism_canonical": "alt_5p_ss",
    "source": "SNAF (Li 2024)",
    "readout": "IFN-g (41% CD8+) + MS",
    "label": "positive",
    "tier": "functional-scorable",
    "confidence": "high",
    "notes": "Data S1 junction E3.2_33963931-E4.2; Fig5 C*04_IPD",
    "provenance_grade": "inferred",
    "provenance_ref": "Fig5C label + Data S1-TCGA sheet",
    "evidence_strength": "strong",
    "label_rationale": "effector readout -> strong positive.",
    "junction_id": "chr5:33954504-33963931(-)",
    "junction_mapping_grade": "coords",
    "assay_context": "unspecified",
    "venue_type": "journal",
    "peptide_status": "published-recovered",
}


@pytest.fixture
def row():
    """A single valid row as a plain dict; mutate then wrap with `frame`."""
    return dict(_TEMPLATE)


@pytest.fixture
def frame():
    """Wrap one or more row dicts into the string-typed DataFrame the validator sees."""

    def _frame(*rows):
        return pd.DataFrame(list(rows), dtype=str).fillna("")

    return _frame
