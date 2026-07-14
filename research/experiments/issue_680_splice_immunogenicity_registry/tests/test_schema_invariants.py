"""Two-resolution registry schema (#1086): nullable identity + typed peptide_status.

The registry used to be peptide-keyed: every row needed an amino-acid sequence to
exist, which locked out every coordinate-first source. These tests pin the schema
that replaces that assumption - identity is the coalesced junction-or-peptide key,
and an absent peptide carries an explicit reason code rather than a bare blank.
"""
import pandas as pd
import pytest

from validate_registry import violations


def _fails_with(df, needle):
    """Assert at least one violation mentions `needle`, and return the matches."""
    hits = [v for v in violations(df) if needle in v]
    assert hits, f"expected a violation mentioning {needle!r}, got: {violations(df)}"
    return hits


def test_template_row_is_green(row, frame):
    assert violations(frame(row)) == []


# --- AC1: both identity columns nullable, under at-least-one-non-null -------------


def test_peptide_null_with_junction_id_is_valid(row, frame):
    """A coordinate-first row (Zhao 2025 shape): no sequence, but a published junction."""
    row.update(peptide="", length="", tier="functional-nonscorable",
               peptide_status="published-pending")
    assert violations(frame(row)) == []


def test_junction_id_null_with_peptide_is_valid(row, frame):
    """The existing sequence-first shape: 62 of 97 rows look like this."""
    row.update(junction_id="", junction_mapping_grade="gene-mechanism")
    assert violations(frame(row)) == []


def test_both_identity_columns_null_is_rejected(row, frame):
    row.update(peptide="", length="", junction_id="",
               junction_mapping_grade="none", tier="functional-nonscorable",
               peptide_status="unpublished-idonly")
    _fails_with(frame(row), "at-least-one-non-null")


def test_violation_message_survives_a_null_peptide(row, frame):
    """The row label used to interpolate r['peptide']; a peptide-null row must still
    be nameable in an error message (by its coalesced identity)."""
    row.update(peptide="", length="", tier="functional-nonscorable",
               peptide_status="published-pending", evidence_strength="bogus")
    hits = _fails_with(frame(row), "bad evidence_strength")
    assert "chr5:33954504-33963931(-)" in hits[0]


# --- AC2: the grade <-> junction_id block is replaced, not merely relaxed ----------


@pytest.mark.parametrize("grade", ["gene-mechanism", "none"])
def test_recovered_junction_id_on_a_low_grade_row_is_green(row, frame, grade):
    """Old rule hard-blocked this, which would bar 64% of rows from ever gaining a
    junction. A junction recovered from an authoritative source is now allowed on
    any grade; the grade still records what the *source* published."""
    row.update(junction_mapping_grade=grade, notes="reason recorded")
    assert violations(frame(row)) == []


@pytest.mark.parametrize("grade", ["coords", "event-id"])
def test_high_grade_row_still_requires_a_junction_id(row, frame, grade):
    row.update(junction_mapping_grade=grade, junction_id="")
    _fails_with(frame(row), "empty junction_id")


# --- AC3: peptide_status is a typed null, harmonized with tier + grade -------------


def test_missing_peptide_status_column_is_a_schema_violation(row, frame):
    df = frame(row).drop(columns=["peptide_status"])
    _fails_with(df, "missing required column: peptide_status")


def test_bad_peptide_status_value_is_rejected(row, frame):
    row["peptide_status"] = "published"  # not in the 4-value vocab
    _fails_with(frame(row), "bad peptide_status")


def test_present_peptide_must_be_published_recovered(row, frame):
    row["peptide_status"] = "published-pending"
    _fails_with(frame(row), "peptide present but peptide_status")


@pytest.mark.parametrize("status", ["published-pending", "unpublished-idonly",
                                    "na-junction-level"])
def test_null_peptide_takes_a_null_status(row, frame, status):
    row.update(peptide="", length="", tier="functional-nonscorable",
               peptide_status=status)
    assert violations(frame(row)) == []


def test_null_peptide_may_not_be_published_recovered(row, frame):
    row.update(peptide="", length="", tier="functional-nonscorable")
    _fails_with(frame(row), "peptide absent but peptide_status")


def test_null_peptide_nulls_its_length(row, frame):
    row.update(peptide="", tier="functional-nonscorable",
               peptide_status="published-pending")  # length left at "9"
    _fails_with(frame(row), "non-empty length")


def test_present_peptide_length_must_match_the_sequence(row, frame):
    row["length"] = "10"
    _fails_with(frame(row), "length")


def test_null_peptide_may_not_be_scorable(row, frame):
    """The benchmark keys on sequence: a peptide-null row cannot enter the scored set."""
    row.update(peptide="", length="", peptide_status="published-pending")
    _fails_with(frame(row), "functional-scorable")


def test_null_peptide_needs_real_junction_evidence(row, frame):
    """With no sequence, the junction *is* the row; a gene-name-only grade cannot
    carry it even though AC2 permits a recovered junction_id on such a grade."""
    row.update(peptide="", length="", tier="functional-nonscorable",
               peptide_status="unpublished-idonly",
               junction_mapping_grade="gene-mechanism", notes="reason recorded")
    _fails_with(frame(row), "peptide-null row needs a coords/event-id grade")


# --- AC6: dedup surfaces through the validator, not only via duplicate_keys() ------


def test_duplicate_row_surfaces_as_a_validator_violation(row, frame):
    """duplicate_keys() is unit-tested directly, but the live-registry run only ever
    exercises its no-duplicates branch; pin the wiring into violations()."""
    _fails_with(frame(row, dict(row)), "duplicate identity key")


# --- #1106: `source` must map 1:1 onto studies ------------------------------------
#
# Matched pair: the same two rows, differing only in whether the second row's
# `source` is a second spelling of the first study. Same fixture, one variable
# flipped, opposite expected outcomes - so the guard is shown to fire AND to stay
# quiet. A test that only ever asserts the failure branch cannot distinguish a
# working guard from one that flags everything.


def test_two_spellings_of_one_study_are_rejected(row, frame):
    """The defect this guard exists for: 'IRIS' and 'IRIS (Pan/Xing, PNAS)' both
    match the 'iris' substring key, so every derivation agrees and nothing fires -
    while a study-level group-by silently counts them as two studies."""
    a = dict(row, source="IRIS (Pan/Xing, PNAS)", peptide="AAAAAAAAA")
    b = dict(row, source="IRIS", peptide="CCCCCCCCC")
    hits = _fails_with(frame(a, b), "source key collision")
    assert "'iris'" in hits[0]
    assert "IRIS (Pan/Xing, PNAS)" in hits[0] and "'IRIS'" in hits[0]


def test_one_spelling_per_study_is_green(row, frame):
    """The control. Two rows from two DIFFERENT studies, one spelling each, must not
    trip the guard - otherwise it would flag every multi-study registry, which is
    every registry we have."""
    a = dict(row, source="IRIS (Pan/Xing, PNAS)", peptide="AAAAAAAAA")
    b = dict(row, source="SNAF (Li 2024, Sci Transl Med)", peptide="CCCCCCCCC")
    assert [v for v in violations(frame(a, b)) if "source key collision" in v] == []


def test_collision_is_reported_once_per_study_not_once_per_row(row, frame):
    """Three rows, two spellings, one study -> one violation naming both spellings.
    Per-row reporting would bury the signal under duplicates on a large fold."""
    rows = [
        dict(row, source="Xiong 2025 (GBM)", peptide="AAAAAAAAA"),
        dict(row, source="Xiong 2025", peptide="CCCCCCCCC"),
        dict(row, source="Xiong 2025", peptide="DDDDDDDDD"),
    ]
    hits = _fails_with(frame(*rows), "source key collision")
    assert len(hits) == 1, f"expected one collision violation, got {len(hits)}"


# --- #1120: assay_context is a T-CELL-SOURCE axis; in_vivo_model is the SETTING ----
#
# The Issue originally asked for an "animal model" value inside assay_context. That is
# a category error: an animal is a venue, not a T-cell source. Both of our in-vivo rows
# are human T cells inside an immunodeficient NSG mouse - which has no T cells of its
# own - so an "animal" T-cell source would be false for both, and would overwrite the
# `cloned_tcr` fact on the COL6A3 row. These tests pin the split so it cannot be
# quietly re-merged.


def test_in_vivo_readout_requires_a_model(row, frame):
    """A readout naming an in-vivo animal experiment cannot sit at `none`. This is what
    a hand-edit into registry.tsv looks like, and the derivation cannot produce it."""
    row.update(readout="engineered-TCR IFN-g + cytotoxicity + in vivo",
               assay_context="cloned_tcr", in_vivo_model="none")
    _fails_with(frame(row), "in-vivo animal readout but")


def test_model_without_an_in_vivo_readout_is_rejected(row, frame):
    """The other direction: a model stamped onto a row that reports no in-vivo readout."""
    row.update(in_vivo_model="xenograft")  # template readout has no in-vivo term
    _fails_with(frame(row), "no in-vivo readout")


def test_bad_in_vivo_model_value_is_rejected(row, frame):
    row.update(readout="IFN-g + in vivo", in_vivo_model="mouse")
    _fails_with(frame(row), "bad in_vivo_model")


def test_xenograft_may_not_claim_the_animal_as_the_t_cell_source(row, frame):
    """**The test that would fail if the two axes were re-merged.**

    A xenograft host is immunodeficient - it has no T cells to be the source of - so
    `animal_syngeneic` (the one context in which the animal IS the T-cell source) is
    incoherent with it. This is precisely the row #1120 would originally have created
    for COL6A3, and it must not validate."""
    row.update(readout="engineered-TCR IFN-g + cytotoxicity + in vivo",
               assay_context="animal_syngeneic", in_vivo_model="xenograft")
    _fails_with(frame(row), "requires in_vivo_model 'syngeneic'")


def test_syngeneic_host_may_carry_an_engineered_tcr_source(row, frame):
    """**The combination that proves the two axes are genuinely free to vary.**

    An immunocompetent (syngeneic) host receiving adoptively-transferred cloned MOUSE
    TCR-T cells: the *setting* is syngeneic, the *T-cell source* is an engineered clone.
    Real, common design - and it must validate.

    This is why `in_vivo_model == syngeneic` does NOT imply
    `assay_context == animal_syngeneic`. A PR #1186 review finding proposed adding that
    converse coupling on the grounds that the two are definitionally equivalent; they
    were only equivalent because the *setting* vocabulary had wrongly been written as a
    claim about T-cell origin. Fixing the definition to describe the host alone (which
    is what a setting axis may say) makes this row representable, and makes the converse
    coupling wrong. Enforcing it would forbid a legitimate experiment.
    """
    row.update(readout="cloned-TCR IFN-g + cytotoxicity + in vivo",
               assay_context="cloned_tcr", in_vivo_model="syngeneic")
    assert violations(frame(row)) == []


def test_animal_syngeneic_is_valid_in_a_syngeneic_host(row, frame):
    """The control: the genuine animal-T-cell-source case (a Burbage-2023-shaped mouse
    fold) must validate, or the new context value would be unusable by construction."""
    row.update(readout="IFN-g + cytotoxicity + in vivo",
               assay_context="animal_syngeneic", in_vivo_model="syngeneic")
    assert violations(frame(row)) == []


def test_cloned_tcr_survives_an_in_vivo_readout(row, frame):
    """The regression #1120 would originally have caused: an in-vivo confirmation must
    NOT overwrite the T-cell source. A row can be `cloned_tcr` AND in-vivo-confirmed -
    both of our real in-vivo rows are exactly that."""
    row.update(readout="engineered-TCR IFN-g + cytotoxicity + in vivo",
               assay_context="cloned_tcr", in_vivo_model="xenograft")
    assert violations(frame(row)) == []


def test_live_registry_keeps_col6a3_as_cloned_tcr():
    """Pinned on the real artifact, because this is the exact regression the Issue's
    original ACs would have shipped."""
    from validate_registry import REGISTRY

    df = pd.read_csv(REGISTRY, sep="\t", dtype=str).fillna("")
    col6a3 = df[df["peptide"] == "FLLDGSANV"]
    assert len(col6a3) == 1
    assert col6a3.iloc[0]["assay_context"] == "cloned_tcr"
    assert col6a3.iloc[0]["in_vivo_model"] == "xenograft"


def test_live_registry_in_vivo_rows_are_exactly_the_marked_ones():
    """Both directions on the real registry: exactly the rows whose readout names an
    in-vivo experiment carry a model, and no others."""
    from labeling_constants import IN_VIVO_MARKERS, IN_VIVO_NONE
    from validate_registry import REGISTRY

    df = pd.read_csv(REGISTRY, sep="\t", dtype=str).fillna("")
    marked = df["readout"].str.lower().apply(lambda r: any(m in r for m in IN_VIVO_MARKERS))
    stamped = df["in_vivo_model"] != IN_VIVO_NONE
    assert (marked == stamped).all()
    assert int(stamped.sum()) == 2, "expected exactly 2 in-vivo rows (RCAN1-4, COL6A3)"


# --- AC5: no regression on the live registry --------------------------------------


def test_live_registry_validates_green():
    from validate_registry import REGISTRY

    df = pd.read_csv(REGISTRY, sep="\t", dtype=str).fillna("")
    assert len(df) == 97, "row count changed; update this guard deliberately"
    assert violations(df) == []


def test_venue_map_has_one_key_per_study():
    """The venue map is keyed one-entry-per-study, and the registry now holds one
    `source` string per study, so the two counts must agree. If they drift, either a
    study was folded without a venue key (which the venue_type check catches) or a
    study was re-split (which the collision guard catches) - this pins the pair."""
    from labeling_constants import VENUE_BY_SOURCE_SUBSTR
    from validate_registry import REGISTRY

    df = pd.read_csv(REGISTRY, sep="\t", dtype=str).fillna("")
    assert len(VENUE_BY_SOURCE_SUBSTR) == df["source"].nunique() == 12


def test_live_registry_has_one_source_string_per_study():
    """#1106 in its post-condition form: the registry carried 15 strings for 12
    studies, and no existing check could see it. This pins the invariant on the real
    artifact, so a future fold that re-splits a study fails here rather than silently
    inflating a study-level count."""
    from validate_registry import REGISTRY

    df = pd.read_csv(REGISTRY, sep="\t", dtype=str).fillna("")
    assert df["source"].nunique() == 12, (
        "distinct `source` strings != 12 studies; a fold either added a study "
        "(update this number deliberately) or re-split an existing one (fix the spelling)"
    )
