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


# --- AC5: no regression on the live registry --------------------------------------


def test_live_registry_validates_green():
    from validate_registry import REGISTRY

    df = pd.read_csv(REGISTRY, sep="\t", dtype=str).fillna("")
    assert len(df) == 97, "row count changed; update this guard deliberately"
    assert violations(df) == []


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
