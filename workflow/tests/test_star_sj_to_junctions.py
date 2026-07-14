"""Tests for star_sj_to_junctions.py — STAR SJ.out.tab → raw_junctions.tsv conversion.

Regression test for [Issue #374](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/374):
STAR emits strand=0 (col 4 = 0) when it cannot infer strand from intron motif —
the previous inline awk silently wrote these as strand `.`, leading to wrong-
orientation flanking sequence in downstream `assemble_contigs.py`. The rescue
uses col 5 (intron motif) to recover strand for the 6 canonical/semi-canonical
motifs and drops the non-canonical (motif=0) remainder rather than corrupting it.

STAR SJ.out.tab columns (per the STAR manual):
    1: chrom
    2: intron first base (1-based)
    3: intron last base  (1-based, inclusive)
    4: strand    (0=undefined, 1=+, 2=-)
    5: motif     (0=non-canonical, 1=GT/AG +, 2=CT/AC -, 3=GC/AG +,
                  4=CT/GC -, 5=AT/AC +, 6=GT/AT -)
    6: annotated (0=novel, 1=annotated)  — used by Issue #375
    7: # uniquely-mapping reads crossing
    8: # multi-mapping reads crossing
    9: max overhang
"""

from pathlib import Path

import pytest

from star_sj_to_junctions import convert_sj_to_junctions


def _sj_line(
    chrom="chr22",
    start=101,
    end=200,
    strand=1,
    motif=1,
    annotated=0,
    unique=10,
    multi=0,
    overhang=50,
):
    """Build one SJ.out.tab line (9 tab-separated fields, integers)."""
    return "\t".join(
        str(x) for x in (chrom, start, end, strand, motif, annotated, unique, multi, overhang)
    )


class TestDirectStrand:
    """When STAR sets col 4 directly (1 or 2), no motif lookup is needed."""

    def test_strand_plus_direct(self, tmp_path):
        sj = tmp_path / "SJ.out.tab"
        sj.write_text(_sj_line(strand=1, motif=1) + "\n")
        out = tmp_path / "raw_junctions.tsv"

        convert_sj_to_junctions(sj, out)
        assert out.read_text().strip() == "chr22:101:200:+\t10\t0"

    def test_strand_minus_direct(self, tmp_path):
        sj = tmp_path / "SJ.out.tab"
        sj.write_text(_sj_line(strand=2, motif=2) + "\n")
        out = tmp_path / "raw_junctions.tsv"

        convert_sj_to_junctions(sj, out)
        assert out.read_text().strip() == "chr22:101:200:-\t10\t0"

    @pytest.mark.parametrize("strand,motif,expected", [
        (1, 2, "+"),  # col 4 says +, motif says -; col 4 wins
        (2, 1, "-"),  # col 4 says -, motif says +; col 4 wins
        (1, 0, "+"),  # col 4 says +, motif is non-canonical; col 4 wins
        (2, 0, "-"),  # col 4 says -, motif is non-canonical; col 4 wins
    ])
    def test_direct_strand_takes_priority_over_motif(
        self, tmp_path, strand, motif, expected
    ):
        """Core invariant: when STAR sets col 4 (1 or 2), motif is ignored."""
        sj = tmp_path / "SJ.out.tab"
        sj.write_text(_sj_line(strand=strand, motif=motif) + "\n")
        out = tmp_path / "raw_junctions.tsv"

        convert_sj_to_junctions(sj, out)
        assert out.read_text().strip() == f"chr22:101:200:{expected}\t10\t0"


class TestMotifRescue:
    """When col 4 = 0, col 5 (intron motif) recovers strand for codes 1-6."""

    @pytest.mark.parametrize("motif,expected_strand", [
        (1, "+"),  # GT/AG
        (3, "+"),  # GC/AG
        (5, "+"),  # AT/AC
    ])
    def test_plus_motifs_rescue_to_plus(self, tmp_path, motif, expected_strand):
        sj = tmp_path / "SJ.out.tab"
        sj.write_text(_sj_line(strand=0, motif=motif) + "\n")
        out = tmp_path / "raw_junctions.tsv"

        convert_sj_to_junctions(sj, out)
        assert out.read_text().strip() == f"chr22:101:200:{expected_strand}\t10\t0"

    @pytest.mark.parametrize("motif,expected_strand", [
        (2, "-"),  # CT/AC
        (4, "-"),  # CT/GC
        (6, "-"),  # GT/AT
    ])
    def test_minus_motifs_rescue_to_minus(self, tmp_path, motif, expected_strand):
        sj = tmp_path / "SJ.out.tab"
        sj.write_text(_sj_line(strand=0, motif=motif) + "\n")
        out = tmp_path / "raw_junctions.tsv"

        convert_sj_to_junctions(sj, out)
        assert out.read_text().strip() == f"chr22:101:200:{expected_strand}\t10\t0"

    def test_motif_zero_strand_zero_dropped(self, tmp_path):
        """strand=0 + motif=0 (truly non-canonical) is dropped, not emitted as '.'.

        The whole point of the bug fix — silently emitting these as strand `.`
        causes wrong-orientation flanking sequence downstream.
        """
        sj = tmp_path / "SJ.out.tab"
        sj.write_text(_sj_line(strand=0, motif=0) + "\n")
        out = tmp_path / "raw_junctions.tsv"

        convert_sj_to_junctions(sj, out)
        assert out.read_text() == ""

    def test_motif_out_of_range_dropped(self, tmp_path):
        """Defensive: unknown motif code (>6) is treated as motif=0 (drop)."""
        sj = tmp_path / "SJ.out.tab"
        sj.write_text(_sj_line(strand=0, motif=7) + "\n")
        out = tmp_path / "raw_junctions.tsv"

        convert_sj_to_junctions(sj, out)
        assert out.read_text() == ""


class TestFiltering:
    """Reads-count filter + malformed-line tolerance."""

    def test_zero_unique_reads_skipped(self, tmp_path):
        sj = tmp_path / "SJ.out.tab"
        sj.write_text(_sj_line(unique=0) + "\n")
        out = tmp_path / "raw_junctions.tsv"

        convert_sj_to_junctions(sj, out)
        assert out.read_text() == ""

    def test_skips_short_malformed_lines(self, tmp_path):
        sj = tmp_path / "SJ.out.tab"
        sj.write_text(
            "chr22\t101\t200\t1\n"  # only 4 fields — malformed
            + _sj_line() + "\n"
        )
        out = tmp_path / "raw_junctions.tsv"

        convert_sj_to_junctions(sj, out)
        assert out.read_text().strip() == "chr22:101:200:+\t10\t0"

    def test_skips_blank_lines(self, tmp_path):
        sj = tmp_path / "SJ.out.tab"
        sj.write_text(
            "\n"
            + _sj_line() + "\n"
            "   \n"
        )
        out = tmp_path / "raw_junctions.tsv"

        convert_sj_to_junctions(sj, out)
        assert out.read_text().strip() == "chr22:101:200:+\t10\t0"

    def test_skips_non_integer_count(self, tmp_path):
        sj = tmp_path / "SJ.out.tab"
        sj.write_text(
            "chr22\t101\t200\t1\t1\t0\tnotanum\t0\t50\n"
            + _sj_line() + "\n"
        )
        out = tmp_path / "raw_junctions.tsv"

        convert_sj_to_junctions(sj, out)
        assert out.read_text().strip() == "chr22:101:200:+\t10\t0"


class TestMultipleRecords:
    """Mixed direct + rescued + dropped in one file."""

    def test_mixed_records_preserve_order(self, tmp_path):
        sj = tmp_path / "SJ.out.tab"
        sj.write_text(
            _sj_line(strand=1, motif=1, unique=5) + "\n"            # direct +
            + _sj_line(start=300, end=400, strand=0, motif=2, unique=3) + "\n"  # rescued -
            + _sj_line(start=500, end=600, strand=0, motif=0, unique=7) + "\n"  # dropped
            + _sj_line(start=700, end=800, strand=2, motif=2, unique=2) + "\n"  # direct -
        )
        out = tmp_path / "raw_junctions.tsv"

        convert_sj_to_junctions(sj, out)
        lines = out.read_text().strip().splitlines()
        assert lines == [
            "chr22:101:200:+\t5\t0",
            "chr22:300:400:-\t3\t0",
            "chr22:700:800:-\t2\t0",
        ]


class TestAnnotatedFlagColumn:
    """Issue #375 — col 6 (annotated flag) is surfaced as a 3rd output column."""

    def test_annotated_row_emits_flag_one(self, tmp_path):
        sj = tmp_path / "SJ.out.tab"
        sj.write_text(_sj_line(strand=1, motif=1, annotated=1) + "\n")
        out = tmp_path / "raw_junctions.tsv"

        convert_sj_to_junctions(sj, out)
        assert out.read_text().strip() == "chr22:101:200:+\t10\t1"

    def test_novel_row_emits_flag_zero(self, tmp_path):
        sj = tmp_path / "SJ.out.tab"
        sj.write_text(_sj_line(strand=1, motif=1, annotated=0) + "\n")
        out = tmp_path / "raw_junctions.tsv"

        convert_sj_to_junctions(sj, out)
        assert out.read_text().strip() == "chr22:101:200:+\t10\t0"

    def test_rescued_strand_row_carries_annotated_flag(self, tmp_path):
        # strand=0 rescued via motif, annotated=1 → still emits the flag.
        sj = tmp_path / "SJ.out.tab"
        sj.write_text(_sj_line(strand=0, motif=2, annotated=1) + "\n")
        out = tmp_path / "raw_junctions.tsv"

        convert_sj_to_junctions(sj, out)
        assert out.read_text().strip() == "chr22:101:200:-\t10\t1"

    def test_mixed_annotated_and_novel_rows(self, tmp_path):
        sj = tmp_path / "SJ.out.tab"
        sj.write_text(
            _sj_line(strand=1, motif=1, annotated=1, unique=5) + "\n"
            + _sj_line(start=300, end=400, strand=2, motif=2, annotated=0, unique=3) + "\n"
        )
        out = tmp_path / "raw_junctions.tsv"

        convert_sj_to_junctions(sj, out)
        assert out.read_text().strip().splitlines() == [
            "chr22:101:200:+\t5\t1",
            "chr22:300:400:-\t3\t0",
        ]


class TestReturnValue:
    """The function returns the number of junctions written."""

    def test_returns_count_of_written_records(self, tmp_path):
        sj = tmp_path / "SJ.out.tab"
        sj.write_text(
            _sj_line(strand=1, motif=1, unique=5) + "\n"
            + _sj_line(start=300, end=400, strand=0, motif=2, unique=3) + "\n"
            + _sj_line(start=500, end=600, strand=0, motif=0, unique=7) + "\n"  # dropped
        )
        out = tmp_path / "raw_junctions.tsv"

        assert convert_sj_to_junctions(sj, out) == 2


class TestUniquenessKnob:
    """The multimapper (col 8) semantic is an explicit policy, not a column accident.

    [Issue #1118](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1118).

    Historically this script read `SJ.out.tab` col 7 (uniquely-mapping reads)
    only and never col 8 (multi-mapping), so a junction supported *only* by
    multimappers was discarded outright - while the HISAT2 path counted every
    read. Same FASTQ, different candidate set, decided by nothing but which
    column a script happened to read.

    [Issue #1122](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1122)
    ruled the unique-only semantic **structurally unsound for a matched
    tumor/normal design** (an independent per-read gate on each arm can destroy
    normal support while leaving tumor untouched, manufacturing a false
    `tumor_exclusive` candidate). So the *default* here is count-all - col 7 +
    col 8 - and unique-only survives only as an explicit opt-in knob
    (`alignment.uniqueness_filter.enabled`), matching the HISAT2 lever from
    [Issue #919](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/919).
    """

    # --- default: count all reads (the #1122 semantic) ---

    def test_default_sums_unique_and_multi_reads(self, tmp_path):
        sj = tmp_path / "SJ.out.tab"
        sj.write_text(_sj_line(unique=3, multi=2) + "\n")
        out = tmp_path / "raw_junctions.tsv"

        convert_sj_to_junctions(sj, out)

        # 3 unique + 2 multi = 5. Reading col 7 alone would say 3.
        assert out.read_text().strip() == "chr22:101:200:+\t5\t0"

    def test_default_keeps_a_multimapper_only_junction(self, tmp_path):
        """The bug in one test: col-7-only silently discarded this junction.

        This is the shape that produced the false `tumor_exclusive` in #1122 -
        support that exists in one arm and is annihilated in the other.
        """
        sj = tmp_path / "SJ.out.tab"
        sj.write_text(_sj_line(unique=0, multi=4) + "\n")
        out = tmp_path / "raw_junctions.tsv"

        assert convert_sj_to_junctions(sj, out) == 1
        assert out.read_text().strip() == "chr22:101:200:+\t4\t0"

    def test_default_still_drops_a_genuinely_zero_read_junction(self, tmp_path):
        """Count-all must not resurrect junctions with no support at all."""
        sj = tmp_path / "SJ.out.tab"
        sj.write_text(_sj_line(unique=0, multi=0) + "\n")
        out = tmp_path / "raw_junctions.tsv"

        assert convert_sj_to_junctions(sj, out) == 0
        assert out.read_text() == ""

    # --- opt-in: unique-only (the pre-#1118 behavior, now explicit) ---

    def test_unique_only_counts_col7_and_ignores_col8(self, tmp_path):
        sj = tmp_path / "SJ.out.tab"
        sj.write_text(_sj_line(unique=3, multi=2) + "\n")
        out = tmp_path / "raw_junctions.tsv"

        convert_sj_to_junctions(sj, out, unique_only=True)

        assert out.read_text().strip() == "chr22:101:200:+\t3\t0"

    def test_unique_only_drops_a_multimapper_only_junction(self, tmp_path):
        """Preserves the historical behavior exactly - but only when asked for."""
        sj = tmp_path / "SJ.out.tab"
        sj.write_text(_sj_line(unique=0, multi=4) + "\n")
        out = tmp_path / "raw_junctions.tsv"

        assert convert_sj_to_junctions(sj, out, unique_only=True) == 0
        assert out.read_text() == ""

    # --- the matched-pair control: the SAME input, one flag flipped ---

    def test_knob_is_the_only_difference_between_the_two_semantics(self, tmp_path):
        """Matched-pair control: identical input, opposite expected outcomes.

        If this passes with both arms equal, the knob is not wired to anything.
        """
        sj = tmp_path / "SJ.out.tab"
        sj.write_text(
            _sj_line(start=101, end=200, unique=5, multi=1) + "\n"
            + _sj_line(start=301, end=400, unique=0, multi=6) + "\n"  # multimapper-only
        )
        count_all = tmp_path / "all.tsv"
        unique_only = tmp_path / "uniq.tsv"

        n_all = convert_sj_to_junctions(sj, count_all, unique_only=False)
        n_uniq = convert_sj_to_junctions(sj, unique_only, unique_only=True)

        assert n_all == 2 and n_uniq == 1
        assert count_all.read_text().strip().splitlines() == [
            "chr22:101:200:+\t6\t0",
            "chr22:301:400:+\t6\t0",
        ]
        assert unique_only.read_text().strip().splitlines() == ["chr22:101:200:+\t5\t0"]

    # --- robustness: a short row must not crash the count-all path ---

    def test_row_missing_col8_is_tolerated_as_zero_multimappers(self, tmp_path):
        """STAR always emits 9 columns, but a truncated row must not explode.

        The pre-existing malformed guard only required >= 7 fields, so a 7-field
        row reached the parser. Count-all must treat the absent col 8 as 0
        rather than raise.
        """
        sj = tmp_path / "SJ.out.tab"
        sj.write_text("chr22\t101\t200\t1\t1\t0\t5\n")  # 7 fields: no col 8, no col 9
        out = tmp_path / "raw_junctions.tsv"

        assert convert_sj_to_junctions(sj, out) == 1
        assert out.read_text().strip() == "chr22:101:200:+\t5\t0"

    # --- the audit log must be true in BOTH modes, not just the convenient one ---

    def test_multimapper_only_count_is_reported_under_both_semantics(self, tmp_path, caplog):
        """The audit counter is what makes the two semantics comparable.

        Bot review on [PR #1168](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/1168)
        caught it sitting *after* the `reads <= 0` gate: under `--unique-only` a
        multimapper-only row has `reads == 0` and is dropped by that gate, so the
        counter reported 0 in the exact mode where the number matters - a log line
        reading "0 dropped" while silently dropping them. It is now counted from
        the raw columns, before the gate.
        """
        sj = tmp_path / "SJ.out.tab"
        sj.write_text(
            _sj_line(start=101, end=200, unique=5, multi=1) + "\n"   # normal support
            + _sj_line(start=301, end=400, unique=0, multi=6) + "\n"  # multimapper-only
            + _sj_line(start=501, end=600, unique=0, multi=2) + "\n"  # multimapper-only
            + _sj_line(start=701, end=800, unique=0, multi=0) + "\n"  # no support at all
        )

        import logging

        with caplog.at_level(logging.INFO):
            convert_sj_to_junctions(sj, tmp_path / "all.tsv", unique_only=False)
        assert "Junctions supported only by multimappers: 2 (kept)" in caplog.text

        caplog.clear()
        with caplog.at_level(logging.INFO):
            convert_sj_to_junctions(sj, tmp_path / "uniq.tsv", unique_only=True)
        # The same 2 junctions exist; under unique-only they are DROPPED, not absent.
        assert "Junctions supported only by multimappers: 2 (dropped)" in caplog.text
