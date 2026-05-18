"""Tests for star_sj_to_junctions.py — STAR SJ.out.tab → junctions.tsv conversion.

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
        out = tmp_path / "junctions.tsv"

        convert_sj_to_junctions(sj, out)
        assert out.read_text().strip() == "chr22:101:200:+\t10"

    def test_strand_minus_direct(self, tmp_path):
        sj = tmp_path / "SJ.out.tab"
        sj.write_text(_sj_line(strand=2, motif=2) + "\n")
        out = tmp_path / "junctions.tsv"

        convert_sj_to_junctions(sj, out)
        assert out.read_text().strip() == "chr22:101:200:-\t10"

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
        out = tmp_path / "junctions.tsv"

        convert_sj_to_junctions(sj, out)
        assert out.read_text().strip() == f"chr22:101:200:{expected}\t10"


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
        out = tmp_path / "junctions.tsv"

        convert_sj_to_junctions(sj, out)
        assert out.read_text().strip() == f"chr22:101:200:{expected_strand}\t10"

    @pytest.mark.parametrize("motif,expected_strand", [
        (2, "-"),  # CT/AC
        (4, "-"),  # CT/GC
        (6, "-"),  # GT/AT
    ])
    def test_minus_motifs_rescue_to_minus(self, tmp_path, motif, expected_strand):
        sj = tmp_path / "SJ.out.tab"
        sj.write_text(_sj_line(strand=0, motif=motif) + "\n")
        out = tmp_path / "junctions.tsv"

        convert_sj_to_junctions(sj, out)
        assert out.read_text().strip() == f"chr22:101:200:{expected_strand}\t10"

    def test_motif_zero_strand_zero_dropped(self, tmp_path):
        """strand=0 + motif=0 (truly non-canonical) is dropped, not emitted as '.'.

        The whole point of the bug fix — silently emitting these as strand `.`
        causes wrong-orientation flanking sequence downstream.
        """
        sj = tmp_path / "SJ.out.tab"
        sj.write_text(_sj_line(strand=0, motif=0) + "\n")
        out = tmp_path / "junctions.tsv"

        convert_sj_to_junctions(sj, out)
        assert out.read_text() == ""

    def test_motif_out_of_range_dropped(self, tmp_path):
        """Defensive: unknown motif code (>6) is treated as motif=0 (drop)."""
        sj = tmp_path / "SJ.out.tab"
        sj.write_text(_sj_line(strand=0, motif=7) + "\n")
        out = tmp_path / "junctions.tsv"

        convert_sj_to_junctions(sj, out)
        assert out.read_text() == ""


class TestFiltering:
    """Reads-count filter + malformed-line tolerance."""

    def test_zero_unique_reads_skipped(self, tmp_path):
        sj = tmp_path / "SJ.out.tab"
        sj.write_text(_sj_line(unique=0) + "\n")
        out = tmp_path / "junctions.tsv"

        convert_sj_to_junctions(sj, out)
        assert out.read_text() == ""

    def test_skips_short_malformed_lines(self, tmp_path):
        sj = tmp_path / "SJ.out.tab"
        sj.write_text(
            "chr22\t101\t200\t1\n"  # only 4 fields — malformed
            + _sj_line() + "\n"
        )
        out = tmp_path / "junctions.tsv"

        convert_sj_to_junctions(sj, out)
        assert out.read_text().strip() == "chr22:101:200:+\t10"

    def test_skips_blank_lines(self, tmp_path):
        sj = tmp_path / "SJ.out.tab"
        sj.write_text(
            "\n"
            + _sj_line() + "\n"
            "   \n"
        )
        out = tmp_path / "junctions.tsv"

        convert_sj_to_junctions(sj, out)
        assert out.read_text().strip() == "chr22:101:200:+\t10"

    def test_skips_non_integer_count(self, tmp_path):
        sj = tmp_path / "SJ.out.tab"
        sj.write_text(
            "chr22\t101\t200\t1\t1\t0\tnotanum\t0\t50\n"
            + _sj_line() + "\n"
        )
        out = tmp_path / "junctions.tsv"

        convert_sj_to_junctions(sj, out)
        assert out.read_text().strip() == "chr22:101:200:+\t10"


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
        out = tmp_path / "junctions.tsv"

        convert_sj_to_junctions(sj, out)
        lines = out.read_text().strip().splitlines()
        assert lines == [
            "chr22:101:200:+\t5",
            "chr22:300:400:-\t3",
            "chr22:700:800:-\t2",
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
        out = tmp_path / "junctions.tsv"

        assert convert_sj_to_junctions(sj, out) == 2
