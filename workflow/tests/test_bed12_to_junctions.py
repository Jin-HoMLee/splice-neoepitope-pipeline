"""Tests for bed12_to_junctions.py — regtools BED12 → raw_junctions.tsv conversion.

Regression test for [Issue #370](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/370):
The HISAT2/regtools path previously emitted BED12 chromStart/chromEnd (anchor
outer boundaries) as junction donor/acceptor coords. The correct coords are
derived from blockSizes/blockStarts and refer to the intron interior.
"""

from pathlib import Path

import pytest

from bed12_to_junctions import convert_bed12_to_junctions


# ---------------------------------------------------------------------------
# Fixture: a single annotated junction with anchor padding
# ---------------------------------------------------------------------------
#
# Intron: chr22:[100, 200) on + strand (0-based half-open)
# Anchors: 50 bp on each side
# So regtools emits a BED12 block:
#   chromStart = 50  (donor - leftAnchorLen)
#   chromEnd   = 250 (acceptor + rightAnchorLen)
#   blockSizes  = "50,50"
#   blockStarts = "0,150"  (offset of right anchor from chromStart)
#
# Expected raw_junctions.tsv line:
#   chr22:101:200:+\t10
#         ^^^  ^^^
#       1-based  0-based
#       donor   half-open
#               intron end
# i.e. <chr>:<1-based donor>:<intron-end-exclusive>:<strand>\t<reads>
# This matches the convention consumed by filter_junctions._parse_junction_id,
# which does `start = int(parts[1]) - 1` and `end = int(parts[2])`.

_BED12_LINE = (
    "chr22\t50\t250\tJUNC0001\t10\t+\t"
    "50\t250\t255,0,0\t2\t50,50\t0,150"
)


class TestConvertBed12ToJunctions:
    def test_emits_intron_coords_not_anchor_outers(self, tmp_path):
        bed = tmp_path / "regtools.bed"
        bed.write_text(_BED12_LINE + "\n")
        out = tmp_path / "raw_junctions.tsv"

        convert_bed12_to_junctions(bed, out)

        lines = out.read_text().strip().splitlines()
        assert lines == ["chr22:101:200:+\t10"]

    def test_skips_zero_read_records(self, tmp_path):
        zero_reads = _BED12_LINE.replace("\t10\t+\t", "\t0\t+\t")
        bed = tmp_path / "regtools.bed"
        bed.write_text(zero_reads + "\n")
        out = tmp_path / "raw_junctions.tsv"

        convert_bed12_to_junctions(bed, out)
        assert out.read_text() == ""

    def test_minus_strand(self, tmp_path):
        minus = _BED12_LINE.replace("\t+\t", "\t-\t")
        bed = tmp_path / "regtools.bed"
        bed.write_text(minus + "\n")
        out = tmp_path / "raw_junctions.tsv"

        convert_bed12_to_junctions(bed, out)
        assert out.read_text().strip() == "chr22:101:200:-\t10"

    def test_multiple_records(self, tmp_path):
        second = (
            "chr1\t1000\t1300\tJUNC0002\t5\t+\t"
            "1000\t1300\t255,0,0\t2\t100,100\t0,200"
        )
        bed = tmp_path / "regtools.bed"
        bed.write_text(_BED12_LINE + "\n" + second + "\n")
        out = tmp_path / "raw_junctions.tsv"

        convert_bed12_to_junctions(bed, out)
        lines = out.read_text().strip().splitlines()
        # Second intron: donor = 1000 + 100 = 1100 (0-based) → 1101 (1-based)
        #               end   = 1000 + 200 = 1200 (0-based exclusive)
        assert lines == [
            "chr22:101:200:+\t10",
            "chr1:1101:1200:+\t5",
        ]

    def test_unstranded_record(self, tmp_path):
        unstranded = _BED12_LINE.replace("\t+\t", "\t.\t")
        bed = tmp_path / "regtools.bed"
        bed.write_text(unstranded + "\n")
        out = tmp_path / "raw_junctions.tsv"

        convert_bed12_to_junctions(bed, out)
        assert out.read_text().strip() == "chr22:101:200:.\t10"

    def test_skips_malformed_short_lines(self, tmp_path):
        # Fewer than 12 fields → silently skipped
        bed = tmp_path / "regtools.bed"
        bed.write_text(
            "chr22\t50\t250\tjunc_short\t10\t+\n"  # only 6 fields
            + _BED12_LINE + "\n"
        )
        out = tmp_path / "raw_junctions.tsv"

        convert_bed12_to_junctions(bed, out)
        # Only the well-formed line emits a junction
        assert out.read_text().strip() == "chr22:101:200:+\t10"

    def test_skips_comment_and_blank_lines(self, tmp_path):
        bed = tmp_path / "regtools.bed"
        bed.write_text(
            "# regtools v1.0.0 header comment\n"
            "\n"
            + _BED12_LINE + "\n"
            "   \n"  # whitespace-only
        )
        out = tmp_path / "raw_junctions.tsv"

        convert_bed12_to_junctions(bed, out)
        assert out.read_text().strip() == "chr22:101:200:+\t10"

    def test_annotated_junction_matches_gencode_reference(self, tmp_path):
        """End-to-end: emitted junction parses to coords that match a GENCODE
        BED reference entry — this is the regression scenario that motivated
        Issue #370 (annotated_discarded was 0 because coords didn't match)."""
        from filter_junctions import _load_reference_junctions, _parse_junction_id

        # GENCODE-style reference BED (0-based half-open)
        ref_bed = tmp_path / "gencode.bed"
        ref_bed.write_text("chr22\t100\t200\tref_junc\t0\t+\n")
        ref = _load_reference_junctions(ref_bed)

        bed = tmp_path / "regtools.bed"
        bed.write_text(_BED12_LINE + "\n")
        out = tmp_path / "raw_junctions.tsv"
        convert_bed12_to_junctions(bed, out)

        junction_id = out.read_text().strip().split("\t")[0]
        parsed = _parse_junction_id(junction_id)
        assert parsed in ref, (
            f"Parsed junction {parsed} not in reference {ref}. "
            "If this fails, the BED12 → raw_junctions.tsv conversion is "
            "emitting anchor outer boundaries instead of intron coords."
        )
