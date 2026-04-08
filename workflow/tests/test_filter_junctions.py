"""Tests for filter_junctions.py — junction origin classification logic."""

import pandas as pd
import pytest

from filter_junctions import (
    _build_normal_junction_set,
    _load_reference_junctions,
    _parse_junction_id,
    classify_junctions,
)


# ---------------------------------------------------------------------------
# _parse_junction_id
# ---------------------------------------------------------------------------

class TestParseJunctionId:
    def test_colon_separator(self):
        # Standard format: chr:start:end:strand (1-based start → 0-based)
        assert _parse_junction_id("chr22:101:200:+") == ("chr22", 100, 200, "+")

    def test_minus_strand(self):
        assert _parse_junction_id("chr1:51:150:-") == ("chr1", 50, 150, "-")

    def test_dash_separator_fallback(self):
        # Some tools emit chr-start-end-strand
        assert _parse_junction_id("chr22-101-200-+") == ("chr22", 100, 200, "+")

    def test_invalid_strand_becomes_dot(self):
        result = _parse_junction_id("chr22:101:200:?")
        assert result is not None
        assert result[3] == "."

    def test_too_few_parts_returns_none(self):
        assert _parse_junction_id("chr22:100") is None

    def test_non_numeric_coords_returns_none(self):
        assert _parse_junction_id("chr22:abc:200:+") is None

    def test_start_converts_to_zero_based(self):
        # 1-based start 1 → 0-based start 0
        chrom, start, end, strand = _parse_junction_id("chr22:1:100:+")
        assert start == 0
        assert end == 100


# ---------------------------------------------------------------------------
# _load_reference_junctions
# ---------------------------------------------------------------------------

class TestLoadReferenceJunctions:
    def test_loads_bed_correctly(self, tmp_path):
        bed = tmp_path / "ref.bed"
        # BED6: chrom, start, end, name, score, strand
        bed.write_text(
            "chr22\t100\t200\tjunc1\t0\t+\n"
            "chr22\t300\t400\tjunc2\t0\t-\n"
        )
        ref = _load_reference_junctions(bed)
        assert ("chr22", 100, 200, "+") in ref
        assert ("chr22", 300, 400, "-") in ref
        assert len(ref) == 2

    def test_skips_comment_lines(self, tmp_path):
        bed = tmp_path / "ref.bed"
        bed.write_text(
            "# comment\n"
            "chr22\t100\t200\tjunc1\t0\t+\n"
        )
        ref = _load_reference_junctions(bed)
        assert len(ref) == 1

    def test_skips_short_lines(self, tmp_path):
        bed = tmp_path / "ref.bed"
        bed.write_text("chr22\t100\t200\n")  # only 3 columns
        ref = _load_reference_junctions(bed)
        assert len(ref) == 0


# ---------------------------------------------------------------------------
# _build_normal_junction_set
# ---------------------------------------------------------------------------

class TestBuildNormalJunctionSet:
    def _write_junction_file(self, path, rows):
        """Write a junction TSV: junction_id<TAB>reads."""
        path.write_text("".join(f"{jid}\t{reads}\n" for jid, reads in rows))

    def test_includes_junctions_above_min_reads(self, tmp_path):
        f = tmp_path / "normal.tsv"
        self._write_junction_file(f, [("chr22:101:200:+", 5)])
        ref = frozenset()
        result = _build_normal_junction_set([(f, "normal")], ref, min_reads=2)
        assert ("chr22", 100, 200, "+") in result

    def test_excludes_junctions_below_min_reads(self, tmp_path):
        f = tmp_path / "normal.tsv"
        self._write_junction_file(f, [("chr22:101:200:+", 1)])
        ref = frozenset()
        result = _build_normal_junction_set([(f, "normal")], ref, min_reads=2)
        assert ("chr22", 100, 200, "+") not in result

    def test_excludes_reference_junctions(self, tmp_path):
        f = tmp_path / "normal.tsv"
        self._write_junction_file(f, [("chr22:101:200:+", 5)])
        ref = frozenset({("chr22", 100, 200, "+")})
        result = _build_normal_junction_set([(f, "normal")], ref, min_reads=2)
        assert len(result) == 0

    def test_empty_normal_files_returns_empty_set(self):
        result = _build_normal_junction_set([], frozenset(), min_reads=2)
        assert len(result) == 0


# ---------------------------------------------------------------------------
# classify_junctions (integration)
# ---------------------------------------------------------------------------

class TestClassifyJunctions:
    def _write_junction_file(self, path, rows):
        path.write_text("".join(f"{jid}\t{reads}\n" for jid, reads in rows))

    def _write_manifest(self, path, entries):
        """entries: list of (file_id, sample_type)"""
        lines = ["file_id\tfile_name\tsample_type\tproject_id\n"]
        for file_id, sample_type in entries:
            lines.append(f"{file_id}\t{file_id}.tsv\t{sample_type}\tlocal\n")
        path.write_text("".join(lines))

    def _write_reference_bed(self, path, junctions):
        """junctions: list of (chrom, start, end, strand)"""
        lines = [f"{c}\t{s}\t{e}\tjunc\t0\t{st}\n" for c, s, e, st in junctions]
        path.write_text("".join(lines))

    def test_tumor_specific_labeled_correctly(self, tmp_path):
        files_dir = tmp_path / "files"
        files_dir.mkdir()

        # Junction only in tumor → tumor_specific.
        # Add a low-read noise junction so the high-read one clears the mean filter
        # (classify_junctions keeps junctions with reads > mean).
        tumor_f = files_dir / "tumor.tsv"
        normal_f = files_dir / "normal.tsv"
        self._write_junction_file(tumor_f, [
            ("chr22:201:300:+", 100),  # unannotated, absent in normal
            ("chr22:401:500:+", 1),    # noise — below mean, filtered out
        ])
        self._write_junction_file(normal_f, [])

        manifest = tmp_path / "manifest.tsv"
        self._write_manifest(manifest, [
            ("tumor", "Primary Tumor"),
            ("normal", "Solid Tissue Normal"),
        ])

        ref_bed = tmp_path / "ref.bed"
        self._write_reference_bed(ref_bed, [("chr22", 100, 200, "+")])

        output = tmp_path / "novel.tsv"
        classify_junctions(
            junction_files=[tumor_f, normal_f],
            manifest_path=manifest,
            reference_bed=ref_bed,
            output_path=output,
        )

        df = pd.read_csv(output, sep="\t")
        assert len(df) == 1
        assert df.iloc[0]["junction_origin"] == "tumor_specific"

    def test_patient_specific_labeled_correctly(self, tmp_path):
        files_dir = tmp_path / "files"
        files_dir.mkdir()

        # Same junction in both tumor and normal → patient_specific
        tumor_f = files_dir / "tumor.tsv"
        normal_f = files_dir / "normal.tsv"
        self._write_junction_file(tumor_f, [
            ("chr22:201:300:+", 100),
            ("chr22:401:500:+", 1),  # noise — below mean, filtered out
        ])
        self._write_junction_file(normal_f, [("chr22:201:300:+", 5)])

        manifest = tmp_path / "manifest.tsv"
        self._write_manifest(manifest, [
            ("tumor", "Primary Tumor"),
            ("normal", "Solid Tissue Normal"),
        ])

        ref_bed = tmp_path / "ref.bed"
        self._write_reference_bed(ref_bed, [])

        output = tmp_path / "novel.tsv"
        classify_junctions(
            junction_files=[tumor_f, normal_f],
            manifest_path=manifest,
            reference_bed=ref_bed,
            output_path=output,
        )

        df = pd.read_csv(output, sep="\t")
        assert len(df) == 1
        assert df.iloc[0]["junction_origin"] == "patient_specific"

    def test_annotated_junctions_discarded(self, tmp_path):
        files_dir = tmp_path / "files"
        files_dir.mkdir()

        tumor_f = files_dir / "tumor.tsv"
        normal_f = files_dir / "normal.tsv"
        self._write_junction_file(tumor_f, [
            ("chr22:101:200:+", 100),
            ("chr22:401:500:+", 1),  # noise — keeps annotated junction above mean
        ])
        self._write_junction_file(normal_f, [])

        manifest = tmp_path / "manifest.tsv"
        self._write_manifest(manifest, [
            ("tumor", "Primary Tumor"),
            ("normal", "Solid Tissue Normal"),
        ])

        ref_bed = tmp_path / "ref.bed"
        # Junction chr22:100:200:+ is annotated (0-based = 100, 1-based start = 101)
        self._write_reference_bed(ref_bed, [("chr22", 100, 200, "+")])

        output = tmp_path / "novel.tsv"
        classify_junctions(
            junction_files=[tumor_f, normal_f],
            manifest_path=manifest,
            reference_bed=ref_bed,
            output_path=output,
        )

        df = pd.read_csv(output, sep="\t")
        assert len(df) == 0

    def test_no_normal_sample_all_labeled_tumor_specific(self, tmp_path):
        files_dir = tmp_path / "files"
        files_dir.mkdir()

        tumor_f = files_dir / "tumor.tsv"
        self._write_junction_file(tumor_f, [
            ("chr22:201:300:+", 100),
            ("chr22:401:500:+", 1),  # noise — below mean, filtered out
        ])

        manifest = tmp_path / "manifest.tsv"
        self._write_manifest(manifest, [("tumor", "Primary Tumor")])

        ref_bed = tmp_path / "ref.bed"
        self._write_reference_bed(ref_bed, [])

        output = tmp_path / "novel.tsv"
        classify_junctions(
            junction_files=[tumor_f],
            manifest_path=manifest,
            reference_bed=ref_bed,
            output_path=output,
        )

        df = pd.read_csv(output, sep="\t")
        assert len(df) == 1
        assert df.iloc[0]["junction_origin"] == "tumor_specific"
