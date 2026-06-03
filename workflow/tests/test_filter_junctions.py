"""Tests for filter_junctions.py — junction origin classification logic."""

import pandas as pd
import pytest

from filter_junctions import (
    _build_cds_donor_lookup,
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

    def test_tumor_exclusive_labeled_correctly(self, tmp_path):
        # Files must live at {sample_id}/raw_junctions.tsv so fp.parent.name == sample_id
        tumor_f = tmp_path / "tumor" / "raw_junctions.tsv"
        normal_f = tmp_path / "normal" / "raw_junctions.tsv"
        tumor_f.parent.mkdir()
        normal_f.parent.mkdir()

        # Junction only in tumor → tumor_exclusive.
        # Add a low-read noise junction so the high-read one clears the mean filter
        # (classify_junctions keeps junctions with reads > mean).
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
        assert df.iloc[0]["junction_origin"] == "tumor_exclusive"

    def test_normal_shared_labeled_correctly(self, tmp_path):
        tumor_f = tmp_path / "tumor" / "raw_junctions.tsv"
        normal_f = tmp_path / "normal" / "raw_junctions.tsv"
        tumor_f.parent.mkdir()
        normal_f.parent.mkdir()

        # Same junction in both tumor and normal → normal_shared
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
        assert df.iloc[0]["junction_origin"] == "normal_shared"

    def test_annotated_junctions_discarded(self, tmp_path):
        tumor_f = tmp_path / "tumor" / "raw_junctions.tsv"
        normal_f = tmp_path / "normal" / "raw_junctions.tsv"
        tumor_f.parent.mkdir()
        normal_f.parent.mkdir()

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

    def test_no_normal_sample_all_labeled_tumor_exclusive(self, tmp_path):
        tumor_f = tmp_path / "tumor" / "raw_junctions.tsv"
        tumor_f.parent.mkdir()

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
        assert df.iloc[0]["junction_origin"] == "tumor_exclusive"
        assert "reading_frame" in df.columns

    # --- GTEx pan-tissue filter (Issue #212) ---------------------------------
    # The GTEx blacklist BED is the same BED6 format as the reference, so the
    # existing _write_reference_bed helper writes it. Junction IDs are 1-based
    # (chr22:201:300:+ → 0-based 200,300), so a matching gtex entry is (200,300).

    def test_gtex_pantissue_shared_labeled_correctly(self, tmp_path):
        # Unannotated junction absent from the matched normal but present in the
        # GTEx blacklist → gtex_pantissue_shared.
        tumor_f = tmp_path / "tumor" / "raw_junctions.tsv"
        normal_f = tmp_path / "normal" / "raw_junctions.tsv"
        tumor_f.parent.mkdir()
        normal_f.parent.mkdir()
        self._write_junction_file(tumor_f, [
            ("chr22:201:300:+", 100),  # in GTEx, not in normal → gtex_pantissue_shared
            ("chr22:901:1000:+", 1),   # noise — below mean, filtered out
        ])
        self._write_junction_file(normal_f, [])

        manifest = tmp_path / "manifest.tsv"
        self._write_manifest(manifest, [
            ("tumor", "Primary Tumor"), ("normal", "Solid Tissue Normal"),
        ])
        ref_bed = tmp_path / "ref.bed"
        self._write_reference_bed(ref_bed, [])
        gtex_bed = tmp_path / "gtex.bed"
        self._write_reference_bed(gtex_bed, [("chr22", 200, 300, "+")])

        output = tmp_path / "novel.tsv"
        classify_junctions(
            junction_files=[tumor_f, normal_f],
            manifest_path=manifest,
            reference_bed=ref_bed,
            output_path=output,
            gtex_bed=gtex_bed,
        )
        df = pd.read_csv(output, sep="\t")
        assert len(df) == 1
        assert df.iloc[0]["junction_origin"] == "gtex_pantissue_shared"

    def test_normal_takes_precedence_over_gtex(self, tmp_path):
        # A junction present in BOTH the matched normal AND the GTEx blacklist is
        # labeled normal_shared (matched normal is the higher-confidence, more
        # specific evidence). The GTEx count is the *marginal* contribution beyond
        # the matched normal, so it does not re-count this junction.
        tumor_f = tmp_path / "tumor" / "raw_junctions.tsv"
        normal_f = tmp_path / "normal" / "raw_junctions.tsv"
        tumor_f.parent.mkdir()
        normal_f.parent.mkdir()
        self._write_junction_file(tumor_f, [
            ("chr22:201:300:+", 100),  # in BOTH normal and gtex
            ("chr22:901:1000:+", 1),   # noise
        ])
        self._write_junction_file(normal_f, [("chr22:201:300:+", 5)])

        manifest = tmp_path / "manifest.tsv"
        self._write_manifest(manifest, [
            ("tumor", "Primary Tumor"), ("normal", "Solid Tissue Normal"),
        ])
        ref_bed = tmp_path / "ref.bed"
        self._write_reference_bed(ref_bed, [])
        gtex_bed = tmp_path / "gtex.bed"
        self._write_reference_bed(gtex_bed, [("chr22", 200, 300, "+")])

        output = tmp_path / "novel.tsv"
        classify_junctions(
            junction_files=[tumor_f, normal_f],
            manifest_path=manifest,
            reference_bed=ref_bed,
            output_path=output,
            gtex_bed=gtex_bed,
        )
        df = pd.read_csv(output, sep="\t")
        assert len(df) == 1
        assert df.iloc[0]["junction_origin"] == "normal_shared"

    def test_gtex_and_tumor_exclusive_coexist(self, tmp_path):
        # Two unannotated junctions absent from the normal: one in the GTEx
        # blacklist (→ gtex_pantissue_shared, excluded downstream) and one in
        # neither (→ tumor_exclusive, the only prediction candidate).
        tumor_f = tmp_path / "tumor" / "raw_junctions.tsv"
        normal_f = tmp_path / "normal" / "raw_junctions.tsv"
        tumor_f.parent.mkdir()
        normal_f.parent.mkdir()
        self._write_junction_file(tumor_f, [
            ("chr22:201:300:+", 100),  # in GTEx → gtex_pantissue_shared
            ("chr22:301:400:+", 100),  # in neither → tumor_exclusive
            ("chr22:901:1000:+", 1),   # noise
        ])
        self._write_junction_file(normal_f, [])

        manifest = tmp_path / "manifest.tsv"
        self._write_manifest(manifest, [
            ("tumor", "Primary Tumor"), ("normal", "Solid Tissue Normal"),
        ])
        ref_bed = tmp_path / "ref.bed"
        self._write_reference_bed(ref_bed, [])
        gtex_bed = tmp_path / "gtex.bed"
        self._write_reference_bed(gtex_bed, [("chr22", 200, 300, "+")])

        output = tmp_path / "novel.tsv"
        classify_junctions(
            junction_files=[tumor_f, normal_f],
            manifest_path=manifest,
            reference_bed=ref_bed,
            output_path=output,
            gtex_bed=gtex_bed,
        )
        df = pd.read_csv(output, sep="\t")
        origins = dict(zip(df["junction_id"], df["junction_origin"]))
        assert origins["chr22:201:300:+"] == "gtex_pantissue_shared"
        assert origins["chr22:301:400:+"] == "tumor_exclusive"
        # Only the tumor_exclusive junction is a candidate downstream.
        assert (df["junction_origin"] == "tumor_exclusive").sum() == 1


# ---------------------------------------------------------------------------
# classify_junctions — junction_filter_stats.tsv (Issue #214)
# ---------------------------------------------------------------------------

class TestClassifyJunctionsStats:
    """Per-sample funnel counts surfaced for report.tsv (tumor samples only)."""

    def _write_junction_file(self, path, rows):
        path.write_text("".join(f"{jid}\t{reads}\n" for jid, reads in rows))

    def _write_manifest(self, path, entries):
        lines = ["file_id\tfile_name\tsample_type\tproject_id\n"]
        for file_id, sample_type in entries:
            lines.append(f"{file_id}\t{file_id}.tsv\t{sample_type}\tlocal\n")
        path.write_text("".join(lines))

    def _write_reference_bed(self, path, junctions):
        lines = [f"{c}\t{s}\t{e}\tjunc\t0\t{st}\n" for c, s, e, st in junctions]
        path.write_text("".join(lines))

    def test_stats_tsv_schema_and_categories(self, tmp_path):
        # Tumor sample with one annotated, one unannotated tumor_exclusive,
        # one unannotated normal_shared. Plus a noise junction below the mean
        # filter so the high-read ones survive.
        tumor_f = tmp_path / "tumor" / "raw_junctions.tsv"
        normal_f = tmp_path / "normal" / "raw_junctions.tsv"
        tumor_f.parent.mkdir()
        normal_f.parent.mkdir()

        self._write_junction_file(tumor_f, [
            ("chr22:101:200:+", 100),  # annotated → discarded
            ("chr22:201:300:+", 100),  # unannotated, in normal → normal_shared
            ("chr22:301:400:+", 100),  # unannotated, NOT in normal → tumor_exclusive
            ("chr22:901:1000:+", 1),   # noise — below mean, filtered out
        ])
        self._write_junction_file(normal_f, [("chr22:201:300:+", 5)])

        manifest = tmp_path / "manifest.tsv"
        self._write_manifest(manifest, [
            ("tumor", "Primary Tumor"),
            ("normal", "Solid Tissue Normal"),
        ])

        ref_bed = tmp_path / "ref.bed"
        self._write_reference_bed(ref_bed, [("chr22", 100, 200, "+")])

        output = tmp_path / "novel.tsv"
        stats_output = tmp_path / "junction_filter_stats.tsv"
        classify_junctions(
            junction_files=[tumor_f, normal_f],
            manifest_path=manifest,
            reference_bed=ref_bed,
            output_path=output,
            stats_output_path=stats_output,
        )

        stats = pd.read_csv(stats_output, sep="\t")
        assert set(stats.columns) == {"sample_id", "sample_type", "category", "count"}

        # Long-format per tumor sample: 6 funnel rows (Issue #214 + the
        # gtex_pantissue_shared row from Issue #212, always emitted) +
        # 4 descriptive distribution rows (Issue #215) = 10 total.
        assert len(stats) == 10
        assert (stats["sample_id"] == "tumor").all()
        assert (stats["sample_type"] == "Primary Tumor").all()

        by_cat = dict(zip(stats["category"], stats["count"]))
        # junctions_raw is BEFORE the mean-reads filter — all 4 rows
        assert by_cat["junctions_raw"] == 4
        # 1 noise junction was below the mean and got filtered out
        assert by_cat["mean_reads_filtered"] == 1
        # After mean filter, 3 junctions survive: 1 annotated, 1 normal_shared, 1 tumor_exclusive
        assert by_cat["annotated_discarded"] == 1
        assert by_cat["normal_shared"] == 1
        assert by_cat["tumor_exclusive"] == 1
        # No GTEx blacklist supplied → the row is still emitted with count 0.
        assert by_cat["gtex_pantissue_shared"] == 0

    def test_funnel_reconciles_arithmetically(self, tmp_path):
        """junctions_raw must equal the sum of the downstream partition categories."""
        tumor_f = tmp_path / "tumor" / "raw_junctions.tsv"
        normal_f = tmp_path / "normal" / "raw_junctions.tsv"
        tumor_f.parent.mkdir()
        normal_f.parent.mkdir()

        # Mix of high-read and low-read junctions to exercise the mean filter
        self._write_junction_file(tumor_f, [
            ("chr22:101:200:+", 100),  # annotated, high read
            ("chr22:201:300:+", 100),  # unannotated, high read, in normal → normal_shared
            ("chr22:301:400:+", 100),  # unannotated, high read → tumor_exclusive
            ("chr22:701:800:+", 1),    # noise → mean_reads_filtered
            ("chr22:801:900:+", 1),    # noise → mean_reads_filtered
        ])
        self._write_junction_file(normal_f, [("chr22:201:300:+", 5)])

        manifest = tmp_path / "manifest.tsv"
        self._write_manifest(manifest, [
            ("tumor", "Primary Tumor"),
            ("normal", "Solid Tissue Normal"),
        ])

        ref_bed = tmp_path / "ref.bed"
        self._write_reference_bed(ref_bed, [("chr22", 100, 200, "+")])

        output = tmp_path / "novel.tsv"
        stats_output = tmp_path / "junction_filter_stats.tsv"
        classify_junctions(
            junction_files=[tumor_f, normal_f],
            manifest_path=manifest,
            reference_bed=ref_bed,
            output_path=output,
            stats_output_path=stats_output,
        )

        stats = pd.read_csv(stats_output, sep="\t")
        by_cat = dict(zip(stats["category"], stats["count"]))
        # The whole point of the funnel: raw == sum of downstream buckets
        assert by_cat["junctions_raw"] == (
            by_cat["mean_reads_filtered"]
            + by_cat["annotated_discarded"]
            + by_cat["normal_shared"]
            + by_cat["gtex_pantissue_shared"]
            + by_cat["tumor_exclusive"]
        )

    def test_gtex_pantissue_shared_counted_and_reconciles(self, tmp_path):
        """With a GTEx blacklist supplied, the gtex_pantissue_shared row carries the
        marginal count and the funnel (now 6 buckets) still reconciles."""
        tumor_f = tmp_path / "tumor" / "raw_junctions.tsv"
        normal_f = tmp_path / "normal" / "raw_junctions.tsv"
        tumor_f.parent.mkdir()
        normal_f.parent.mkdir()
        self._write_junction_file(tumor_f, [
            ("chr22:101:200:+", 100),  # annotated → discarded
            ("chr22:201:300:+", 100),  # in normal → normal_shared
            ("chr22:301:400:+", 100),  # in GTEx, not normal → gtex_pantissue_shared
            ("chr22:401:500:+", 100),  # in neither → tumor_exclusive
            ("chr22:901:1000:+", 1),   # noise → mean_reads_filtered
        ])
        self._write_junction_file(normal_f, [("chr22:201:300:+", 5)])

        manifest = tmp_path / "manifest.tsv"
        self._write_manifest(manifest, [
            ("tumor", "Primary Tumor"), ("normal", "Solid Tissue Normal"),
        ])
        ref_bed = tmp_path / "ref.bed"
        self._write_reference_bed(ref_bed, [("chr22", 100, 200, "+")])
        gtex_bed = tmp_path / "gtex.bed"
        self._write_reference_bed(gtex_bed, [("chr22", 300, 400, "+")])

        output = tmp_path / "novel.tsv"
        stats_output = tmp_path / "junction_filter_stats.tsv"
        classify_junctions(
            junction_files=[tumor_f, normal_f],
            manifest_path=manifest,
            reference_bed=ref_bed,
            output_path=output,
            stats_output_path=stats_output,
            gtex_bed=gtex_bed,
        )
        stats = pd.read_csv(stats_output, sep="\t")
        by_cat = dict(zip(stats["category"], stats["count"]))
        assert by_cat["normal_shared"] == 1
        assert by_cat["gtex_pantissue_shared"] == 1
        assert by_cat["tumor_exclusive"] == 1
        assert by_cat["junctions_raw"] == (
            by_cat["mean_reads_filtered"]
            + by_cat["annotated_discarded"]
            + by_cat["normal_shared"]
            + by_cat["gtex_pantissue_shared"]
            + by_cat["tumor_exclusive"]
        )

    def test_stats_tsv_omits_normal_samples(self, tmp_path):
        tumor_f = tmp_path / "tumor" / "raw_junctions.tsv"
        normal_f = tmp_path / "normal" / "raw_junctions.tsv"
        tumor_f.parent.mkdir()
        normal_f.parent.mkdir()

        self._write_junction_file(tumor_f, [
            ("chr22:201:300:+", 100),
            ("chr22:401:500:+", 1),
        ])
        self._write_junction_file(normal_f, [
            ("chr22:201:300:+", 5),
            ("chr22:601:700:+", 5),
        ])

        manifest = tmp_path / "manifest.tsv"
        self._write_manifest(manifest, [
            ("tumor", "Primary Tumor"),
            ("normal", "Solid Tissue Normal"),
        ])

        ref_bed = tmp_path / "ref.bed"
        self._write_reference_bed(ref_bed, [])

        output = tmp_path / "novel.tsv"
        stats_output = tmp_path / "junction_filter_stats.tsv"
        classify_junctions(
            junction_files=[tumor_f, normal_f],
            manifest_path=manifest,
            reference_bed=ref_bed,
            output_path=output,
            stats_output_path=stats_output,
        )

        stats = pd.read_csv(stats_output, sep="\t")
        # No rows for the normal sample even though it has 2 raw junctions
        assert "normal" not in set(stats["sample_id"])

    def test_multi_tumor_samples_each_get_their_own_rows(self, tmp_path):
        tumor1 = tmp_path / "tumor1" / "raw_junctions.tsv"
        tumor2 = tmp_path / "tumor2" / "raw_junctions.tsv"
        tumor1.parent.mkdir()
        tumor2.parent.mkdir()

        self._write_junction_file(tumor1, [
            ("chr22:201:300:+", 100),
            ("chr22:401:500:+", 1),
        ])
        self._write_junction_file(tumor2, [
            ("chr22:301:400:+", 100),
            ("chr22:501:600:+", 1),
        ])

        manifest = tmp_path / "manifest.tsv"
        self._write_manifest(manifest, [
            ("tumor1", "Primary Tumor"),
            ("tumor2", "Primary Tumor"),
        ])

        ref_bed = tmp_path / "ref.bed"
        self._write_reference_bed(ref_bed, [])

        output = tmp_path / "novel.tsv"
        stats_output = tmp_path / "junction_filter_stats.tsv"
        classify_junctions(
            junction_files=[tumor1, tumor2],
            manifest_path=manifest,
            reference_bed=ref_bed,
            output_path=output,
            stats_output_path=stats_output,
        )

        stats = pd.read_csv(stats_output, sep="\t")
        # 10 categories per tumor sample (6 funnel incl. gtex_pantissue_shared +
        # 4 distribution) × 2 samples = 20 rows
        assert len(stats) == 20
        assert set(stats["sample_id"]) == {"tumor1", "tumor2"}

    def test_stats_tsv_optional(self, tmp_path):
        """Existing callers that don't pass stats_output_path still work."""
        tumor_f = tmp_path / "tumor" / "raw_junctions.tsv"
        tumor_f.parent.mkdir()
        self._write_junction_file(tumor_f, [
            ("chr22:201:300:+", 100),
            ("chr22:401:500:+", 1),
        ])
        manifest = tmp_path / "manifest.tsv"
        self._write_manifest(manifest, [("tumor", "Primary Tumor")])
        ref_bed = tmp_path / "ref.bed"
        self._write_reference_bed(ref_bed, [])

        output = tmp_path / "novel.tsv"
        # No stats_output_path — should not raise
        classify_junctions(
            junction_files=[tumor_f],
            manifest_path=manifest,
            reference_bed=ref_bed,
            output_path=output,
        )
        assert output.exists()

    def test_stats_tsv_includes_read_distribution(self, tmp_path):
        """Issue #215 — 4 descriptive distribution rows on the raw read counts.

        These rows let the Scientist sanity-check the silent per-file mean
        threshold against the sample's actual read-count distribution.
        """
        tumor_f = tmp_path / "tumor" / "raw_junctions.tsv"
        tumor_f.parent.mkdir()
        # Read counts: 1, 5, 10, 100. min=1, median=7.5, mean=29.0, max=100.
        self._write_junction_file(tumor_f, [
            ("chr22:101:200:+", 1),
            ("chr22:201:300:+", 5),
            ("chr22:301:400:+", 10),
            ("chr22:401:500:+", 100),
        ])
        manifest = tmp_path / "manifest.tsv"
        self._write_manifest(manifest, [("tumor", "Primary Tumor")])
        ref_bed = tmp_path / "ref.bed"
        self._write_reference_bed(ref_bed, [])

        output = tmp_path / "novel.tsv"
        stats_output = tmp_path / "junction_filter_stats.tsv"
        classify_junctions(
            junction_files=[tumor_f],
            manifest_path=manifest,
            reference_bed=ref_bed,
            output_path=output,
            stats_output_path=stats_output,
        )

        stats = pd.read_csv(stats_output, sep="\t")
        by_cat = dict(zip(stats["category"], stats["count"]))
        assert by_cat["min_reads"] == 1
        assert by_cat["median_reads"] == 7.5
        assert by_cat["mean_reads"] == 29.0
        assert by_cat["max_reads"] == 100


# ---------------------------------------------------------------------------
# _build_cds_donor_lookup
# ---------------------------------------------------------------------------

def _write_gtf(path, records):
    """Write minimal GTF records for testing the CDS donor lookup."""
    lines = []
    for r in records:
        tx_type = r.get("tx_type", "protein_coding")
        tx_id = r.get("tx_id", "ENST001")
        attrs = (
            f'gene_id "{r["gene_id"]}"; transcript_id "{tx_id}"; '
            f'transcript_type "{tx_type}";'
        )
        lines.append(
            f'{r["chrom"]}\tHAVANA\t{r["feature"]}\t{r["start"]}\t{r["end"]}'
            f'\t.\t{r["strand"]}\t{r["frame"]}\t{attrs}\n'
        )
    path.write_text("".join(lines))


class TestBuildCdsDonorLookup:
    def test_plus_strand_clean_boundary(self, tmp_path):
        # exon len=3, frame=0 → phase=0 → frame_offset=0; donor_coord=end0=102
        gtf = tmp_path / "test.gtf"
        _write_gtf(gtf, [{"chrom": "chr1", "feature": "CDS", "start": 100, "end": 102,
                           "strand": "+", "frame": 0, "gene_id": "G1"}])
        lookup = _build_cds_donor_lookup(gtf)
        assert ("chr1", 102, "+") in lookup
        assert lookup[("chr1", 102, "+")] == {0}

    def test_plus_strand_phase1_gives_frame_offset_2(self, tmp_path):
        # exon len=4, frame=0 → phase=1 → frame_offset=2; donor=103
        gtf = tmp_path / "test.gtf"
        _write_gtf(gtf, [{"chrom": "chr1", "feature": "CDS", "start": 100, "end": 103,
                           "strand": "+", "frame": 0, "gene_id": "G1"}])
        lookup = _build_cds_donor_lookup(gtf)
        assert lookup[("chr1", 103, "+")] == {2}

    def test_plus_strand_phase2_gives_frame_offset_1(self, tmp_path):
        # exon len=5, frame=0 → phase=2 → frame_offset=1; donor=104
        gtf = tmp_path / "test.gtf"
        _write_gtf(gtf, [{"chrom": "chr1", "feature": "CDS", "start": 100, "end": 104,
                           "strand": "+", "frame": 0, "gene_id": "G1"}])
        lookup = _build_cds_donor_lookup(gtf)
        assert lookup[("chr1", 104, "+")] == {1}

    def test_nonzero_gtf_frame_shifts_phase(self, tmp_path):
        # exon len=4, frame=1 → phase=(4-1)%3=0 → frame_offset=0
        gtf = tmp_path / "test.gtf"
        _write_gtf(gtf, [{"chrom": "chr1", "feature": "CDS", "start": 100, "end": 103,
                           "strand": "+", "frame": 1, "gene_id": "G1"}])
        lookup = _build_cds_donor_lookup(gtf)
        assert lookup[("chr1", 103, "+")] == {0}

    def test_minus_strand_donor_at_start0(self, tmp_path):
        # − strand: donor_coord = start0 = int(start)-1 = 99; end0=102
        gtf = tmp_path / "test.gtf"
        _write_gtf(gtf, [{"chrom": "chr1", "feature": "CDS", "start": 100, "end": 102,
                           "strand": "-", "frame": 0, "gene_id": "G1"}])
        lookup = _build_cds_donor_lookup(gtf)
        assert ("chr1", 99, "-") in lookup
        assert ("chr1", 102, "-") not in lookup

    def test_multi_transcript_same_frame_single_entry(self, tmp_path):
        # Two transcripts, same donor and same frame → set with one element
        gtf = tmp_path / "test.gtf"
        _write_gtf(gtf, [
            {"chrom": "chr1", "feature": "CDS", "start": 100, "end": 102,
             "strand": "+", "frame": 0, "gene_id": "G1", "tx_id": "T1"},
            {"chrom": "chr1", "feature": "CDS", "start": 100, "end": 102,
             "strand": "+", "frame": 0, "gene_id": "G1", "tx_id": "T2"},
        ])
        lookup = _build_cds_donor_lookup(gtf)
        assert lookup[("chr1", 102, "+")] == {0}

    def test_two_genes_different_frames_at_same_donor_kept_as_union(self, tmp_path):
        # Two genes, same donor coord 102 (end0), different frame offsets → both in set
        # G1: start=100 (1-based) → start0=99, end0=102, len=3, phase=0, offset=0
        # G2: start=99  (1-based) → start0=98, end0=102, len=4, phase=1, offset=2
        gtf = tmp_path / "test.gtf"
        _write_gtf(gtf, [
            {"chrom": "chr1", "feature": "CDS", "start": 100, "end": 102,
             "strand": "+", "frame": 0, "gene_id": "G1"},
            {"chrom": "chr1", "feature": "CDS", "start": 99, "end": 102,
             "strand": "+", "frame": 0, "gene_id": "G2"},
        ])
        lookup = _build_cds_donor_lookup(gtf)
        assert lookup[("chr1", 102, "+")] == {0, 2}

    def test_non_protein_coding_excluded(self, tmp_path):
        gtf = tmp_path / "test.gtf"
        _write_gtf(gtf, [{"chrom": "chr1", "feature": "CDS", "start": 100, "end": 102,
                           "strand": "+", "frame": 0, "gene_id": "G1",
                           "tx_type": "nonsense_mediated_decay"}])
        lookup = _build_cds_donor_lookup(gtf)
        assert ("chr1", 102, "+") not in lookup

    def test_non_cds_feature_excluded(self, tmp_path):
        gtf = tmp_path / "test.gtf"
        _write_gtf(gtf, [{"chrom": "chr1", "feature": "exon", "start": 100, "end": 102,
                           "strand": "+", "frame": 0, "gene_id": "G1"}])
        lookup = _build_cds_donor_lookup(gtf)
        assert len(lookup) == 0


# ---------------------------------------------------------------------------
# classify_junctions — reading_frame column
# ---------------------------------------------------------------------------

class TestClassifyJunctionsReadingFrame:
    def _write_junction_file(self, path, rows):
        path.write_text("".join(f"{jid}\t{reads}\n" for jid, reads in rows))

    def _write_manifest(self, path, entries):
        lines = ["file_id\tfile_name\tsample_type\tproject_id\n"]
        for fid, stype in entries:
            lines.append(f"{fid}\t{fid}.tsv\t{stype}\tlocal\n")
        path.write_text("".join(lines))

    def _read_output(self, path):
        return pd.read_csv(path, sep="\t", dtype={"reading_frame": str}, keep_default_na=False)

    def test_reading_frame_annotated_when_gtf_supplied(self, tmp_path):
        # Junction chr1:201:300:+ → donor (0-based start) = 200
        # GTF CDS: start=198 (1-based) → start0=197, end0=200, len=3, frame=0
        # phase=0, frame_offset=0
        tumor_f = tmp_path / "tumor" / "raw_junctions.tsv"
        tumor_f.parent.mkdir()
        self._write_junction_file(tumor_f, [("chr1:201:300:+", 100), ("chr1:401:500:+", 1)])
        manifest = tmp_path / "manifest.tsv"
        self._write_manifest(manifest, [("tumor", "Primary Tumor")])
        ref_bed = tmp_path / "ref.bed"
        ref_bed.write_text("")
        gtf = tmp_path / "test.gtf"
        _write_gtf(gtf, [{"chrom": "chr1", "feature": "CDS", "start": 198, "end": 200,
                           "strand": "+", "frame": 0, "gene_id": "G1"}])
        output = tmp_path / "novel.tsv"
        classify_junctions(
            junction_files=[tumor_f], manifest_path=manifest,
            reference_bed=ref_bed, output_path=output, gencode_gtf=gtf,
        )
        df = self._read_output(output)
        assert df.iloc[0]["reading_frame"] == "0"

    def test_reading_frame_empty_when_no_cds_match(self, tmp_path):
        tumor_f = tmp_path / "tumor" / "raw_junctions.tsv"
        tumor_f.parent.mkdir()
        self._write_junction_file(tumor_f, [("chr1:201:300:+", 100), ("chr1:401:500:+", 1)])
        manifest = tmp_path / "manifest.tsv"
        self._write_manifest(manifest, [("tumor", "Primary Tumor")])
        ref_bed = tmp_path / "ref.bed"
        ref_bed.write_text("")
        gtf = tmp_path / "empty.gtf"
        gtf.write_text("")
        output = tmp_path / "novel.tsv"
        classify_junctions(
            junction_files=[tumor_f], manifest_path=manifest,
            reference_bed=ref_bed, output_path=output, gencode_gtf=gtf,
        )
        df = self._read_output(output)
        assert df.iloc[0]["reading_frame"] == ""

    def test_reading_frame_empty_when_no_gtf(self, tmp_path):
        tumor_f = tmp_path / "tumor" / "raw_junctions.tsv"
        tumor_f.parent.mkdir()
        self._write_junction_file(tumor_f, [("chr1:201:300:+", 100), ("chr1:401:500:+", 1)])
        manifest = tmp_path / "manifest.tsv"
        self._write_manifest(manifest, [("tumor", "Primary Tumor")])
        ref_bed = tmp_path / "ref.bed"
        ref_bed.write_text("")
        output = tmp_path / "novel.tsv"
        classify_junctions(
            junction_files=[tumor_f], manifest_path=manifest,
            reference_bed=ref_bed, output_path=output,
        )
        df = self._read_output(output)
        assert df.iloc[0]["reading_frame"] == ""

    def test_multi_frame_union_written_sorted(self, tmp_path):
        # Two genes share donor coord 200 (end0); different frame offsets → union
        # G1: start=99 (1-based) → start0=98, end0=200, len=102, phase=0, offset=0
        # G2: start=100 (1-based) → start0=99, end0=200, len=101, phase=2, offset=1
        tumor_f = tmp_path / "tumor" / "raw_junctions.tsv"
        tumor_f.parent.mkdir()
        self._write_junction_file(tumor_f, [("chr1:201:300:+", 100), ("chr1:401:500:+", 1)])
        manifest = tmp_path / "manifest.tsv"
        self._write_manifest(manifest, [("tumor", "Primary Tumor")])
        ref_bed = tmp_path / "ref.bed"
        ref_bed.write_text("")
        gtf = tmp_path / "test.gtf"
        _write_gtf(gtf, [
            {"chrom": "chr1", "feature": "CDS", "start": 99, "end": 200,
             "strand": "+", "frame": 0, "gene_id": "G1"},
            {"chrom": "chr1", "feature": "CDS", "start": 100, "end": 200,
             "strand": "+", "frame": 0, "gene_id": "G2"},
        ])
        output = tmp_path / "novel.tsv"
        classify_junctions(
            junction_files=[tumor_f], manifest_path=manifest,
            reference_bed=ref_bed, output_path=output, gencode_gtf=gtf,
        )
        df = self._read_output(output)
        assert set(df.iloc[0]["reading_frame"].split(",")) == {"0", "1"}

    def test_minus_strand_reading_frame_annotated(self, tmp_path):
        # − strand: donor_coord = junction.end = cds_exon_start0
        # Junction chr1:100:201:- → end=201, donor_coord=201
        # CDS: start=202 (1-based) → start0=201, end0=300, len=99, frame=0
        # phase=(99-0)%3=0, offset=(-0)%3=0 → reading_frame "0"
        tumor_f = tmp_path / "tumor" / "raw_junctions.tsv"
        tumor_f.parent.mkdir()
        self._write_junction_file(tumor_f, [("chr1:100:201:-", 100), ("chr1:400:500:-", 1)])
        manifest = tmp_path / "manifest.tsv"
        self._write_manifest(manifest, [("tumor", "Primary Tumor")])
        ref_bed = tmp_path / "ref.bed"
        ref_bed.write_text("")
        gtf = tmp_path / "test.gtf"
        _write_gtf(gtf, [
            {"chrom": "chr1", "feature": "CDS", "start": 202, "end": 300,
             "strand": "-", "frame": 0, "gene_id": "G1"},
        ])
        output = tmp_path / "novel.tsv"
        classify_junctions(
            junction_files=[tumor_f], manifest_path=manifest,
            reference_bed=ref_bed, output_path=output, gencode_gtf=gtf,
        )
        df = self._read_output(output)
        matched = df[df["junction_origin"] == "tumor_exclusive"]
        assert matched.iloc[0]["reading_frame"] == "0"
