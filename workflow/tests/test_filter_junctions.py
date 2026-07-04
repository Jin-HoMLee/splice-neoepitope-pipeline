"""Tests for filter_junctions.py — junction origin classification logic."""

import pandas as pd
import pytest

from filter_junctions import (
    _build_cds_donor_lookup,
    _build_normal_junction_sources,
    _load_reference_junctions,
    _parse_junction_id,
    _read_junction_file,
    classify_junctions,
)


# ---------------------------------------------------------------------------
# _read_junction_file (Issue #375 — optional 3rd annotated-flag column)
# ---------------------------------------------------------------------------

class TestReadJunctionFile:
    def test_two_column_legacy_annotated_none(self, tmp_path):
        # HISAT2/regtools path (bed12_to_junctions.py) emits only 2 columns.
        f = tmp_path / "raw_junctions.tsv"
        f.write_text("chr22:101:200:+\t10\nchr22:301:400:-\t5\n")
        rows = _read_junction_file(f)
        assert rows == [
            ("chr22:101:200:+", 10.0, None),
            ("chr22:301:400:-", 5.0, None),
        ]

    def test_three_column_star_annotated_parsed(self, tmp_path):
        # STAR path now carries the annotated flag as a 3rd column.
        f = tmp_path / "raw_junctions.tsv"
        f.write_text("chr22:101:200:+\t10\t1\nchr22:301:400:-\t5\t0\n")
        rows = _read_junction_file(f)
        assert rows == [
            ("chr22:101:200:+", 10.0, 1),
            ("chr22:301:400:-", 5.0, 0),
        ]

    def test_three_column_non_integer_flag_falls_back_to_none(self, tmp_path):
        f = tmp_path / "raw_junctions.tsv"
        f.write_text("chr22:101:200:+\t10\tnotanint\n")
        rows = _read_junction_file(f)
        assert rows == [("chr22:101:200:+", 10.0, None)]

    def test_skips_comment_and_blank_lines(self, tmp_path):
        f = tmp_path / "raw_junctions.tsv"
        f.write_text("# header\n\nchr22:101:200:+\t10\t1\n")
        rows = _read_junction_file(f)
        assert rows == [("chr22:101:200:+", 10.0, 1)]


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
# _build_normal_junction_sources
# ---------------------------------------------------------------------------

class TestBuildNormalJunctionSources:
    def _write_junction_file(self, path, rows):
        """Write a junction TSV: junction_id<TAB>reads."""
        path.write_text("".join(f"{jid}\t{reads}\n" for jid, reads in rows))

    def test_includes_junctions_above_min_reads(self, tmp_path):
        f = tmp_path / "normal.tsv"
        self._write_junction_file(f, [("chr22:101:200:+", 5)])
        ref = frozenset()
        result = _build_normal_junction_sources([(f, "normal", "Solid Tissue Normal")], ref, min_reads=2)
        assert ("chr22", 100, 200, "+") in result
        # provenance: the sample_type that contained it is recorded
        assert result[("chr22", 100, 200, "+")] == {"Solid Tissue Normal"}

    def test_excludes_junctions_below_min_reads(self, tmp_path):
        f = tmp_path / "normal.tsv"
        self._write_junction_file(f, [("chr22:101:200:+", 1)])
        ref = frozenset()
        result = _build_normal_junction_sources([(f, "normal", "Solid Tissue Normal")], ref, min_reads=2)
        assert ("chr22", 100, 200, "+") not in result

    def test_excludes_reference_junctions(self, tmp_path):
        f = tmp_path / "normal.tsv"
        self._write_junction_file(f, [("chr22:101:200:+", 5)])
        ref = frozenset({("chr22", 100, 200, "+")})
        result = _build_normal_junction_sources([(f, "normal", "Solid Tissue Normal")], ref, min_reads=2)
        assert len(result) == 0

    def test_empty_normal_files_returns_empty(self):
        result = _build_normal_junction_sources([], frozenset(), min_reads=2)
        assert len(result) == 0

    def test_junction_in_two_normal_types_records_both(self, tmp_path):
        solid = tmp_path / "solid.tsv"
        blood = tmp_path / "blood.tsv"
        self._write_junction_file(solid, [("chr22:101:200:+", 5)])
        self._write_junction_file(blood, [("chr22:101:200:+", 4)])
        result = _build_normal_junction_sources(
            [(solid, "n1", "Solid Tissue Normal"), (blood, "n2", "Blood Derived Normal")],
            frozenset(), min_reads=2,
        )
        assert result[("chr22", 100, 200, "+")] == {"Solid Tissue Normal", "Blood Derived Normal"}


# ---------------------------------------------------------------------------
# classify_junctions — STAR annotated-flag cross-check (Issue #375)
# ---------------------------------------------------------------------------

class TestAnnotatedFlagCrossCheck:
    """The col-6 flag is a diagnostic cross-check against GENCODE membership.

    It is WARNING-only — it must never change the classification, which stays
    authoritative on ref_junctions (GENCODE BED) membership.
    """

    def _write_3col(self, path, rows):
        """rows: list of (junction_id, reads, annotated_flag)."""
        path.write_text("".join(f"{jid}\t{reads}\t{ann}\n" for jid, reads, ann in rows))

    def _write_manifest(self, path, entries):
        lines = ["file_id\tfile_name\tsample_type\tproject_id\n"]
        for fid, stype in entries:
            lines.append(f"{fid}\t{fid}.tsv\t{stype}\tlocal\n")
        path.write_text("".join(lines))

    def _write_reference_bed(self, path, junctions):
        path.write_text("".join(f"{c}\t{s}\t{e}\tjunc\t0\t{st}\n" for c, s, e, st in junctions))

    def test_warns_when_star_annotated_but_not_in_bed(self, tmp_path, caplog):
        # STAR flags chr22:201:300:+ as annotated (1), but it is NOT in our BED
        # → disagreement → WARNING. Classification still uses BED membership, so
        # the junction is unannotated and survives as tumor_exclusive.
        tumor_f = tmp_path / "tumor" / "raw_junctions.tsv"
        tumor_f.parent.mkdir()
        self._write_3col(tumor_f, [
            ("chr22:201:300:+", 100, 1),  # STAR-annotated, absent from BED
            ("chr22:401:500:+", 1, 0),    # noise — below mean, filtered out
        ])

        manifest = tmp_path / "manifest.tsv"
        self._write_manifest(manifest, [("tumor", "Primary Tumor")])
        ref_bed = tmp_path / "ref.bed"
        self._write_reference_bed(ref_bed, [])  # empty → nothing is in our BED

        output = tmp_path / "novel.tsv"
        with caplog.at_level("WARNING"):
            classify_junctions(
                junction_files=[tumor_f],
                manifest_path=manifest,
                reference_bed=ref_bed,
                output_path=output,
            )

        # WARNING fired naming the disagreeing junction.
        assert any(
            "chr22:201:300:+" in rec.message
            for rec in caplog.records
            if rec.levelname == "WARNING"
        )
        # Classification is UNCHANGED — BED membership stays authoritative.
        df = pd.read_csv(output, sep="\t")
        assert len(df) == 1
        assert df.iloc[0]["junction_origin"] == "tumor_exclusive"

    def test_warns_when_in_bed_but_not_star_annotated(self, tmp_path, caplog):
        # STAR flags chr22:101:200:+ as novel (0), but it IS in our BED
        # → disagreement → WARNING. BED membership wins: it is discarded.
        tumor_f = tmp_path / "tumor" / "raw_junctions.tsv"
        tumor_f.parent.mkdir()
        self._write_3col(tumor_f, [
            ("chr22:101:200:+", 100, 0),  # STAR-novel, but present in BED
            ("chr22:401:500:+", 1, 0),    # noise — keeps the other above mean
        ])

        manifest = tmp_path / "manifest.tsv"
        self._write_manifest(manifest, [("tumor", "Primary Tumor")])
        ref_bed = tmp_path / "ref.bed"
        self._write_reference_bed(ref_bed, [("chr22", 100, 200, "+")])

        output = tmp_path / "novel.tsv"
        with caplog.at_level("WARNING"):
            classify_junctions(
                junction_files=[tumor_f],
                manifest_path=manifest,
                reference_bed=ref_bed,
                output_path=output,
            )

        assert any(
            "chr22:101:200:+" in rec.message
            for rec in caplog.records
            if rec.levelname == "WARNING"
        )
        # Annotated per BED → discarded, regardless of STAR's flag.
        df = pd.read_csv(output, sep="\t")
        assert len(df) == 0

    def test_no_warning_when_flags_agree(self, tmp_path, caplog):
        # STAR-annotated AND in BED → agreement; STAR-novel AND not in BED →
        # agreement. No cross-check WARNING should fire.
        tumor_f = tmp_path / "tumor" / "raw_junctions.tsv"
        tumor_f.parent.mkdir()
        self._write_3col(tumor_f, [
            ("chr22:101:200:+", 100, 1),  # annotated, in BED → agree
            ("chr22:201:300:+", 100, 0),  # novel, not in BED → agree
            ("chr22:401:500:+", 1, 0),    # noise
        ])

        manifest = tmp_path / "manifest.tsv"
        self._write_manifest(manifest, [("tumor", "Primary Tumor")])
        ref_bed = tmp_path / "ref.bed"
        self._write_reference_bed(ref_bed, [("chr22", 100, 200, "+")])

        output = tmp_path / "novel.tsv"
        with caplog.at_level("WARNING"):
            classify_junctions(
                junction_files=[tumor_f],
                manifest_path=manifest,
                reference_bed=ref_bed,
                output_path=output,
            )

        assert not any(
            "annotated-flag" in rec.message.lower()
            or "star" in rec.message.lower()
            for rec in caplog.records
            if rec.levelname == "WARNING"
        )

    def test_legacy_two_column_no_crosscheck_warning(self, tmp_path, caplog):
        # 2-column (HISAT2) input → annotated is None → cross-check is skipped,
        # no spurious WARNING.
        tumor_f = tmp_path / "tumor" / "raw_junctions.tsv"
        tumor_f.parent.mkdir()
        tumor_f.write_text("chr22:201:300:+\t100\nchr22:401:500:+\t1\n")

        manifest = tmp_path / "manifest.tsv"
        self._write_manifest(manifest, [("tumor", "Primary Tumor")])
        ref_bed = tmp_path / "ref.bed"
        self._write_reference_bed(ref_bed, [])

        output = tmp_path / "novel.tsv"
        with caplog.at_level("WARNING"):
            classify_junctions(
                junction_files=[tumor_f],
                manifest_path=manifest,
                reference_bed=ref_bed,
                output_path=output,
            )

        assert not any(
            "annotated-flag" in rec.message.lower()
            for rec in caplog.records
            if rec.levelname == "WARNING"
        )


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
        # Issue #940: normal_source is empty for non-normal_shared rows.
        assert pd.isna(df.iloc[0]["normal_source"]) or df.iloc[0]["normal_source"] == ""

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
        # Issue #940: the exclusion carries its provenance.
        assert df.iloc[0]["normal_source"] == "Solid Tissue Normal"

    def test_blood_derived_normal_is_used_and_labeled(self, tmp_path):
        """Issue #940: a Blood Derived Normal still subtracts (conservative
        self-antigen filter), and the resulting normal_shared row is labeled with
        that source so it is distinguishable from a solid-tissue exclusion."""
        tumor_f = tmp_path / "tumor" / "raw_junctions.tsv"
        solid_f = tmp_path / "solid" / "raw_junctions.tsv"
        blood_f = tmp_path / "blood" / "raw_junctions.tsv"
        for f in (tumor_f, solid_f, blood_f):
            f.parent.mkdir()

        # Tumor junctions clearing the mean filter: one shared with blood only,
        # one with solid only, one with BOTH (exercises the sorted comma-join).
        # Plus a noise junction below the mean.
        self._write_junction_file(tumor_f, [
            ("chr22:201:300:+", 100),   # shared with blood
            ("chr22:601:700:+", 100),   # shared with solid
            ("chr22:701:800:+", 100),   # shared with BOTH
            ("chr22:401:500:+", 1),     # noise, filtered out
        ])
        self._write_junction_file(solid_f, [("chr22:601:700:+", 5), ("chr22:701:800:+", 5)])
        self._write_junction_file(blood_f, [("chr22:201:300:+", 5), ("chr22:701:800:+", 5)])

        manifest = tmp_path / "manifest.tsv"
        self._write_manifest(manifest, [
            ("tumor", "Primary Tumor"),
            ("solid", "Solid Tissue Normal"),
            ("blood", "Blood Derived Normal"),
        ])
        ref_bed = tmp_path / "ref.bed"
        self._write_reference_bed(ref_bed, [])

        output = tmp_path / "novel.tsv"
        stats = tmp_path / "stats.tsv"
        classify_junctions(
            junction_files=[tumor_f, solid_f, blood_f],
            manifest_path=manifest,
            reference_bed=ref_bed,
            output_path=output,
            stats_output_path=stats,
        )

        df = pd.read_csv(output, sep="\t")
        src = dict(zip(df["junction_id"], df["normal_source"]))
        assert src["chr22:201:300:+"] == "Blood Derived Normal"
        assert src["chr22:601:700:+"] == "Solid Tissue Normal"
        # a junction in both normals carries both, sorted + comma-joined
        assert src["chr22:701:800:+"] == "Blood Derived Normal,Solid Tissue Normal"
        # all three removed as normal_shared (blood is used, not dropped)
        assert (df["junction_origin"] == "normal_shared").sum() == 3

        # Stats funnel carries the per-source breakdown (a both-normals junction
        # counts under each source, so the per-source rows can sum to > normal_shared),
        # and the funnel normal_shared total counts distinct junctions.
        sdf = pd.read_csv(stats, sep="\t")
        cats = dict(zip(sdf["category"], sdf["count"]))
        assert cats["normal_shared"] == 3
        assert cats["normal_shared:Blood Derived Normal"] == 2   # chr22:201 + chr22:701
        assert cats["normal_shared:Solid Tissue Normal"] == 2    # chr22:601 + chr22:701

    def test_blood_only_manifest_subtracts_and_labels(self, tmp_path):
        """Issue #940 AC: a Blood-Derived-Normal-only manifest (patient_002's shape:
        a CD3+ PBMC as the sole normal) still subtracts, and the exclusion is
        explicitly attributed to blood - not silent, not dropped."""
        tumor_f = tmp_path / "tumor" / "raw_junctions.tsv"
        blood_f = tmp_path / "blood" / "raw_junctions.tsv"
        for f in (tumor_f, blood_f):
            f.parent.mkdir()

        self._write_junction_file(tumor_f, [
            ("chr22:201:300:+", 100),   # shared with blood -> normal_shared (blood)
            ("chr22:601:700:+", 100),   # absent in blood -> tumor_exclusive
            ("chr22:401:500:+", 1),     # noise, filtered out
        ])
        self._write_junction_file(blood_f, [("chr22:201:300:+", 5)])

        manifest = tmp_path / "manifest.tsv"
        self._write_manifest(manifest, [
            ("tumor", "Primary Tumor"),
            ("blood", "Blood Derived Normal"),
        ])
        ref_bed = tmp_path / "ref.bed"
        self._write_reference_bed(ref_bed, [])

        output = tmp_path / "novel.tsv"
        classify_junctions(
            junction_files=[tumor_f, blood_f],
            manifest_path=manifest,
            reference_bed=ref_bed,
            output_path=output,
        )

        df = pd.read_csv(output, sep="\t")
        src = dict(zip(df["junction_id"], df["normal_source"]))
        origin = dict(zip(df["junction_id"], df["junction_origin"]))
        assert origin["chr22:201:300:+"] == "normal_shared"
        assert src["chr22:201:300:+"] == "Blood Derived Normal"
        assert origin["chr22:601:700:+"] == "tumor_exclusive"

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

    def test_gtex_applies_when_no_normal_sample(self, tmp_path):
        # GTEx as the SOLE filter (no matched normal): every other GTEx test passes
        # a normal sample, so the [tumor only, no normal] + GTEx-active path was
        # untested. A no-normal patient is structurally the stacking path with an
        # empty normal set, so GTEx must still partition tumor junctions. One in the
        # blacklist (→ gtex_pantissue_shared) and one in neither (→ tumor_exclusive).
        tumor_f = tmp_path / "tumor" / "raw_junctions.tsv"
        tumor_f.parent.mkdir()
        self._write_junction_file(tumor_f, [
            ("chr22:201:300:+", 100),  # in GTEx → gtex_pantissue_shared
            ("chr22:301:400:+", 100),  # in neither → tumor_exclusive
            ("chr22:901:1000:+", 1),   # noise — below mean, filtered out
        ])

        manifest = tmp_path / "manifest.tsv"
        self._write_manifest(manifest, [("tumor", "Primary Tumor")])
        ref_bed = tmp_path / "ref.bed"
        self._write_reference_bed(ref_bed, [])
        gtex_bed = tmp_path / "gtex.bed"
        self._write_reference_bed(gtex_bed, [("chr22", 200, 300, "+")])

        output = tmp_path / "novel.tsv"
        classify_junctions(
            junction_files=[tumor_f],
            manifest_path=manifest,
            reference_bed=ref_bed,
            output_path=output,
            gtex_bed=gtex_bed,
        )
        df = pd.read_csv(output, sep="\t")
        origins = dict(zip(df["junction_id"], df["junction_origin"]))
        assert origins["chr22:201:300:+"] == "gtex_pantissue_shared"
        assert origins["chr22:301:400:+"] == "tumor_exclusive"
        # GTEx still removed the blacklisted junction with no normal present.
        assert (df["junction_origin"] == "tumor_exclusive").sum() == 1

    def test_gtex_disabled_junction_falls_through_to_tumor_exclusive(self, tmp_path):
        # Contract: when gtex_bed is omitted (filter disabled), a junction that WOULD
        # match a GTEx blacklist on disk is still labeled tumor_exclusive — only the
        # passed argument enables filtering, never the mere existence of the file.
        tumor_f = tmp_path / "tumor" / "raw_junctions.tsv"
        normal_f = tmp_path / "normal" / "raw_junctions.tsv"
        tumor_f.parent.mkdir()
        normal_f.parent.mkdir()
        self._write_junction_file(tumor_f, [
            ("chr22:201:300:+", 100),  # would be gtex_pantissue_shared if the BED were passed
            ("chr22:901:1000:+", 1),   # noise
        ])
        self._write_junction_file(normal_f, [])

        manifest = tmp_path / "manifest.tsv"
        self._write_manifest(manifest, [
            ("tumor", "Primary Tumor"), ("normal", "Solid Tissue Normal"),
        ])
        ref_bed = tmp_path / "ref.bed"
        self._write_reference_bed(ref_bed, [])
        # The BED file exists and contains the junction, but is NOT passed.
        gtex_bed = tmp_path / "gtex.bed"
        self._write_reference_bed(gtex_bed, [("chr22", 200, 300, "+")])

        output = tmp_path / "novel.tsv"
        classify_junctions(
            junction_files=[tumor_f, normal_f],
            manifest_path=manifest,
            reference_bed=ref_bed,
            output_path=output,
            gtex_bed=None,
        )
        df = pd.read_csv(output, sep="\t")
        assert len(df) == 1
        assert df.iloc[0]["junction_origin"] == "tumor_exclusive"

    def _write_3col_junction_file(self, path, rows):
        """rows: list of (junction_id, reads, annotated_flag) — STAR 3-column."""
        path.write_text("".join(f"{jid}\t{reads}\t{ann}\n" for jid, reads, ann in rows))

    def test_star_3col_input_with_gtex_active(self, tmp_path, caplog):
        # Merge-seam regression (Issue #375 × Issue #211/#212): 3-column STAR input
        # AND an active gtex_bed exercised together in one classify_junctions call,
        # the exact interplay reconciled when this branch was merged onto main. Proves
        # (a) the annotated flag flows through the GTEx classification path, (b) the
        # WARNING-only cross-check counters stay independent of the gtex counter, and
        # (c) a cross-check disagreement fires alongside gtex labeling without
        # perturbing either origin assignment.
        tumor_f = tmp_path / "tumor" / "raw_junctions.tsv"
        normal_f = tmp_path / "normal" / "raw_junctions.tsv"
        tumor_f.parent.mkdir()
        normal_f.parent.mkdir()
        self._write_3col_junction_file(tumor_f, [
            ("chr22:101:200:+", 100, 1),  # flag=1, in ref BED → agree → annotated, discarded
            ("chr22:201:300:+", 100, 0),  # flag=0, in GTEx → gtex_pantissue_shared (agree)
            ("chr22:301:400:+", 100, 1),  # flag=1 but NOT in ref → DISAGREE → warn; tumor_exclusive
            ("chr22:901:1000:+", 1, 0),   # noise — below mean, filtered out
        ])
        self._write_junction_file(normal_f, [])

        manifest = tmp_path / "manifest.tsv"
        self._write_manifest(manifest, [
            ("tumor", "Primary Tumor"), ("normal", "Solid Tissue Normal"),
        ])
        ref_bed = tmp_path / "ref.bed"
        self._write_reference_bed(ref_bed, [("chr22", 100, 200, "+")])
        gtex_bed = tmp_path / "gtex.bed"
        self._write_reference_bed(gtex_bed, [("chr22", 200, 300, "+")])

        output = tmp_path / "novel.tsv"
        with caplog.at_level("WARNING"):
            classify_junctions(
                junction_files=[tumor_f, normal_f],
                manifest_path=manifest,
                reference_bed=ref_bed,
                output_path=output,
                gtex_bed=gtex_bed,
            )

        df = pd.read_csv(output, sep="\t")
        origins = dict(zip(df["junction_id"], df["junction_origin"]))
        # Annotated junction discarded (BED membership authoritative), not emitted.
        assert "chr22:101:200:+" not in origins
        # GTEx and tumor_exclusive partition the two surviving unannotated junctions.
        assert origins["chr22:201:300:+"] == "gtex_pantissue_shared"
        assert origins["chr22:301:400:+"] == "tumor_exclusive"
        assert (df["junction_origin"] == "gtex_pantissue_shared").sum() == 1
        assert (df["junction_origin"] == "tumor_exclusive").sum() == 1
        # The cross-check WARNING fired for the one disagreeing junction, and the
        # GTEx-labeled junction (flag=0, novel, agreeing) did NOT trigger one.
        warnings = [r.message for r in caplog.records if r.levelname == "WARNING"]
        assert any("chr22:301:400:+" in m for m in warnings)
        assert not any("chr22:201:300:+" in m and "annotated-flag" in m.lower()
                       for m in warnings)


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
        # 4 descriptive distribution rows (Issue #215) + 1 per-normal-source
        # breakdown row (Issue #940, one normal type removed a junction) = 11 total.
        assert len(stats) == 11
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
        # Issue #940: per-source breakdown row for the one normal type used.
        assert by_cat["normal_shared:Solid Tissue Normal"] == 1

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
