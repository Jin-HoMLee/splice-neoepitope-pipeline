# =============================================================================
# Rule module: FASTQ alignment — STAR or HISAT2
# =============================================================================
#
# Provides a FASTQ-based input pathway as an alternative to GDC downloads.
# The aligner is selected via config.alignment.aligner:
#
#   "star"   — Full accuracy, ~32 GB RAM (default; same aligner used by GDC)
#   "hisat2" — Lower memory (~8 GB), suitable for laptops / small servers
#
# Both aligners produce junction quantification TSVs in the same format so
# all downstream rules are aligner-agnostic.
#
# To use:
#   1. Set `data_source: "fastq"` in config/config.yaml
#   2. Set `alignment.aligner: "star"` or `"hisat2"`
#   3. Populate config/samples.tsv with patient_id, sample_id, fastq paths
#   4. Run the pipeline
#
# =============================================================================

import csv
import os
from pathlib import Path


# Only define these rules when running in fastq mode
if config.get("data_source") == "fastq":

    # ── Shared TSV reader ────────────────────────────────────────────────────

    def _read_samples_tsv(samples_tsv, patient_id=None):
        """Read rows from the samples TSV, skipping blank/comment lines.

        Args:
            samples_tsv: path to the TSV (string or Path)
            patient_id:  if given, return only rows matching this patient

        Returns a list of dicts (one per sample row).
        """
        samples_path = Path(samples_tsv)
        if not samples_path.exists():
            return []
        rows = []
        with samples_path.open() as f:
            for row in csv.DictReader(f, delimiter="\t"):
                pid = (row.get("patient_id") or "").strip()
                if not pid or pid.startswith("#"):
                    continue
                row["patient_id"] = pid
                rows.append(row)
        if patient_id is not None:
            rows = [r for r in rows if r["patient_id"] == patient_id]
        return rows


    # ── Shared manifest + checkpoint ─────────────────────────────────────────

    rule create_alignment_manifest:
        """Create a manifest TSV from the samples configuration.

        Mirrors the GDC manifest format so downstream rules use the same
        filtering logic regardless of data source.
        """
        output:
            manifest=os.path.join(OUT["raw_data"], "{patient_id}", "manifest.tsv"),
        log:
            os.path.join(OUT["logs"], "alignment", "{patient_id}_manifest.log"),
        params:
            samples_tsv=config["samples_tsv"],
        conda:
            "../envs/python.yaml"
        run:
            output_path = Path(output.manifest)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            with output_path.open("w") as fout:
                fout.write("file_id\tfile_name\tsample_type\tproject_id\n")
                for row in _read_samples_tsv(params.samples_tsv, patient_id=wildcards.patient_id):
                    sid = row["sample_id"]
                    sample_type = row.get("sample_type", "Unknown")
                    fout.write(f"{sid}\t{sid}.tsv\t{sample_type}\t{wildcards.patient_id}\n")


    checkpoint alignment_complete:
        """Checkpoint that triggers after all samples are aligned.

        Replaces the GDC download checkpoint for fastq mode. Both STAR and
        HISAT2 produce junction TSVs at the same output path so this
        checkpoint is aligner-agnostic.
        """
        input:
            manifest=os.path.join(OUT["raw_data"], "{patient_id}", "manifest.tsv"),
            samples=lambda wildcards: expand(
                os.path.join(OUT["raw_data"], wildcards.patient_id, "files", "{sample}.tsv"),
                sample=[s["sample_id"] for s in _read_samples_tsv(config["samples_tsv"], wildcards.patient_id)],
            ),
        output:
            done=os.path.join(OUT["raw_data"], "{patient_id}", "download.done"),
            data_dir=directory(os.path.join(OUT["raw_data"], "{patient_id}", "files")),
        shell:
            "touch {output.done}"


    # ── Shared FASTQ input functions ─────────────────────────────────────────
    # These are aligner-agnostic: they look up fastq paths from the TSV only.

    def _get_fastq1(wildcards):
        for s in _read_samples_tsv(config["samples_tsv"], wildcards.patient_id):
            if s["sample_id"] == wildcards.sample:
                return s["fastq1"]
        raise ValueError(f"Sample not found: {wildcards.sample} (patient: {wildcards.patient_id})")


    def _get_fastq2(wildcards):
        for s in _read_samples_tsv(config["samples_tsv"], wildcards.patient_id):
            if s["sample_id"] == wildcards.sample:
                fastq2 = s.get("fastq2", "")
                return [fastq2] if fastq2 else []
        return []


    # Shared output paths for both align rules — defined once to avoid drift.
    _JUNCTION_OUTPUT = os.path.join(OUT["raw_data"], "{patient_id}", "files", "{sample}.tsv")
    _JUNCTION_DONE   = os.path.join(OUT["raw_data"], "{patient_id}", "files", "{sample}.done")

    # ── HISAT2 ───────────────────────────────────────────────────────────────

    if config.get("alignment", {}).get("aligner") == "hisat2":

        # Configurable index directory — allows test (chr22) and production
        # (full genome) to maintain separate indices without overwriting.
        _HISAT2_INDEX_DIR = config.get("alignment", {}).get(
            "hisat2_index_dir", "resources/hisat2_index"
        )

        rule hisat2_index:
            """Build HISAT2 genome index (one-time; reused for all samples).

            Requires ~8 GB RAM for the full human genome, or ~1 GB for a
            single-chromosome test reference (e.g. chr22).
            """
            input:
                genome=config["reference"]["genome_fasta"],
            output:
                index_dir=directory(_HISAT2_INDEX_DIR),
                done=touch(os.path.join(_HISAT2_INDEX_DIR, "index.done")),
            log:
                os.path.join(OUT["logs"], "hisat2", "hisat2_index.log"),
            threads: 8
            resources:
                mem_mb=8000,
            conda:
                "../envs/hisat2.yaml"
            params:
                index_prefix=os.path.join(_HISAT2_INDEX_DIR, "genome"),
            shell:
                """
                mkdir -p {output.index_dir}
                hisat2-build \\
                    -p {threads} \\
                    {input.genome} \\
                    {params.index_prefix} \\
                    2>&1 | tee {log}
                """


        rule hisat2_align:
            """Run HISAT2 alignment on a single sample.

            Input:  FASTQ files (single-end or paired-end)
            Output: Junction quantification TSV in pipeline-compatible format

            Uses regtools to extract junctions from the BAM file.
            """
            input:
                index_dir=_HISAT2_INDEX_DIR,
                index_done=os.path.join(_HISAT2_INDEX_DIR, "index.done"),
                fastq1=_get_fastq1,
                fastq2=_get_fastq2,
            output:
                junctions=_JUNCTION_OUTPUT,
                done=touch(_JUNCTION_DONE),
            log:
                os.path.join(OUT["logs"], "hisat2", "{patient_id}_{sample}_align.log"),
            params:
                index_prefix=os.path.join(_HISAT2_INDEX_DIR, "genome"),
                output_prefix=lambda wildcards: os.path.join(
                    OUT["raw_data"], wildcards.patient_id, "files", wildcards.sample
                ),
            threads: 8
            resources:
                mem_mb=8000,
            conda:
                "../envs/hisat2.yaml"
            shell:
                """
                mkdir -p $(dirname {output.junctions})

                if [[ -n "{input.fastq2}" ]]; then
                    FASTQ_ARGS="-1 {input.fastq1} -2 {input.fastq2}"
                else
                    FASTQ_ARGS="-U {input.fastq1}"
                fi

                hisat2 \\
                    -p {threads} \\
                    -x {params.index_prefix} \\
                    $FASTQ_ARGS \\
                    2>> {log} | \\
                    samtools sort -@ {threads} -o {params.output_prefix}.bam - 2>> {log}

                samtools index {params.output_prefix}.bam 2>> {log}

                regtools junctions extract \\
                    -s XS \\
                    -a 8 \\
                    -m 50 \\
                    -M 500000 \\
                    -o {params.output_prefix}_junctions.bed \\
                    {params.output_prefix}.bam \\
                    2>> {log}

                awk -F'\\t' '{{
                    if ($5 > 0) print $1":"$2":"$3":"$6"\\t"$5
                }}' {params.output_prefix}_junctions.bed > {output.junctions}

                rm -f {params.output_prefix}.bam {params.output_prefix}.bam.bai {params.output_prefix}_junctions.bed

                echo "Extracted $(wc -l < {output.junctions}) junctions from HISAT2 output" >> {log}
                """


    # ── STAR ─────────────────────────────────────────────────────────────────

    else:  # default aligner: star

        _STAR_INDEX_DIR = config.get("alignment", {}).get(
            "star_index_dir", "resources/star_index"
        )

        rule star_index:
            """Build STAR genome index (one-time; reused for all samples).

            Requires ~32 GB RAM for the full human genome.
            """
            input:
                genome=config["reference"]["genome_fasta"],
                gtf=config["reference"]["gencode_gtf"],
            output:
                index_dir=directory(_STAR_INDEX_DIR),
                done=touch(os.path.join(_STAR_INDEX_DIR, "index.done")),
            log:
                os.path.join(OUT["logs"], "star", "star_index.log"),
            threads: 8
            resources:
                mem_mb=32000,
            conda:
                "../envs/star.yaml"
            shell:
                """
                GTF_FILE="{input.gtf}"
                if [[ "$GTF_FILE" == *.gz ]]; then
                    gunzip -c "$GTF_FILE" > resources/temp_annotation.gtf
                    GTF_FILE="resources/temp_annotation.gtf"
                fi

                STAR \\
                    --runMode genomeGenerate \\
                    --runThreadN {threads} \\
                    --genomeDir {output.index_dir} \\
                    --genomeFastaFiles {input.genome} \\
                    --sjdbGTFfile "$GTF_FILE" \\
                    --sjdbOverhang 100 \\
                    2>&1 | tee {log}

                if [[ -f "resources/temp_annotation.gtf" ]]; then
                    rm resources/temp_annotation.gtf
                fi
                """


        rule star_align:
            """Run STAR alignment on a single sample.

            Input:  FASTQ files (single-end or paired-end)
            Output: Junction quantification TSV in GDC-compatible format
            """
            input:
                index_dir=_STAR_INDEX_DIR,
                index_done=os.path.join(_STAR_INDEX_DIR, "index.done"),
                fastq1=_get_fastq1,
                fastq2=_get_fastq2,
            output:
                junctions=_JUNCTION_OUTPUT,
                done=touch(_JUNCTION_DONE),
            log:
                os.path.join(OUT["logs"], "star", "{patient_id}_{sample}_align.log"),
            params:
                output_prefix=lambda wildcards: os.path.join(
                    OUT["raw_data"], wildcards.patient_id, "files", f"{wildcards.sample}_"
                ),
            threads: 8
            resources:
                mem_mb=32000,
            conda:
                "../envs/star.yaml"
            shell:
                """
                READCMD=""
                if [[ "{input.fastq1}" == *.gz ]]; then
                    READCMD="--readFilesCommand zcat"
                fi

                FASTQ_FILES="{input.fastq1}"
                if [[ -n "{input.fastq2}" ]]; then
                    FASTQ_FILES="{input.fastq1} {input.fastq2}"
                fi

                STAR \\
                    --runMode alignReads \\
                    --runThreadN {threads} \\
                    --genomeDir {input.index_dir} \\
                    --readFilesIn $FASTQ_FILES \\
                    $READCMD \\
                    --outFileNamePrefix {params.output_prefix} \\
                    --outSAMtype None \\
                    --outSJfilterReads Unique \\
                    --outSJfilterCountUniqueMin 1 1 1 1 \\
                    --outSJfilterCountTotalMin 1 1 1 1 \\
                    2>&1 | tee {log}

                awk -F'\\t' '{{
                    strand = ".";
                    if ($4 == 1) strand = "+";
                    else if ($4 == 2) strand = "-";
                    if ($7 > 0) print $1":"$2":"$3":"strand"\\t"$7
                }}' {params.output_prefix}SJ.out.tab > {output.junctions}

                echo "Converted $(wc -l < {output.junctions}) junctions from STAR output"
                """
