# =============================================================================
# Rule module: FASTQ alignment — STAR or HISAT2
# =============================================================================
#
# The aligner is selected via config.alignment.aligner:
#
#   "star"   — Full accuracy, ~32 GB RAM
#   "hisat2" — Lower memory (~8 GB), suitable for laptops / small servers
#
# Both aligners produce junction quantification TSVs in the same format so
# all downstream rules are aligner-agnostic.
#
# To use:
#   1. Set `alignment.aligner: "star"` or `"hisat2"` in config/config.yaml
#   2. Populate config/samples.tsv with patient_id, sample_id, fastq paths
#   3. Run the pipeline
#
# =============================================================================

# Helpers (_read_samples_tsv, _local_fastq) are defined in common.smk,
# which is included before this file by the Snakefile.

# ── Shared manifest + checkpoint ─────────────────────────────────────────────

rule create_alignment_manifest:
    """Create a manifest TSV from the samples configuration."""
    output:
        manifest=os.path.join(_RES, "{patient_id}", "alignment", "manifest.tsv"),
    log:
        os.path.join(_LOGS, "{patient_id}", "alignment", "manifest.log"),
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


# ── Shared FASTQ input functions ─────────────────────────────────────────────
# Aligner-agnostic: look up fastq paths from the TSV only.

def _get_fastq1(wildcards):
    for s in _read_samples_tsv(config["samples_tsv"], wildcards.patient_id):
        if s["sample_id"] == wildcards.sample:
            return _local_fastq(s["fastq1"], wildcards.patient_id, wildcards.sample)
    raise ValueError(f"Sample not found: {wildcards.sample} (patient: {wildcards.patient_id})")


def _get_fastq2(wildcards):
    for s in _read_samples_tsv(config["samples_tsv"], wildcards.patient_id):
        if s["sample_id"] == wildcards.sample:
            fastq2 = (s.get("fastq2") or "").strip()
            return [_local_fastq(fastq2, wildcards.patient_id, wildcards.sample)] if fastq2 else []
    return []


# Shared output paths for both align rules — defined once to avoid drift.
_JUNCTION_OUTPUT = os.path.join(_RES, "{patient_id}", "alignment", "{sample}", "raw_junctions.tsv")
_JUNCTION_DONE   = os.path.join(_RES, "{patient_id}", "alignment", "{sample}", "done")

# ── HISAT2 ───────────────────────────────────────────────────────────────────

if config.get("alignment", {}).get("aligner") == "hisat2":

    # Per-sample strandness support — translates samples.tsv `strandness` column
    # (biological direction: `unstranded`/`forward`/`reverse`) into the HISAT2
    # --rna-strandness flag value (`F`/`R`/`FR`/`RF` or empty). Pure helpers live
    # in workflow/scripts/strandness.py and are unit-tested in test_strandness.py.
    # `srcdir()` is unavailable in Snakemake 8 — use `workflow.basedir` (which
    # points at the Snakefile's directory, i.e. the repo root here) instead.
    import sys
    sys.path.insert(0, os.path.join(workflow.basedir, "workflow", "scripts"))
    from strandness import get_strandness_from_row
    from hisat2_command import build_read_args

    def _get_hisat2_read_args(wildcards, input):
        """Resolve the HISAT2 read-input args (`-U` vs `-1/-2`) for one sample.

        Mirrors `_get_hisat2_strandness`: delegates the single-end / paired-end
        selection to the pure `build_read_args` helper (unit-tested in
        test_hisat2_command.py). Reads the already-resolved input paths so the
        emitted command uses exactly the staged files. `input.fastq2` is a
        (possibly empty) list from `_get_fastq2`.
        """
        fastq2 = input.fastq2[0] if input.fastq2 else ""
        return build_read_args(input.fastq1, fastq2)

    def _get_hisat2_strandness(wildcards):
        """Resolve the HISAT2 --rna-strandness flag value for a single sample.

        Returns the empty string for samples without a `strandness` column entry
        or with `unstranded` — the rule then omits the flag entirely (preserves
        backward compat with samples.tsv files predating this column).
        """
        for s in _read_samples_tsv(config["samples_tsv"], wildcards.patient_id):
            if s["sample_id"] == wildcards.sample:
                return get_strandness_from_row(s)
        return ""

    # Configurable index directory — allows test (chr22) and production
    # (full genome) to maintain separate indices without overwriting.
    _HISAT2_INDEX_DIR = config.get("alignment", {}).get(
        "hisat2_index_dir", "indices/hisat2"
    )
    _HISAT2_PREBUILT_URL = config.get("alignment", {}).get("hisat2_prebuilt_url", "")

    if _HISAT2_PREBUILT_URL:
        # genome_tran index: hg38 (UCSC naming, chr-prefix) + GENCODE splice sites baked in.
        # Files unpack as genome_tran.N.ht2 directly into the index dir.
        _HISAT2_INDEX_PREFIX = os.path.join(_HISAT2_INDEX_DIR, "genome_tran")

        rule hisat2_download_index:
            """Download pre-built HISAT2 hg38 index (UCSC naming, genome_tran, with splice sites)."""
            output:
                index_dir=directory(_HISAT2_INDEX_DIR),
                done=touch(os.path.join(_HISAT2_INDEX_DIR, "index.done")),
            log:
                os.path.join(_SHARED_LOG, "alignment", "hisat2_index.log"),
            params:
                url=_HISAT2_PREBUILT_URL,
            resources:
                mem_mb=1000,
            shell:
                """
                set -euo pipefail
                mkdir -p {output.index_dir}
                ( curl --fail -L --no-progress-meter \
                    --retry 3 --retry-connrefused --retry-delay 5 \
                    {params.url} \
                    | tar -xz --strip-components=1 -C {output.index_dir} ) \
                    2>&1 | tee {log}
                """

    else:
        _HISAT2_INDEX_PREFIX = os.path.join(_HISAT2_INDEX_DIR, "genome")

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
                os.path.join(_SHARED_LOG, "alignment", "hisat2_index.log"),
            threads: config.get("alignment", {}).get("threads", 8)
            resources:
                mem_mb=8000,
            conda:
                "../envs/hisat2.yaml"
            params:
                index_prefix=_HISAT2_INDEX_PREFIX,
            shell:
                """
                set -euo pipefail
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
        Output: Junction quantification TSV, BAM, BAI, and junction BED

        Uses regtools to extract junctions from the BAM file.
        BAM/BAI/BED are retained on disk; the run_cloud_gpu.sh script
        uploads the entire results directory to GCS after the pipeline
        completes.
        """
        input:
            index_dir=_HISAT2_INDEX_DIR,
            index_done=os.path.join(_HISAT2_INDEX_DIR, "index.done"),
            fastq1=ancient(_get_fastq1),
            fastq2=ancient(_get_fastq2),
        output:
            junctions=_JUNCTION_OUTPUT,
            done=touch(_JUNCTION_DONE),
            bam=os.path.join(_RES, "{patient_id}", "alignment", "{sample}", "{sample}.bam"),
            bai=os.path.join(_RES, "{patient_id}", "alignment", "{sample}", "{sample}.bam.bai"),
            bed=os.path.join(_RES, "{patient_id}", "alignment", "{sample}", "{sample}_junctions.bed"),
        log:
            os.path.join(_LOGS, "{patient_id}", "alignment", "{sample}_align.log"),
        params:
            index_prefix=_HISAT2_INDEX_PREFIX,
            strandness=_get_hisat2_strandness,
            read_args=_get_hisat2_read_args,
        threads: config.get("alignment", {}).get("threads", 8)
        resources:
            mem_mb=20000,
        conda:
            "../envs/hisat2.yaml"
        shell:
            """
            set -euo pipefail
            mkdir -p $(dirname {output.junctions})

            # Fail fast if the index is missing — avoids samtools hanging on a
            # broken pipe, which would prevent Snakemake from writing the error
            # to pipeline.log and leave the orchestrator polling indefinitely.
            if [[ ! -f "{params.index_prefix}.1.ht2" ]]; then
                echo "ERROR: HISAT2 index not found at {params.index_prefix}.*.ht2" | tee -a {log}
                exit 1
            fi

            STRANDNESS_ARGS=""
            if [[ -n "{params.strandness}" ]]; then
                STRANDNESS_ARGS="--rna-strandness {params.strandness}"
            fi

            hisat2 \\
                -p {threads} \\
                -x {params.index_prefix} \\
                $STRANDNESS_ARGS \\
                {params.read_args} \\
                2>> {log} | \\
                samtools sort -@ {threads} -m 1G -o {output.bam} - 2>> {log}

            samtools index {output.bam} 2>> {log}

            regtools junctions extract \\
                -s XS \\
                -a 8 \\
                -m 50 \\
                -M 500000 \\
                -o {output.bed} \\
                {output.bam} \\
                2>> {log}

            # regtools BED12 cols 2-3 are anchor outer boundaries, NOT intron
            # donor/acceptor — see Issue #370. The helper derives the real
            # intron coords from blockSizes/blockStarts.
            python workflow/scripts/bed12_to_junctions.py \\
                --input {output.bed} \\
                --output {output.junctions} \\
                2>> {log}

            echo "Extracted $(wc -l < {output.junctions}) junctions from HISAT2 output" >> {log}
            """


# ── STAR ─────────────────────────────────────────────────────────────────────

elif config.get("alignment", {}).get("aligner") == "star":

    _STAR_INDEX_DIR = config.get("alignment", {}).get(
        "star_index_dir", "indices/star"
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
            os.path.join(_SHARED_LOG, "alignment", "star_index.log"),
        threads: config.get("alignment", {}).get("threads", 8)
        resources:
            mem_mb=32000,
        conda:
            "../envs/star.yaml"
        shell:
            """
            set -euo pipefail
            GTF_FILE="{input.gtf}"
            if [[ "$GTF_FILE" == *.gz ]]; then
                TEMP_GTF=$(mktemp)
                trap 'rm -f "$TEMP_GTF"' EXIT
                gunzip -c "$GTF_FILE" > "$TEMP_GTF"
                GTF_FILE="$TEMP_GTF"
            fi

            STAR \\
                --runMode genomeGenerate \\
                --runThreadN {threads} \\
                --genomeDir {output.index_dir} \\
                --genomeFastaFiles {input.genome} \\
                --sjdbGTFfile "$GTF_FILE" \\
                --sjdbOverhang 100 \\
                2>&1 | tee {log}
            """


    rule star_align:
        """Run STAR alignment on a single sample.

        Input:  FASTQ files (single-end or paired-end)
        Output: Junction quantification TSV in pipeline-compatible format
        """
        input:
            index_dir=_STAR_INDEX_DIR,
            index_done=os.path.join(_STAR_INDEX_DIR, "index.done"),
            fastq1=ancient(_get_fastq1),
            fastq2=ancient(_get_fastq2),
        output:
            junctions=_JUNCTION_OUTPUT,
            done=touch(_JUNCTION_DONE),
        log:
            os.path.join(_LOGS, "{patient_id}", "alignment", "{sample}_align.log"),
        params:
            output_prefix=lambda wildcards: os.path.join(
                _RES, wildcards.patient_id, "alignment", wildcards.sample, "star_"
            ),
        threads: config.get("alignment", {}).get("threads", 8)
        resources:
            mem_mb=32000,
        conda:
            "../envs/star.yaml"
        shell:
            """
            set -euo pipefail
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
                --twopassMode Basic \\
                --limitSjdbInsertNsj 2000000 \\
                2>&1 | tee {log}

            # SJ.out.tab col 4=0 means STAR couldn't infer strand; rescue from
            # col 5 (intron motif) where possible, drop truly non-canonical —
            # see Issue #374. The helper also handles the col 7 > 0 filter
            # and emits the same raw_junctions.tsv format as the HISAT2 path.
            python workflow/scripts/star_sj_to_junctions.py \\
                --input {params.output_prefix}SJ.out.tab \\
                --output {output.junctions} \\
                2>&1 | tee -a {log}
            """

else:
    _aligner = config.get("alignment", {}).get("aligner", "<unset>")
    raise ValueError(
        f"Unknown aligner '{_aligner}'. "
        "Set config.alignment.aligner to 'star' or 'hisat2'."
    )
