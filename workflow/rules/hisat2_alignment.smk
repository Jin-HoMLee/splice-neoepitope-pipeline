# =============================================================================
# Rule module: HISAT2 alignment — Low-memory alternative to STAR
# =============================================================================
#
# This module provides a lightweight alternative to STAR for users with limited
# RAM (e.g., laptops or small servers). HISAT2 requires only ~8 GB RAM compared
# to STAR's 32 GB, making it accessible to more users.
#
# Trade-offs vs STAR:
#   - Lower memory (~8 GB vs ~32 GB)
#   - Smaller index (~8 GB vs ~30 GB)
#   - Slightly lower sensitivity for novel junctions
#   - Still produces compatible junction quantification output
#
# To use HISAT2 instead of STAR:
#   1. Set `data_source: "fastq"` in config/config.yaml
#   2. Set `alignment.aligner: "hisat2"` in config/config.yaml
#   3. Run the pipeline
#
# =============================================================================

import os
from pathlib import Path


# Only define these rules when running in fastq mode with hisat2
if config.get("data_source") == "fastq" and config.get("alignment", {}).get("aligner") == "hisat2":

    # Configurable index directory — allows test (chr22) and production (full genome)
    # to maintain separate indices without overwriting each other.
    # Override via config: alignment.hisat2_index_dir: "resources/test/hisat2_index"
    _HISAT2_INDEX_DIR = config.get("alignment", {}).get("hisat2_index_dir", "resources/hisat2_index")

    rule hisat2_index:
        """Build HISAT2 genome index for local alignment.

        This rule only runs once; the index is reused for all samples.
        HISAT2 index generation requires ~8 GB RAM for the full human genome,
        or ~1 GB for a single-chromosome test reference (e.g. chr22).
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


    def get_local_samples_hisat2(patient_id=None):
        """Read samples from the local samples TSV, optionally filtered by patient_id."""
        import csv
        samples_file = Path(config["samples_tsv"])
        if not samples_file.exists():
            return []
        with samples_file.open() as f:
            rows = [r for r in csv.DictReader(f, delimiter="\t")
                    if not r["patient_id"].startswith("#")]
        if patient_id is not None:
            rows = [r for r in rows if r["patient_id"] == patient_id]
        return rows


    def get_fastq1_hisat2(wildcards):
        """Get path to read 1 FASTQ for a sample."""
        for sample in get_local_samples_hisat2(wildcards.patient_id):
            if sample["sample_id"] == wildcards.sample:
                return sample["fastq1"]
        raise ValueError(f"Sample not found: {wildcards.sample} (patient: {wildcards.patient_id})")


    def get_fastq2_hisat2(wildcards):
        """Get path to read 2 FASTQ for a sample (or empty list if single-end)."""
        for sample in get_local_samples_hisat2(wildcards.patient_id):
            if sample["sample_id"] == wildcards.sample:
                fastq2 = sample.get("fastq2", "")
                return [fastq2] if fastq2 else []
        return []


    rule hisat2_align:
        """Run HISAT2 alignment on a single sample to generate junction quantification.
        
        Input: FASTQ files (single-end or paired-end)
        Output: Junction quantification TSV in pipeline-compatible format
        
        Uses regtools to extract junctions from the BAM file.
        """
        input:
            index_dir=_HISAT2_INDEX_DIR,
            index_done=os.path.join(_HISAT2_INDEX_DIR, "index.done"),
            fastq1=get_fastq1_hisat2,
            fastq2=get_fastq2_hisat2,
        output:
            junctions=os.path.join(OUT["raw_data"], "{patient_id}", "files", "{sample}.tsv"),
            done=touch(os.path.join(OUT["raw_data"], "{patient_id}", "files", "{sample}.done")),
        log:
            os.path.join(OUT["logs"], "hisat2", "{patient_id}_{sample}_align.log"),
        params:
            index_prefix=os.path.join(_HISAT2_INDEX_DIR, "genome"),
            output_prefix=lambda wildcards: os.path.join(OUT["raw_data"], wildcards.patient_id, "files", wildcards.sample),
        threads: 8
        resources:
            mem_mb=8000,  # HISAT2 alignment needs only ~8 GB
        conda:
            "../envs/hisat2.yaml"
        shell:
            """
            # Create output directory
            mkdir -p $(dirname {output.junctions})

            # Determine if paired-end
            if [[ -n "{input.fastq2}" ]]; then
                FASTQ_ARGS="-1 {input.fastq1} -2 {input.fastq2}"
            else
                FASTQ_ARGS="-U {input.fastq1}"
            fi

            # Run HISAT2 alignment and pipe to BAM
            hisat2 \\
                -p {threads} \\
                -x {params.index_prefix} \\
                $FASTQ_ARGS \\
                2>> {log} | \\
                samtools sort -@ {threads} -o {params.output_prefix}.bam - 2>> {log}

            # Index the BAM
            samtools index {params.output_prefix}.bam 2>> {log}
            
            # Extract junctions using regtools
            regtools junctions extract \\
                -s XS \\
                -a 8 \\
                -m 50 \\
                -M 500000 \\
                -o {params.output_prefix}_junctions.bed \\
                {params.output_prefix}.bam \\
                2>> {log}
            
            # Convert regtools BED to pipeline format
            # regtools BED: chrom, start, end, name, score(reads), strand
            # Pipeline format: junction_id(chr:start:end:strand) \\t mapped_reads
            awk -F'\\t' '{{
                if ($5 > 0) print $1":"$2":"$3":"$6"\\t"$5
            }}' {params.output_prefix}_junctions.bed > {output.junctions}
            
            # Clean up intermediate files
            rm -f {params.output_prefix}.bam {params.output_prefix}.bam.bai {params.output_prefix}_junctions.bed
            
            echo "Extracted $(wc -l < {output.junctions}) junctions from HISAT2 output" >> {log}
            """


    rule create_hisat2_manifest:
        """Create a manifest TSV from the local samples configuration.

        This mirrors the GDC manifest format so downstream rules can use the
        same filtering logic regardless of data source.
        """
        output:
            manifest=os.path.join(OUT["raw_data"], "{patient_id}", "manifest.tsv"),
        log:
            os.path.join(OUT["logs"], "hisat2", "{patient_id}_manifest.log"),
        params:
            samples_tsv=config["samples_tsv"],
        conda:
            "../envs/python.yaml"
        run:
            import csv
            from pathlib import Path

            samples_file = Path(params.samples_tsv)
            output_path = Path(output.manifest)
            output_path.parent.mkdir(parents=True, exist_ok=True)

            with samples_file.open() as fin, output_path.open("w") as fout:
                reader = csv.DictReader(fin, delimiter="\t")
                fout.write("file_id\tfile_name\tsample_type\tproject_id\n")
                for row in reader:
                    if row["patient_id"] != wildcards.patient_id:
                        continue
                    sample_id = row["sample_id"]
                    sample_type = row.get("sample_type", "Unknown")
                    fout.write(f"{sample_id}\t{sample_id}.tsv\t{sample_type}\t{wildcards.patient_id}\n")


    # Checkpoint to collect all aligned samples
    checkpoint hisat2_alignment_complete:
        """Checkpoint that triggers after all local samples are aligned.

        This replaces the GDC download checkpoint for local mode.
        """
        input:
            manifest=os.path.join(OUT["raw_data"], "{patient_id}", "manifest.tsv"),
            samples=lambda wildcards: expand(
                os.path.join(OUT["raw_data"], wildcards.patient_id, "files", "{sample}.tsv"),
                sample=[s["sample_id"] for s in get_local_samples_hisat2(wildcards.patient_id)]
            ),
        output:
            done=os.path.join(OUT["raw_data"], "{patient_id}", "download.done"),
            data_dir=directory(os.path.join(OUT["raw_data"], "{patient_id}", "files")),
        shell:
            "touch {output.done}"
