# =============================================================================
# Rule module: Local STAR alignment — Alternative to GDC download
# =============================================================================
#
# This module provides an alternative data input pathway for users who cannot
# access controlled TCGA data from the GDC portal. Instead of downloading
# pre-computed splice junction quantification files, users can align their own
# RNA-Seq FASTQ files locally using STAR.
#
# STAR is the same aligner used by the GDC to process TCGA data, ensuring
# compatibility with downstream analysis steps.
#
# To use local alignment mode with STAR:
#   1. Set `data_source: "local"` in config/config.yaml
#   2. Set `alignment.aligner: "star"` in config/config.yaml
#   3. Provide a samples TSV file listing FASTQ paths
#   4. Run the pipeline
#
# Note: STAR requires ~32 GB RAM. For lower-memory systems, use HISAT2 instead.
#
# =============================================================================

import os


# Only define these rules when running in local mode with star
if config.get("data_source") == "local" and config.get("alignment", {}).get("aligner") == "star":

    rule star_index:
        """Build STAR genome index for local alignment.
        
        This rule only runs once; the index is reused for all samples.
        Note: STAR index generation requires ~32 GB RAM for the human genome.
        """
        input:
            genome=config["reference"]["genome_fasta"],
            gtf=config["reference"]["gencode_gtf"],
        output:
            index_dir=directory("resources/star_index"),
            done=touch("resources/star_index/index.done"),
        log:
            os.path.join(OUT["logs"], "star", "star_index.log"),
        threads: 8
        resources:
            mem_mb=32000,  # STAR needs ~32 GB for human genome indexing
        conda:
            "../envs/star.yaml"
        shell:
            """
            # Handle gzipped GTF
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
            
            # Clean up temp file if created
            if [[ -f "resources/temp_annotation.gtf" ]]; then
                rm resources/temp_annotation.gtf
            fi
            """


    def get_local_samples():
        """Read samples from the local samples TSV file."""
        import csv
        samples_file = Path(config["samples_tsv"])
        if not samples_file.exists():
            return []
        with samples_file.open() as f:
            reader = csv.DictReader(f, delimiter="\t")
            return list(reader)


    def get_fastq1(wildcards):
        """Get path to read 1 FASTQ for a sample."""
        for sample in get_local_samples():
            if sample["sample_id"] == wildcards.sample:
                return sample["fastq1"]
        raise ValueError(f"Sample not found: {wildcards.sample}")


    def get_fastq2(wildcards):
        """Get path to read 2 FASTQ for a sample (or empty list if single-end)."""
        for sample in get_local_samples():
            if sample["sample_id"] == wildcards.sample:
                fastq2 = sample.get("fastq2", "")
                return [fastq2] if fastq2 else []
        return []


    rule star_align:
        """Run STAR alignment on a single sample to generate junction quantification.
        
        Input: FASTQ files (single-end or paired-end)
        Output: Junction quantification TSV in GDC-compatible format
        """
        input:
            index_dir="resources/star_index",
            index_done="resources/star_index/index.done",
            fastq1=get_fastq1,
            fastq2=get_fastq2,
        output:
            junctions=os.path.join(OUT["raw_data"], "local", "files", "{sample}.tsv"),
            done=touch(os.path.join(OUT["raw_data"], "local", "files", "{sample}.done")),
        log:
            os.path.join(OUT["logs"], "star", "{sample}_align.log"),
        params:
            output_prefix=lambda w: os.path.join(OUT["raw_data"], "local", "files", f"{w.sample}_"),
        threads: 8
        resources:
            mem_mb=32000,
        conda:
            "../envs/star.yaml"
        shell:
            """
            # Determine read files command for gzipped input
            READCMD=""
            if [[ "{input.fastq1}" == *.gz ]]; then
                READCMD="--readFilesCommand zcat"
            fi
            
            # Determine if paired-end
            FASTQ_FILES="{input.fastq1}"
            if [[ -n "{input.fastq2}" ]]; then
                FASTQ_FILES="{input.fastq1} {input.fastq2}"
            fi
            
            # Run STAR
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
            
            # Convert STAR SJ.out.tab to pipeline format
            # STAR SJ.out.tab: chrom, start, end, strand(0/1/2), motif, annotated, unique, multi, overhang
            # Pipeline format: junction_id(chr:start:end:strand) \\t mapped_reads
            awk -F'\\t' '{{
                strand = ".";
                if ($4 == 1) strand = "+";
                else if ($4 == 2) strand = "-";
                if ($7 > 0) print $1":"$2":"$3":"strand"\\t"$7
            }}' {params.output_prefix}SJ.out.tab > {output.junctions}
            
            echo "Converted $(wc -l < {output.junctions}) junctions from STAR output"
            """


    rule create_local_manifest:
        """Create a manifest TSV from the local samples configuration.
        
        This mirrors the GDC manifest format so downstream rules can use the
        same filtering logic regardless of data source.
        """
        output:
            manifest=os.path.join(OUT["raw_data"], "local", "manifest.tsv"),
        log:
            os.path.join(OUT["logs"], "star", "manifest.log"),
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
                    sample_id = row["sample_id"]
                    sample_type = row.get("sample_type", "Unknown")
                    fout.write(f"{sample_id}\t{sample_id}.tsv\t{sample_type}\tlocal\n")


    # Checkpoint to collect all aligned samples
    checkpoint local_alignment_complete:
        """Checkpoint that triggers after all local samples are aligned.
        
        This replaces the GDC download checkpoint for local mode.
        """
        input:
            manifest=os.path.join(OUT["raw_data"], "local", "manifest.tsv"),
            samples=lambda w: expand(
                os.path.join(OUT["raw_data"], "local", "files", "{sample}.tsv"),
                sample=[s["sample_id"] for s in get_local_samples()]
            ),
        output:
            done=os.path.join(OUT["raw_data"], "local", "download.done"),
            data_dir=directory(os.path.join(OUT["raw_data"], "local", "files")),
        shell:
            "touch {output.done}"
