# =============================================================================
# Rule module: Step 3 — Assemble nucleotide contigs around splice junctions
# =============================================================================

_FLANK_NT = 3 * (max(config["translation"]["peptide_lengths"]) - 1)

rule assemble_contigs:
    """For each novel splice junction, extract a contig centred on the junction
    using bedtools getfasta. Flank size is derived from the maximum configured
    peptide length: flank_nt = 3 * (max_length - 1).
    Contigs that contain soft-clipped regions are excluded.
    Output: FASTA file of contigs per patient."""
    input:
        novel_junctions=rules.filter_junctions.output.novel_junctions,
        genome_fasta=config["reference"]["genome_fasta"],
    output:
        contigs_fasta=os.path.join(
            _RES, "{patient_id}", "contigs", "contigs.fa"
        ),
        stats=os.path.join(
            _RES, "{patient_id}", "contigs", "contig_assemble_stats.tsv"
        ),
    log:
        os.path.join(_LOGS, "{patient_id}", "assemble_contigs", "assemble.log"),
    params:
        upstream_nt=_FLANK_NT,
        downstream_nt=_FLANK_NT,
    conda:
        "../envs/biotools.yaml"
    script:
        "../scripts/assemble_contigs.py"
