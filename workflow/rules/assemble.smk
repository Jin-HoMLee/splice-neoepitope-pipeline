# =============================================================================
# Rule module: Step 3 — Assemble 50 nt nucleotide contigs around splice junctions
# =============================================================================

rule assemble_contigs:
    """For each novel splice junction, extract a 50 nt contig centred on the
    junction (26 nt upstream + 24 nt downstream) using bedtools getfasta.
    Contigs that contain soft-clipped regions are excluded.
    Output: FASTA file of contigs per cancer type."""
    input:
        novel_junctions=rules.filter_junctions.output.novel_junctions,
        genome_fasta=config["reference"]["genome_fasta"],
    output:
        contigs_fasta=os.path.join(
            OUT["contigs"], "{cancer_type}", "contigs.fa"
        ),
    log:
        os.path.join(OUT["logs"], "assemble", "{cancer_type}_assemble.log"),
    params:
        upstream_nt=config["assembly"]["upstream_nt"],
        downstream_nt=config["assembly"]["downstream_nt"],
    conda:
        "../envs/biotools.yaml"
    script:
        "../scripts/assemble_contigs.py"
