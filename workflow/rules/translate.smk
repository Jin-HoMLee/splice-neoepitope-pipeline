# =============================================================================
# Rule module: Step 4 — Extract junction-spanning 9-mers from contigs
# =============================================================================

rule translate_peptides:
    """Extract junction-spanning 9-mer peptides directly from 50 nt contigs.
    For each contig, all 27 nt windows satisfying the spanning condition
    (first codon fully upstream, last codon fully downstream) are translated
    in all three reading frames.
    Output: TSV file of junction-spanning 9-mers per patient."""
    input:
        contigs_fasta=rules.assemble_contigs.output.contigs_fasta,
    output:
        peptides_tsv=os.path.join(
            OUT["peptides"], "{patient_id}", "peptides.tsv"
        ),
    log:
        os.path.join(OUT["logs"], "translate", "{patient_id}_translate.log"),
    params:
        upstream_nt=config["assembly"]["upstream_nt"],
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/translate_peptides.py"
