# =============================================================================
# Rule module: Step 4 — Extract junction-spanning peptides from contigs
# =============================================================================

rule translate_peptides:
    """Extract junction-spanning peptides from junction contigs.
    For each configured peptide length, all windows satisfying the spanning
    condition (first codon fully upstream, last codon fully downstream) are
    translated in all three reading frames.
    Output: TSV file of junction-spanning peptides per patient."""
    input:
        contigs_fasta=rules.assemble_contigs.output.contigs_fasta,
    output:
        peptides_tsv=os.path.join(
            _RES, "{patient_id}", "peptides", "peptides.tsv"
        ),
    log:
        os.path.join(_LOGS, "{patient_id}", "translate_peptides", "translate.log"),
    params:
        upstream_nt=_FLANK_NT,
        peptide_lengths=config["translation"]["peptide_lengths"],
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/translate_peptides.py"
