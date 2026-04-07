# =============================================================================
# Rule module: Step 4 — In-silico translation of contigs to 16-mer peptides
# =============================================================================

rule translate_peptides:
    """Translate each 50 nt contig in all three reading frames to produce up to
    three 16-mer peptides per contig.  Stop codons truncate the peptide;
    start codons (M) are noted in the sequence but not used to truncate.
    Output: FASTA file of 16-mer peptides per cancer type."""
    input:
        contigs_fasta=rules.assemble_contigs.output.contigs_fasta,
    output:
        peptides_fasta=os.path.join(
            OUT["peptides"], "{cancer_type}", "peptides.fa"
        ),
    log:
        os.path.join(OUT["logs"], "translate", "{cancer_type}_translate.log"),
    params:
        reading_frames=config["translation"]["reading_frames"],
        peptide_length=config["translation"]["peptide_length"],
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/translate_peptides.py"
