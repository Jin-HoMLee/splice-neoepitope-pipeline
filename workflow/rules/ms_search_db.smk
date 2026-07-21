# =============================================================================
# Rule module: Optional - MS search database (junction-spanning ORF FASTA)
# =============================================================================
#
# Decoupled 30-nt contig assembly (reuses assemble_contigs.py) feeding the
# ORF-stretch FASTA emitter (orf_fasta_from_contigs.py). Optional target;
# NOT part of `rule all` and does not touch the production MHCflurry k-mer
# path (translate_peptides.smk stays on translation.peptide_lengths / 27-nt
# flanks). See config["ms_search_db"] for the 30-nt flank + min peptide len.

rule assemble_contigs_wide:
    """Re-run contig assembly at a wider 30-nt flank (independent of the
    production translation.peptide_lengths-derived flank) so every class-I
    8-11mer straddling the breakpoint is covered for the MS search path."""
    input:
        novel_junctions=rules.filter_junctions.output.novel_junctions,
        genome_fasta=config["reference"]["genome_fasta"],
    output:
        contigs_fasta=os.path.join(
            _RES, "{patient_id}", "ms_search_db", "contigs_wide.fa"
        ),
    log:
        os.path.join(_LOGS, "{patient_id}", "ms_search_db", "assemble_contigs_wide.log"),
    params:
        upstream_nt=config["ms_search_db"]["flank_nt"],
        downstream_nt=config["ms_search_db"]["flank_nt"],
    conda:
        "../envs/biotools.yaml"
    script:
        "../scripts/assemble_contigs.py"


rule emit_orf_fasta:
    """Translate each wide contig in all three frames and keep the stop-free
    stretch crossing the breakpoint, writing one FASTA record per
    junction/frame for a nonspecific MS search (Issue #1176)."""
    input:
        contigs_fasta=rules.assemble_contigs_wide.output.contigs_fasta,
    output:
        orf_fasta=os.path.join(
            _RES, "{patient_id}", "ms_search_db", "junction_orf.fasta"
        ),
    log:
        os.path.join(_LOGS, "{patient_id}", "ms_search_db", "emit_orf_fasta.log"),
    params:
        flank_nt=config["ms_search_db"]["flank_nt"],
        min_peptide_len=config["ms_search_db"]["min_peptide_len"],
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/orf_fasta_from_contigs.py"
