# =============================================================================
# Rule module: Step 6 — TCRdock structural validation (optional, GPU only)
# =============================================================================
#
# Predicts the 3D structure of the TCR-peptide-MHC ternary complex for the
# top-ranked neoepitope candidate(s).
#
# Requirements
# ------------
#   - Linux x64 with NVIDIA GPU (AlphaFold v2 backend)
#   - Docker with access to the configured TCRdock image at
#     config[tcrdock][docker_image]
#
# This rule is skipped entirely when config[tcrdock][enabled] is false,
# so local / macOS arm64 / CPU-only runs are unaffected.
#
# Outputs
# -------
#   results/predictions/{patient_id}/tcrdock/top_candidate.pdb   — predicted ternary complex
#   results/predictions/{patient_id}/tcrdock/docking_scores.tsv  — docking geometry metrics

_TCRDOCK_ENABLED = config.get("tcrdock", {}).get("enabled", False)


if _TCRDOCK_ENABLED:

    rule run_tcrdock:
        """Run TCRdock on the top strong-binding neoepitope candidate.

        Selects the top N candidates by IC50 from MHCflurry predictions,
        builds a TCRdock input TSV using the configured (or fallback) HLA
        allele and TCR sequences, runs TCRdock, and writes the predicted
        PDB structure and docking geometry scores.
        """
        input:
            predictions_tsv=rules.run_mhcflurry.output.predictions_tsv,
        output:
            pdb=os.path.join(
                OUT["predictions"], "{patient_id}", "tcrdock", "top_candidate.pdb"
            ),
            scores_tsv=os.path.join(
                OUT["predictions"], "{patient_id}", "tcrdock", "docking_scores.tsv"
            ),
        log:
            os.path.join(OUT["logs"], "tcrdock", "{patient_id}_tcrdock.log"),
        params:
            n_candidates=config["tcrdock"]["n_candidates"],
            docker_image=config["tcrdock"]["docker_image"],
            fallback_hla=config["tcrdock"]["fallback_hla"],
            fallback_tcr=config["tcrdock"]["fallback_tcr"],
        conda:
            "../envs/python.yaml"
        script:
            "../scripts/run_tcrdock.py"

    rule generate_report_with_structure:
        """Generate HTML report including embedded Mol* 3D structure viewer.

        Extends the standard report with an interactive TCR-peptide-MHC
        structure for the top candidate, with the PDB inlined as text so
        the report remains fully self-contained.
        """
        input:
            novel_junctions=rules.filter_junctions.output.novel_junctions,
            predictions_tsv=rules.run_mhcflurry.output.predictions_tsv,
            contigs_fasta=rules.assemble_contigs.output.contigs_fasta,
            pdb=rules.run_tcrdock.output.pdb,
            scores_tsv=rules.run_tcrdock.output.scores_tsv,
            hla_qc=_hla_qc_input,
        output:
            report_html=os.path.join(
                OUT["reports"], "{patient_id}", "report.html"
            ),
        log:
            os.path.join(OUT["logs"], "report", "{patient_id}_report.log"),
        params:
            ic50_strong=config["mhcflurry"]["ic50_strong"],
        conda:
            "../envs/python.yaml"
        script:
            "../scripts/generate_report.py"
