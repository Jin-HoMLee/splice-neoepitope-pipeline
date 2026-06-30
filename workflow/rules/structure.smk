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
#   results/{patient_id}/tcrdock/top_candidate.pdb   — predicted ternary complex
#   results/{patient_id}/tcrdock/docking_scores.tsv  — docking geometry metrics

_TCRDOCK_ENABLED = config.get("tcrdock", {}).get("enabled", False)
_HLA_ENABLED = config.get("hla", {}).get("enabled", False)


def _vdjdb_panel_input(wildcards):
    """Attach the per-patient VDJdb panel (Issue #204) as an input only when
    HLA typing is enabled — the fetch_vdjdb_panel rule that produces it is
    gated on config[hla][enabled]. When absent, run_tcrdock.py falls back to
    the DMF5 TCR (Issue #205)."""
    if _HLA_ENABLED:
        return {"vdjdb_panel": os.path.join(
            _RES, wildcards.patient_id, "tcr_panel", "vdjdb", "panel.tsv"
        )}
    return {}


if _TCRDOCK_ENABLED:

    rule run_tcrdock:
        """Run TCRdock on the top strong-binding neoepitope candidate.

        Selects the top N candidates ranked by genotype_presentation_score
        (quality-gated by presentation_percentile_weak) from MHCflurry predictions,
        picks the highest-confidence HLA-matched TCR from the VDJdb panel for each
        candidate's allele (DMF5 fallback when the panel has no match), runs
        TCRdock, and writes the predicted PDB structure and docking geometry scores.
        """
        input:
            unpack(_vdjdb_panel_input),
            # Calibrated TSV (apply_calibrator) — positions the calibrator rule
            # pre-TCRdock in the DAG (Issue #709). TCRdock ranks by the unchanged
            # genotype_presentation_score; the added columns are carried through.
            predictions_tsv=rules.apply_calibrator.output.calibrated_tsv,
        output:
            pdb=os.path.join(
                _RES, "{patient_id}", "tcrdock", "top_candidate.pdb"
            ),
            scores_tsv=os.path.join(
                _RES, "{patient_id}", "tcrdock", "docking_scores.tsv"
            ),
        log:
            os.path.join(_LOGS, "{patient_id}", "structure", "tcrdock.log"),
        params:
            n_candidates=config["tcrdock"]["n_candidates"],
            docker_image=config["tcrdock"]["docker_image"],
            fallback_hla=config["tcrdock"]["fallback_hla"],
            fallback_tcr=config["tcrdock"]["fallback_tcr"],
            presentation_percentile_weak=config["mhcflurry"]["presentation_percentile_weak"],
        conda:
            "../envs/python.yaml"
        script:
            "../scripts/run_tcrdock.py"
