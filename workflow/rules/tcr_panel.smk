# =============================================================================
# Rule module: TCR panel construction (Issue #204)
# =============================================================================
#
# Per-patient: filters VDJdb for HLA-matched paired α/β TCRs, reconstructs
# full chain sequences via stitchr, and writes a reference TCR panel.
#
# Outputs
# -------
#   results/{patient_id}/tcr_panel/vdjdb/panel.tsv      — top-10 per allele
#   results/{patient_id}/tcr_panel/vdjdb/panel_qc.tsv   — per-allele coverage
#
# Gating
# ------
# Gated on config[hla][enabled]. When HLA typing is disabled, alleles.tsv
# isn't produced and this rule is excluded from the DAG. Independent of
# config[tcrdock][enabled] — the panel is a reference set, not a GPU output.

_HLA_ENABLED = config.get("hla", {}).get("enabled", False)


if _HLA_ENABLED:

    rule fetch_vdjdb_panel:
        """Build the per-patient VDJdb TCR panel."""
        input:
            vdjdb_tsv = f"resources/vdjdb/{config['tcrdock']['vdjdb_release']}/vdjdb_full.txt",
            vdjdb_sentinel = f"resources/vdjdb/{config['tcrdock']['vdjdb_release']}/.download.done",
            imgt_sentinel = "resources/imgt_germlines/.download.done",
            alleles_tsv = os.path.join(_RES, "{patient_id}", "hla_typing", "alleles.tsv"),
        output:
            panel = os.path.join(_RES, "{patient_id}", "tcr_panel", "vdjdb", "panel.tsv"),
            qc = os.path.join(_RES, "{patient_id}", "tcr_panel", "vdjdb", "panel_qc.tsv"),
        log:
            os.path.join(_LOGS, "{patient_id}", "tcr_panel", "fetch_vdjdb_panel.log"),
        params:
            min_score = config["tcrdock"]["vdjdb_min_score"],
            panel_size = config["tcrdock"]["vdjdb_panel_size"],
        conda:
            "../envs/vdjdb.yaml"
        script:
            "../scripts/fetch_vdjdb_panel.py"
