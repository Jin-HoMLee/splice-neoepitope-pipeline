# Issue #204 — `fetch_vdjdb_panel` Snakemake rule + stitchr/VDJdb setup

**Date:** 2026-05-20
**Issue:** [Issue #204](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/204) (VDJdb panel rule)
**Parent epic:** [Issue #86](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/86) (HLA-matched TCR panel)
**Spinoff filed:** [Issue #435](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/435) (flatten `predictions/` wrapper)
**Author:** Developer

## Goal

Add a new Snakemake rule `fetch_vdjdb_panel` that produces, per patient, a **reference TCR panel** of up to 10 HLA-matched α/β TCRs per patient HLA allele — sourced from VDJdb, reconstructed to full chain sequences via stitchr. The panel is the foundation for sub-issue [#205](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/205) (HLA-matched TCR selection) and sub-issue [#206](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/206) (report surfacing).

## Empirical pre-step (completed)

Before writing any code, the issue's pre-implementation step required exact-match counts per patient allele. Results in [Issue #204 comment](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/204#issuecomment-4502381646):

- Pin: VDJdb release `2026-05-16` (Blooming May)
- 8,901 paired α/β rows globally pass `HomoSapiens + MHCI + score≥2`
- **patient_001:** 3/6 alleles empty (`A*31:01`, `A*26:01`, `B*15:63`), 2 low-coverage, 1 ok
- **patient_002:** 0 empty, 3 ok, 2 low-coverage
- Conclusion: ship exact-match-only as designed; supertype fallback motivated as high-priority follow-up

## Architecture

```
┌──────────────────────────────┐
│ download_vdjdb_release       │ ← one-time, sentinel + SHA256 check
│ (download.smk)               │
└─────────────┬────────────────┘
              │ vdjdb_full.txt
              ▼
┌──────────────────────────────┐
│ download_imgt_germlines      │ ← one-time, stitchrdl
│ (download.smk)               │
└─────────────┬────────────────┘
              │ IMGT HUMAN data
              ▼
results/{p}/hla_typing/alleles.tsv ──┐
                                     ▼
                       ┌─────────────────────────────┐
                       │ fetch_vdjdb_panel           │ ← per-patient, gated on hla.enabled
                       │ (tcr_panel.smk, NEW file)   │
                       │ env: vdjdb.yaml (NEW)       │
                       │ script: fetch_vdjdb_panel.py│
                       └─────────────────────────────┘
                                     │
                                     ▼
                results/{p}/tcr_panel/vdjdb/
                  ├── panel.tsv
                  └── panel_qc.tsv
```

## File-by-file changes

### 1. `workflow/envs/vdjdb.yaml` (NEW)

```yaml
name: vdjdb
channels: [conda-forge, bioconda]
dependencies:
  - python=3.11
  - pandas
  - requests
  - pip
  - pip:
    - stitchr>=1.1.0
    - IMGTgeneDL>=0.7.0
```

Dedicated env per the established pattern (`optitype.yaml`, `hisat2.yaml`, `star.yaml`). Local-runnable on M1 (no GPU/Docker dependency).

### 2. `workflow/rules/download.smk` — two new rules

**`download_vdjdb_release`** — pulls `vdjdb-{release}.zip` from `antigenomics/vdjdb-db` GitHub release, verifies SHA256, extracts `vdjdb_full.txt`, writes sentinel.

- Inputs: none (uses `config[tcrdock][vdjdb_release]` + `config[tcrdock][vdjdb_sha256]`)
- Outputs:
  - `resources/vdjdb/{release}/vdjdb_full.txt`
  - `resources/vdjdb/{release}/.download.done`
- Uses `vdjdb.yaml` env (requests + python).

**`download_imgt_germlines`** — runs `stitchrdl --species HUMAN`.

- Inputs: none
- Outputs: `resources/imgt_germlines/.download.done`
- Uses `vdjdb.yaml` env (provides `stitchrdl` via IMGTgeneDL).

**Idempotency:** Sentinel-gated; re-runs are no-ops. Old release dirs preserved on version bump (manual cleanup).

### 3. `workflow/rules/tcr_panel.smk` (NEW file)

One rule, `fetch_vdjdb_panel`. Per-patient, gated on `config[hla][enabled]`.

```python
rule fetch_vdjdb_panel:
    input:
        vdjdb_tsv = lambda wc: f"resources/vdjdb/{config['tcrdock']['vdjdb_release']}/vdjdb_full.txt",
        vdjdb_sentinel = lambda wc: f"resources/vdjdb/{config['tcrdock']['vdjdb_release']}/.download.done",
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
```

Add `include: "workflow/rules/tcr_panel.smk"` to the top-level `Snakefile`.

### 4. `workflow/scripts/fetch_vdjdb_panel.py` (NEW)

Lazy-imports pandas + stitchr at first use (per [PR #428](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/428) lazy-import convention). Pipeline:

1. **Load alleles** from `alleles.tsv` (existing schema: `allele1`, `allele2` columns × 3 loci → 6 alleles, deduplicated).
2. **Load + filter VDJdb** (`vdjdb_full.txt`): `species == "HomoSapiens"`, `mhc.class == "MHCI"`, `vdjdb.score >= min_score`.
3. **Normalize `mhc.a` to 4-digit**: truncate to first two `:`-separated fields (`HLA-A*02:01:110` → `HLA-A*02:01`; `HLA-A*02` excluded as 2-digit).
4. **For each patient allele:**
   - Exact-match filter to get all candidate rows.
   - Sort by `vdjdb.score DESC`, then `meta.subject.id ASC` (deterministic tiebreak).
   - Iterate top-down: call `stitchr` on each candidate → `alpha_seq`, `beta_seq`. On stitch failure: log ERROR (with allele + donor ID + V/J/CDR3), skip row, continue.
   - Stop when `panel_size` successes accumulated OR candidate pool exhausted (whichever first). Resulting count goes into `panel_qc.tsv.n_in_panel`.
5. **Write `panel.tsv`** — 11 columns: `allele, va_gene, ja_gene, cdr3a, vb_gene, jb_gene, cdr3b, alpha_seq, beta_seq, vdjdb_score, vdjdb_donor_id`.
6. **Write `panel_qc.tsv`** — 4 columns: `allele, n_exact_matches, n_in_panel, panel_status` (where `panel_status ∈ {ok, low_coverage, empty}`).
7. **Log summary**: `"VDJdb panel built for {patient} — {n_alleles} alleles: {n_ok} ok, {n_low} low_coverage, {n_empty} empty."`

### 5. `config/config.yaml` — additions to existing `tcrdock:` block

```yaml
tcrdock:
  enabled: false  # existing
  # New keys for VDJdb panel fetch (Issue #204):
  vdjdb_release: "2026-05-16"
  vdjdb_sha256: "0dce79ec55c109000da10b7bc72300e352ffb7df92d5de30682de20bab35a366"
  vdjdb_min_score: 2
  vdjdb_panel_size: 10
```

### 6. `research/glossary.md` — new entry

Add `## K` section between `## I` and `## L`:

```markdown
## K

**KIR** — Killer-cell Immunoglobulin-like Receptor. Family of NK-cell receptors whose dominant ligands are HLA-C alleles; HLA-C–KIR engagement inhibits NK cytotoxicity under "missing-self" surveillance. Relevant here because HLA-C is studied primarily through this NK biology rather than T-cell presentation, which explains why TCR repertoire databases (VDJdb, McPAS-TCR) have thin paired α/β coverage for C-locus alleles. See [Issue #204](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/204) empirical coverage check. *Domain: bio.*
```

## Output layout (forward-compatible)

```
results/{patient_id}/
  tcr_panel/                        ← NEW top-level, source-keyed
    vdjdb/
      panel.tsv                     ← this PR
      panel_qc.tsv                  ← this PR
    trust4/                         ← future, parent epic Issue #24
    protcr/                         ← future, parent epic Issue #24
  predictions/                      ← will flatten via Issue #435
    tcrdock/
      ...
```

**Why `tcr_panel/<source>/` (not `predictions/tcrdock/vdjdb_panel.tsv`):** TCR panels are *candidate sets* (input to TCRdock), not prediction outputs. The source-keyed shape makes [Issue #24](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/24) (TRUST4 + ProTCR) drop in as a sibling subdir with the same `panel.tsv` + `panel_qc.tsv` filenames — uniform schema across sources.

## DAG re-fire behavior

| Change | Re-fires `fetch_vdjdb_panel`? | Re-fires downloads? |
|---|---|---|
| New patient | ✅ for that patient | ❌ (sentinel exists) |
| Bump `vdjdb_release` | ✅ all patients (new release dir → new input path) | ✅ new release downloaded; old dir kept |
| Bump `vdjdb_min_score` or `vdjdb_panel_size` | ✅ all patients (params change) | ❌ |
| `alleles.tsv` changes | ✅ that patient | ❌ |
| Re-run, no changes | ❌ (mtime no-op) | ❌ |

## Error handling

| Failure mode | Behavior |
|---|---|
| VDJdb URL 404 / network error | Rule fails loudly; no sentinel written |
| VDJdb SHA256 mismatch | Rule fails with explicit mismatch — strong reproducibility guarantee |
| `stitchrdl` IMGT download fails | Rule fails; manual recovery by re-running `download_imgt_germlines` |
| `alleles.tsv` missing | Cannot happen if `hla.enabled = true` (Snakemake DAG enforces); if false, rule excluded from DAG |
| Allele has 0 matches | WARNING log; `panel_qc.tsv` row with `panel_status=empty`; no rows in `panel.tsv` for that allele; **no exception** (expected case per pre-step) |
| Allele has 1–9 matches | WARNING log; `panel_status=low_coverage` |
| stitchr fails for a single row | ERROR log; skip row; if drops allele below top-N, fall back to next-best entries |
| All entries fail to stitch | Rule fails with summary count |
| VDJdb schema changes | Script reads columns by name; missing required column raises `KeyError` with the offending column |

## Testing

### Unit tests (`workflow/tests/test_fetch_vdjdb_panel.py`)

- `test_normalize_allele_to_4digit` — table-driven, ~6 cases
- `test_filter_paired_alpha_beta` — synthetic 5-row VDJdb fixture
- `test_filter_homosapiens_mhci_score` — species/class/score combinations
- `test_exact_match_per_allele` — `HLA-A*02:01` matches `HLA-A*02:01:110`, not `HLA-A*02:02`
- `test_top_n_deterministic_tiebreak` — tied scores resolved by `meta.subject.id` lex order
- `test_panel_status_classification` — `0 → empty`, `1–9 → low_coverage`, `≥10 → ok`
- `test_empty_allele_emits_warning` — `caplog` assertion on warning line
- `test_stitchr_smoke` — real `stitchrdl`-cached IMGT data + DMF5 (`TRAV12-2 / TRAJ21 / CAVNFGGGKLI`); assert stitched α ends with CDR3 and begins with TRAV12-2 leader (lightweight sanity, not full equality)
- `test_output_schema_columns` — `panel.tsv` has the 11 spec'd columns
- `test_qc_schema_columns` — `panel_qc.tsv` has 4 columns

### Integration

- **Dry-run (CI):** `snakemake -n --use-conda --configfile config/config.yaml -- results/test/patient_001/tcr_panel/vdjdb/panel.tsv` (using the `--` terminator per CLAUDE.md's Snakemake 8 gotcha) — piggybacks on the existing `pipeline-snakemake-dry-run` CI gate.
- **Real run (local M1):** end-to-end on the chr22 test dataset post-merge.

### Fixtures

- `workflow/tests/fixtures/vdjdb_mini.tsv` — synthetic ~10-row VDJdb table (35 columns) checked in.
- IMGT germline data — lazy-downloaded by pytest fixture on first run, cached. Marked `@pytest.mark.network`.

## Acceptance criteria (from Issue #204)

- [x] Empirical pre-implementation step done — exact-match counts surfaced for patient_001 + patient_002 ([comment](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/204#issuecomment-4502381646))
- [ ] `stitchr` installed (via `vdjdb.yaml` env) and IMGT data accessible
- [ ] VDJdb release `2026-05-16` pinned in config and downloaded reproducibly (SHA256-verified)
- [ ] `fetch_vdjdb_panel` rule produces valid panel for patient_001 (6 alleles)
- [ ] `panel_qc.tsv` records per-allele match counts and `panel_status`
- [ ] Rare/zero-match alleles emit WARNING lines as specified
- [ ] Pytest covers parse + filter + stitchr smoke test on a known TCR

## Out of scope (per Issue #204)

- **Supertype / 2-digit matching fallback** — deferred; the QC logging above is the empirical input that decides whether a follow-up issue is needed (likely yes given patient_001's 3/6 empty alleles).
- CDR3 clustering for diversity-based selection — deferred.
- Wiring panel into `report.html` — sub-issue [#206](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/206).
- TCR selection from panel → TCRdock input — sub-issue [#205](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/205).
- Flattening `results/{p}/predictions/` wrapper — [Issue #435](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/435).

## Known limitations

- **IMGT germline data has no version pin.** `stitchrdl` always pulls the current IMGT release. This is a reproducibility gap that follows from stitchr's design. Worth tracking as a separate follow-up; not blocking for this PR.
- **3/6 patient_001 alleles will hit DMF5 fallback** under exact-match-only behavior. Surfaced cleanly via `panel_status=empty` for the eventual supertype-fallback follow-up.
