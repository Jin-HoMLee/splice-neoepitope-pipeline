# Issue #566 — ASNEO cross-check (open-only, MHCflurry-standardized)

**Parent Issue:** [#566](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/566) (ASNEO cross-check) — itself a single-caller slice of the broader [#679](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/679) open caller benchmark; parked from the [#546](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/546) ASNEO desk eval (verdict: component reuse + cross-check ⚠️).

**Status:** 🟢 prepped — env built + validated, ASNEO architecture source-verified, option-B patcher + turnkey chr22 runner written & tested. ⛔ **The run itself is VM-bound** (needs STAR; not usable on the arm64 laptop, project policy = STAR is VM-only) — same constraint as [#636](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/636). No concordance numbers yet.

## Goal

Run **ASNEO** (our closest published peer: RNA-seq → splice junctions → neoepitopes) as an orthogonal second opinion on the same input, and compare its splice-neoepitope calls against our pipeline. The MHC step is **held constant on MHCflurry** so the comparison isolates *junction-detection + translation* differences rather than confounding them with MHC-predictor differences.

## Scope decision — MHC path: MHCflurry-swap, layered (option B, 2026-06-23)

1. **Primary (MHC-agnostic):** compare the *peptide candidate sets* — junctions detected, frames translated, peptides emitted (Jaccard overlap + divergence breakdown).
2. **Secondary (ranking):** score *both* candidate sets through the **same MHCflurry** presentation predictor; compare ranking concordance.
3. **Native NetMHCpan/NetCTLpan path demoted off the critical path** to optional supplementary (only if a license is ever obtained).

Rationale: isolates the junction-detection signal · open-only/reproducible (Apache-2.0) · consistent with #679 · removes the NetMHCpan license gate.

## Verified facts about ASNEO (source-checked 2026-06-23)

| Item | Finding | Source |
|------|---------|--------|
| Canonical repo | **`bm2-lab/ASNEO`** (Qi Liu's bm2-lab, Tongji — the *Aging* 2020 authors; pushed 2021-04-09). `zzb23/ASNEO` is the original-author 2019 mirror, same code. | `gh repo view` |
| Genome build | **hg19 / GRCh37** (`hg19.fa` or `GRCh37.fa`) | README |
| Input format | **STAR `SJ.out.tab`** (+ optional indexed BAM) | README + `test/run_ASNEO.sh` |
| Invocation | `python ASNEO.py -j <SJ.out.tab> -a <HLA,...> -g hg19.fa` | `test/run_ASNEO.sh` |
| Python deps | python3 + sj2psi, pysam, pandas, biopython, sklearn, xgboost; bedtools | README |
| **Biopython pin** | **< 1.80** — `ASNEO.py:11-12` imports `Bio.SubsMat` (removed 1.80) + `pairwise2` (removed 1.84) | `ASNEO.py:11-12` |
| Bundled binaries | `src/software.tar.gz` → `netMHCpan-4.0`, `netCTLpan-1.1`, `pepmatch_db_x86_64` (non-redistributable) | `ASNEO.py:446-448` |
| **Candidate peptides** | written to `putative_peptide.txt` at `ASNEO.py:252`, right before the netMHCpan call (`:256`) — **the MHCflurry-swap interception point** | `ASNEO.py:204,252,463` |
| Candidate set is **normal-subtracted** | `putative_peptide.txt` = `peps - norm_peps` (k-mers minus `Norm_Protein.fasta` proteome, `:246-249`) — i.e. ASNEO's junction-derived novel peptides *after* its self-filter. The right artifact to compare. | `ASNEO.py:246-253` |
| **Allele-independent under option B** | ASNEO's `-a` HLA arg is consumed *only* by the bypassed netMHCpan/netCTLpan steps → candidate generation needs no real HLA (use a placeholder); real HLA enters only at the downstream MHCflurry step | `ASNEO.py:256,378+` |
| MHC call sites (bypassed) | netMHCpan `:218`, netCTLpan `:280`, netMHCpan `:290`, pepmatch wild-type-match `:270` (used only in the bypassed `ProcessNeo`) | `ASNEO.py` |
| `tmp/` is deleted at end | `main()` `shutil.rmtree(path['tdir'])` — so the patch copies `putative_peptide.txt` to the outdir before cleanup | `ASNEO.py:485` |

## Liftover plan — re-align, don't lift coordinates

ASNEO consumes `SJ.out.tab` natively, so the robust path is to **generate a native hg19 `SJ.out.tab`** by re-aligning the same FASTQs with STAR against an hg19 index — *not* to liftover our GRCh38 junction coordinates (cross-build splice-junction liftover is error-prone). Our pipeline's GRCh38 calls stay as-is on our side; ASNEO gets its own hg19 alignment. Both originate from identical reads, so the comparison is fair.

- **First input = chr22 PoC** (smoke): build a chr22-only STAR hg19 index (single-chr index fits in laptop RAM, unlike the whole-genome >8 GB build), align the existing chr22 test FASTQs, feed the resulting `SJ.out.tab` to ASNEO.
- Then a patient (heavier; VM-bound).

## Files & how to run

| File | Role |
|------|------|
| `asneo_env.yml` | Conda env (biopython=1.79, python=3.8 + bedtools); **built & import-validated** locally. Bundled MHC binaries not installed (option B). (At the experiment root, not `env/` — that collides with the `ENV/` gitignore rule on a case-insensitive FS.) |
| `apply_optionB_patch.py` | String-anchored patcher: makes ASNEO stop after `putative_peptide.txt`, copy it to the outdir, repoint `bedtools` to PATH (skips `software.tar.gz`), and skip `ProcessNeo`. Self-verifying (each anchor must match once). **Tested against the real `ASNEO.py` — applies + compiles clean.** |
| `run_asneo_chr22.sh` | Turnkey chr22 PoC runner (VM): fetch hg19 chr22 → STAR hg19 chr22 index → align `SRR9143066` → patch+run ASNEO → `putative_peptide.txt`. `bash -n` clean. |
| `notebook.ipynb` | *(to add at run time)* MHCflurry concordance: score ASNEO's + our candidate peptides via the same predictor; primary call-concordance + secondary ranking. |

**To run (on the VM):** `conda env create -f asneo_env.yml` (once), then `conda activate asneo && bash research/experiments/issue_566_asneo_crosscheck/run_asneo_chr22.sh` from the repo root. STAR must be on PATH.

## Outputs index

`outputs/` — (empty; populated on first run: ASNEO `putative_peptide.txt`, candidate-set concordance tables, MHCflurry ranking comparison, figures).

## Open sub-questions

- ~~Keep ASNEO's `pepmatch` self-proteome filter as part of its "call"?~~ **Resolved:** the candidate set already has the normal-proteome subtraction applied (`peps - norm_peps`, `:246-249`); the later `pepmatch` (`:270`) is a *wild-type-match* step inside the bypassed `ProcessNeo`, not part of candidate definition. `putative_peptide.txt` is the right comparison artifact.
- Patient input choice + HLA typing source for the **secondary** MHCflurry step (the primary call-concordance is MHC-agnostic, so this only gates the ranking comparison).
- Confirm chr22-only genome FASTA + ASNEO's whole-genome `data/iRefSeq.bed` interplay at run time (chr22 junctions → chr22 isoforms only; expected fine, verify on first run).

## Cross-experiment deps

- Reuses our existing chr22 STAR/HISAT2 test fixtures and MHCflurry setup. No write into other experiments' `outputs/`.
