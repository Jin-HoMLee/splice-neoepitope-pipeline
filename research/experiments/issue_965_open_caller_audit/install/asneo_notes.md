# ASNEO - open-only recipe pointer (already engineered by #566)

ASNEO (`bm2-lab/ASNEO`, Apache-2.0) is an open-only **GO**. The open-only install decision and the NetMHCpan->MHCflurry swap were already built and locally validated in [#566](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/566) - do not re-derive. Everything lives in [`../../issue_566_asneo_crosscheck/`](../../issue_566_asneo_crosscheck/).

## What #566 already delivers (reuse verbatim)

- **`asneo_env.yml`** - conda env (python 3.8, biopython 1.79 pinned <1.80, bedtools 2.31.1, sj2psi), **built + import-validated on osx-arm64**. The `asneo` conda env already exists on this box. The bundled non-redistributable binaries are deliberately NOT installed.
- **`apply_optionB_patch.py`** - tested, string-anchored patcher: ASNEO stops after writing its normal-subtracted candidate peptides (`putative_peptide.txt`, `ASNEO.py:252`), copies them out before `tmp/` cleanup, repoints bedtools to PATH so `src/software.tar.gz` (netMHCpan-4.0 / netCTLpan-1.1 / pepmatch) is **never extracted**, and skips `ProcessNeo`. Applies + compiles clean against the real `ASNEO.py`.
- **`run_asneo_chr22.sh`** - turnkey chr22 PoC runner (STAR hg19 -> `SJ.out.tab` -> patched ASNEO).

## Open-only run (this issue)

```bash
conda activate asneo
git clone https://github.com/bm2-lab/ASNEO && \
  python research/experiments/issue_566_asneo_crosscheck/apply_optionB_patch.py ASNEO/ASNEO.py
python ASNEO/ASNEO.py -j <SJ.out.tab> -a HLA-A02:01 -g <hg19.fa> -o <outdir>
#   -> <outdir>/putative_peptide.txt   (junction-derived, normal-subtracted candidates)
# then score putative_peptide.txt with MHCflurry (the still-to-write concordance notebook)
```

## Where it runs
- **Caller-proper (patched ASNEO): local arm64 CPU**, given a pre-made hg19 `SJ.out.tab`. `-a HLA` is a placeholder here - under option B it is consumed only by the bypassed MHC steps; real HLA enters at the downstream MHCflurry step.
- **Front-end alignment (STAR -> `SJ.out.tab`): VM-bound** (STAR's human index build exceeds 8 GB; STAR is VM-only per project policy). A chr22-only STAR index can fit 8 GB but the clean laptop path is a pre-generated junction table.

## Local smoke (this issue, PASS)

`install/asneo_smoke.sh` runs the patched (option-B, no netMHCpan) caller on **ASNEO's own bundled test `SJ.out.tab`** (`test/SRR2660032.SJ.out.tab`) subset to chr22, against the small hg19 chr22 FASTA - no STAR, no hand-crafted input. On arm64 CPU in the `asneo` env it runs all 11 stages end to end (exit 0) and emits junction-derived, normal-subtracted candidate peptides:

- **Default thresholds** (`--reads 10 --psi 0.1 -l 9`): 6194 chr22 junctions -> 11 pass filters -> 1 novel isoform -> 0 nine-mers (the test data has low chr22 coverage; a scale artifact, not a failure - the pipeline still writes its `putative_peptide.txt`).
- **Relaxed thresholds** (`--reads 2 --psi 0.05 -l 8,9,10,11`): 6194 -> 110 junctions -> 60 novel isoforms -> **800 candidate peptides** (`outputs/asneo_smoke_peptides_head.txt`). Confirms the full translate -> k-mer -> normal-subtract -> write tail runs and emits the artifact that would feed MHCflurry.

So the open-only ASNEO caller-proper is now demonstrated running locally on arm64, not just env-validated. The bundled netMHCpan/netCTLpan binaries are never extracted (option-B); the only remaining open-only work is the MHCflurry scoring of `putative_peptide.txt`.

## Remaining open-only work
Only the MHCflurry scoring of `putative_peptide.txt` (the #566 concordance notebook, tracked by #848). No license or install work remains.

## Openness
ASNEO code: Apache-2.0. MHCflurry: Apache-2.0. The vendored netMHCpan-4.0 / netCTLpan-1.1 / pepmatch binaries are non-redistributable and are **bypassed, never extracted or installed**.
