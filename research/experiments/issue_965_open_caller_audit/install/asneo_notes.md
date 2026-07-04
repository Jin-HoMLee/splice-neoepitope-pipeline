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

## Remaining open-only work
Only the MHCflurry scoring of `putative_peptide.txt` (the #566 concordance notebook, not yet written). No license or install work remains.

## Openness
ASNEO code: Apache-2.0. MHCflurry: Apache-2.0. The vendored netMHCpan-4.0 / netCTLpan-1.1 / pepmatch binaries are non-redistributable and are **bypassed, never extracted or installed**.
