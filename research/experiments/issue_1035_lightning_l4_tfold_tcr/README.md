# Issue 1035 spike: free Lightning-AI GPU as a $0 path to AF3-class TCR-pMHC co-folding

Turnkey recipe for running **tFold-TCR** on a **free Lightning-AI GPU** at a hard $0 ceiling.
This is the desk-prep + live-setup record; the prediction timing/VRAM numbers are captured on the live run (see [Live-run findings](#live-run-findings-2026-07-20---complete)).

Tracks [Issue #1035](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1035); feeds the backend re-eval parked in [Issue #601](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/601).

## Freshness check (2026-07-20) - premise mostly holds, with one correction

Lightning-AI free tier, per its [billing FAQ](https://lightning.ai/docs/overview/faq/billing) + confirmed live at signup:

- **15 credits/month, no credit card required** - but **phone verification IS required** (so "$0" is not zero-friction).
- **L4 IS selectable, but pricier than the Issue assumed.** The machine picker offers CPU (free), **T4 $0.19/hr**, **L4 $1.58/hr**, **L40S $2.89/hr**, **RTX PRO 6000 $4.64/hr**, **A100 $2.19/hr**, **H100 $4.50/hr**, **H200 $6.53/hr** - all free-tier selectable (no card; drawn from the monthly credits). This *corrects* the Issue's premise twice: A100/H100/H200 are **not** paid-gated, and **L4 is $1.58/hr, not the assumed ~$0.60-0.70** -> **~9.5 free L4-hr/mo, not 21-25**.
  - Credit->hour (15 credits/mo, 1 credit ~ $1): T4 ~79 hr, **L4 ~9.5 hr**, **L40S ~5 hr**, A100 ~7 hr, H100 ~3 hr.
- **The right pick is L40S/A100, and the L4 is moot** - see the [card-fit finding](#live-run-findings-2026-07-20---complete): tFold-TCR peaks at **26 GB/complex**, so it **OOMs the 24 GB L4** (and the 16 GB T4). L40S (48 GB, bf16 + FP8) is both capable and has the VRAM headroom; A100 (40 GB) also fits.
- tFold-TCR runs in **seconds/complex** on the GPU, so L40S's ~5 free hr/mo is ample (~2,100 complexes/mo at $0).

## The two open gaps #601 flagged

| Gap | Status | Finding |
|-----|--------|---------|
| Weight/code license | **CLOSED** | Repo is [PolyForm Noncommercial 1.0.0](https://github.com/TencentAI4S/tfold/blob/master/LICENSE) - same terms for code and weights. Fine for our research/portfolio use; bars any commercial productization. |
| VRAM / bf16-FP8 wall | **CLOSED - no wall on L40S** | Free L40S runs Python 3.12 + **torch 2.8.0+cu128, CUDA available** - bf16/FP8-capable. **Peak VRAM = 26.0 GB/complex** (bundled example): fits L40S (48 GB) with headroom, but **exceeds L4's 24 GB** (would OOM). |

## Environment on Lightning (CORRECTED - do NOT `conda env create`)

**Lightning Studios allow only ONE conda env (the default `cloudspace`).**
`conda env create -f environment.yaml` is **blocked** ("max 1 environment") - so install tFold + its deps into the **default env** instead.

```bash
git clone https://github.com/TencentAI4S/tfold.git
cd tfold
pip install .                                                   # installs the tfold package into the default env
pip install termcolor biopython ml-collections dm-tree modelcif # tFold deps NOT already present
```

Notes learned live:
- The default env is **Python 3.12** with **torch 2.8.0+cu128** already installed - tFold imports and runs on 3.12 despite its `environment.yaml` pinning `python=3.8`.
- **Drop tFold's old version pins.** `numpy==1.21.2` / `biopython==1.79` are Python-3.8-era and have no 3.12 wheels; install the deps **unpinned** so pip picks 3.12-compatible builds.
- **`deepspeed` is NOT needed** for inference (it never imports at predict time) - skip it; it is the slowest/heaviest pin.
- `numpy`, `scipy`, `torch` are already in the default env - do not reinstall.

## Weights

Auto-downloaded on first run into the tFold cache (`esm_ppi_650m_tcr` = 2.43 GB, then the `tfold_tcr_trunk`).
**The default mirror is slow (~1.2 MB/s observed)** - if it crawls, the weights are also on **Zenodo** (and Google Drive), usually much faster; fetch manually into `./checkpoints`.

## Invocation

```bash
python projects/tfold_tcr/predict.py \
  --json  examples/tcr_pmhc_example.json \
  --output predictions_example/ \
  --model_version Complex
```

## Input

A TCR-pMHC complex needs five chains, each `{"id": ..., "sequence": ...}`: `A` TCR-alpha, `B` TCR-beta, `M` MHC-I heavy, `N` beta-2 microglobulin, `P` peptide.
`examples/tcr_pmhc_example.json` (bundled in the repo) is the **runnability probe input** - guaranteed well-formed, and it isolates "does the tool run on the free GPU" from input-construction risk (and does not burn free hours on a malformed hand-built file).

### Why NOT a hand-built "real patient" complex (yet)

AC2 asks for "one real patient structure", but our only local run (`results/patient_001_test/`) is **HLA-A\*02:01-negative** (alleles A\*26:01/A\*33:03/B\*15:01/B\*18:20/C\*03:02/C\*07:01).
Our pipeline's fallback TCR (DMF5) is **A\*02:01-restricted**, so pairing it with this patient's neoepitopes is a *biologically invalid* complex.
A valid real-patient complex needs an **HLA-matched TCR** for one of those alleles - a science-pairing decision (**Scientist's call**), not a Dev fabrication.
So: runnability + timing + VRAM come from the bundled example now; the biologically-real complex is a Sci-coordinated follow-up (see handoff).
A schema template for the eventual real complex is at `input/tcr_pmhc.template.json` (placeholders, do not run as-is).

## Live-run findings (2026-07-20) - COMPLETE

Machine: free **L40S** (48 GB, Ada). Operator: Jin-Ho. Runnability probe = bundled `tcr_pmhc_example.json` (PDB 6zkw), which isolates "does it run on the free GPU" from input-construction risk.

Setup notes: signup needed no card (phone verification only); env is Python 3.12.11, torch 2.8.0+cu128, CUDA True; `pip install .` into the default env + deps `termcolor biopython ml-collections dm-tree modelcif`, no deepspeed.

**Timing + VRAM (single TCR-pMHC complex):**

| Metric | CPU (free tier) | L40S (free-credit) |
|---|---|---|
| Wall-clock / complex | 4,740 s (79 min) | **8.4 s** |
| Speedup vs CPU | 1x | **~567x** |
| Peak VRAM | n/a (system RAM) | **26.0 GB** (26,001 MiB) |
| GPU utilisation | n/a | 100% during inference |
| Model load (cached weights) | - | ~16 s |

- **bf16/FP8 path engaged:** the `torch.cuda.amp.custom_fwd` autocast path fires (mixed precision on CUDA).
- **Weights cache survives a machine swap:** the CPU-run download persisted to the Studio filesystem; the L40S re-run loaded from cache (no 2.43 GB re-download), so the 8.4 s is pure inference, directly comparable to the CPU number.
- **Output:** a valid 5-chain complex PDB (`6zkw_E_D_A_B_C_Complex.pdb`).

**Cost (the $0 check):**

- Free tier = **15 credits/month (1 credit ~ $1)**; CPU Studios are **free** (the 2h46m CPU run billed $0.00). Balance at probe time: **13.69 credits**.
- The GPU probe itself cost a **fraction of one credit** (<1 min of L40S). At 8.4 s/complex, ~5 free L40S-hr/mo = **~2,100 complexes/month at $0**.
- **A GPU Studio bills per-second while ON, even idle** - unlike the free CPU box. Stop it immediately after a run. Real new money is spent only via "Add Credits"; staying within the monthly 15 is $0.

**Card-fit finding (load-bearing for the verdict):** tFold-TCR peaks at **26 GB** for one complex, so it does **NOT** fit the **24 GB L4** the Issue is premised on (nor the 16 GB T4). The minimum free card that runs it is **L40S (48 GB)** or **A100 (40 GB)**; H100/H200 also fit. Other AF3-class tools (Chai-1, Boltz-2) may fit 24 GB - that is tool-specific and untested here.

## Verdict (feeds Issue #601)

Free Lightning-AI **L40S clears the AF3-class hardware bar at a true $0**: bf16/FP8 + 48 GB + 8.4 s/complex, comfortably inside the 15-credit/mo free budget (~2,100 complexes/mo). The #601 park rationale ("wait until we can afford an L4-class Ada GPU, >=24 GB, bf16+FP8") no longer holds - that hardware is free-selectable today.

**One amendment to the park's framing:** tFold-TCR specifically needs **>=~28 GB**, so the working free card is **L40S/A100-40G, not the literal 24 GB L4** the Issue named. The L4 premise was wrong on two counts (price: $1.58 not ~$0.65/hr; VRAM: 24 GB OOMs) - but the *spirit* (a free Ada-class GPU unlocks the AF3-class frontier at $0) is **confirmed via L40S**.

**Recommendation to #601:** re-trigger the AF3-class backend re-eval - the $0 hardware constraint that forced PARK is lifted. Weight/code license was the other #601 gap; for tFold specifically it is PolyForm Noncommercial (fine for our research/portfolio use, bars commercial productization) - see the gaps table above. Note this is a *runnability* result on a curated example, not an accuracy re-eval; the accuracy/DockQ comparison and CDR3-pLDDT reranking transfer remain #601's scientific scope.

## Follow-ups

- **Biologically-real chr22 complex needs an HLA-matched TCR (Sci).** AC2's "real patient structure" clause is deferred: our only local run (`results/patient_001_test/`) is HLA-A\*02:01-negative and our fallback TCR (DMF5) is A\*02:01-restricted, so pairing them is biologically invalid. A valid complex needs a Scientist-selected HLA-matched TCR for one of that patient's alleles. Template at `input/tcr_pmhc.template.json`.
- **Accuracy re-eval + CDR3-pLDDT reranking** stay in [Issue #601](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/601)'s scientific scope (this spike only cleared the hardware gate).

## Sources

- tFold repo: https://github.com/TencentAI4S/tfold (PolyForm Noncommercial 1.0.0)
- tFold-TCR paper: Wu et al., bioRxiv 2025.01.12.632367 (+30.7% DockQ vs AlphaFold-3, >25x faster)
- Lightning-AI free tier: https://lightning.ai/docs/overview/faq/billing , https://lightning.ai/pricing/
