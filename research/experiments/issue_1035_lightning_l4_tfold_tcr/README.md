# Issue 1035 spike: free Lightning-AI GPU as a $0 path to AF3-class TCR-pMHC co-folding

Turnkey recipe for running **tFold-TCR** on a **free Lightning-AI GPU** at a hard $0 ceiling.
This is the desk-prep + live-setup record; the prediction timing/VRAM numbers are captured on the live run (see [Live-run findings](#live-run-findings-2026-07-20)).

Tracks [Issue #1035](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1035); feeds the backend re-eval parked in [Issue #601](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/601).

## Freshness check (2026-07-20) - premise mostly holds, with one correction

Lightning-AI free tier, per its [billing FAQ](https://lightning.ai/docs/overview/faq/billing) + confirmed live at signup:

- **15 credits/month, no credit card required** - but **phone verification IS required** (so "$0" is not zero-friction).
- **No L4 in the current free lineup.** The offered machines + free-hour budgets are: CPU (free), **T4 (79 hr)**, **L40S (5 hr)**, **A100 (3 hr)**, **H100 (3 hr, no wait)**, **H200 (2 hr, 4+ min wait)**.
  - This *corrects* the Issue's filed premise (L4 at ~21-31 hr/mo) and its assumption that A100 is paid-gated: **A100/H100/H200 are all free-tier selectable here.**
  - The right Ada-class pick is **L40S** (48 GB, bf16 + FP8, the most free hours of the capable cards). T4 is 16 GB Turing with no bf16/FP8 - the same wall Kaggle/Colab hit, so avoid it.
- Free hours are **per-GPU** (the advertised "up to 80" is essentially the T4 figure). tFold-TCR runs in seconds/complex, so L40S's ~5 hr is ample for a probe.

## The two open gaps #601 flagged

| Gap | Status | Finding |
|-----|--------|---------|
| Weight/code license | **CLOSED** | Repo is [PolyForm Noncommercial 1.0.0](https://github.com/TencentAI4S/tfold/blob/master/LICENSE) - same terms for code and weights. Fine for our research/portfolio use; bars any commercial productization. |
| VRAM / bf16-FP8 wall | **No platform wall** | The free L40S runs Python 3.12 + **torch 2.8.0+cu128, CUDA available** - bf16/FP8-capable. Peak-VRAM per complex is the last live number (pending). |

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

## Live-run findings (2026-07-20)

Machine: free **L40S** (48 GB, Ada). Operator: Jin-Ho.

- Signup: no card, phone verification required. GPU lineup captured above (no L4; A100/H100/H200 free-selectable).
- Env: Python 3.12.11, torch 2.8.0+cu128, CUDA True. `pip install .` succeeds in the default env; deps `termcolor biopython ml-collections dm-tree modelcif` added; no deepspeed.
- tFold-TCR imports clean and begins ESM-650M weight download (2.43 GB) on first run.
- **Pending:** per-complex wall-clock, peak VRAM, bf16/FP8 path engaged. Fills in after the download completes and prediction runs.

## Operator handoff (remaining)

1. Let the weight download finish, then let `predict.py` run the bundled example to completion.
2. Capture per-complex wall-clock + peak VRAM (`nvidia-smi`) + any bf16/FP8 log line. Paste back to Developer.
3. **Stop the Studio** when done (free L40S budget is ~5 hr and ticks while the Studio is on).
4. Developer lands the timing/VRAM table (AC3) and posts the verdict to [Issue #601](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/601) (AC4): does free L40S clear its HW park, and does it move the P2/P3 call?
5. Developer flags the "biologically-real chr22 complex needs an HLA-matched TCR" piece to Scientist as a follow-up.

## Sources

- tFold repo: https://github.com/TencentAI4S/tfold (PolyForm Noncommercial 1.0.0)
- tFold-TCR paper: Wu et al., bioRxiv 2025.01.12.632367 (+30.7% DockQ vs AlphaFold-3, >25x faster)
- Lightning-AI free tier: https://lightning.ai/docs/overview/faq/billing , https://lightning.ai/pricing/
