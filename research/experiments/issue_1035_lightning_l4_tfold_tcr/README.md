# Issue 1035 spike: free Lightning-AI L4 as a $0 path to AF3-class TCR-pMHC co-folding

Turnkey recipe for running **tFold-TCR** on a **free Lightning-AI L4** at a hard $0 ceiling.
This is the desk-prep half of the spike (assembled without a GPU); the live run is a signup + paste-and-run handoff (see [Operator handoff](#operator-handoff)).

Tracks [Issue #1035](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1035); feeds the backend re-eval parked in [Issue #601](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/601).

## Freshness check (2026-07-20) - premise holds

Lightning-AI free tier, per its [billing FAQ](https://lightning.ai/docs/overview/faq/billing) + [pricing](https://lightning.ai/pricing/):

- **15 credits/month (1 credit = $1), no credit card required**, renewing monthly.
- Free interruptible Studios expose **L4 (Ada, 24 GB, bf16 + FP8)** at ~$0.48-0.70/hr, i.e. **~21-31 free L4-hours/month at $0**.
- Free Studios restart every ~4 hours (interruptible), which suits tFold-TCR's seconds-per-complex batch profile.

The one thing web search cannot confirm is whether the current UI exposes L4 selection with **no card on file** - that is a signup-time check and the **first live AC**.

## The two open gaps #601 flagged

| Gap | Status | Finding |
|-----|--------|---------|
| Weight/code license | **CLOSED** | Repo is [PolyForm Noncommercial 1.0.0](https://github.com/TencentAI4S/tfold/blob/master/LICENSE) - same terms for code and weights. Fine for our research/portfolio use; bars any commercial productization. |
| VRAM / bf16-FP8 wall | **Bounded, resolved by the live probe** | Undocumented in the repo. Architecture is an ESM-650M PPI module + a folding trunk, "seconds per complex", so it is expected to fit the L4's 24 GB comfortably. The peak-VRAM + precision-engaged measurement IS the live probe (AC2). |

## Environment (from the repo `environment.yaml`)

```bash
git clone https://github.com/TencentAI4S/tfold.git
cd tfold
conda env create -f environment.yaml   # name: tfold, python=3.8
conda activate tfold
pip install .                          # or: pip install tfold
```

Key pins: `python=3.8`, `torch>=2.0.0`, `deepspeed==0.12.3`, `biopython==1.79`, `ml-collections==0.1.1`, `dm-tree==0.1.8`, `numpy==1.21.2`, `modelcif==0.9`, `scipy`.
`mmseqs2` is only needed for tFold-Ag (antigen), not tFold-TCR - skip it here.

## Weights

Auto-downloaded on first run, or fetch manually into `./checkpoints` from any of the three mirrors (Tencent Weiyun / Google Drive / Zenodo).
tFold-TCR pulls `esm_ppi_650m_tcr` + `tfold_tcr_trunk`.

## Invocation

```bash
python projects/tfold_tcr/predict.py \
  --json  input/tcr_pmhc.json \
  --output predictions/ \
  --model_version Complex
```

## Input format (5 chains, JSON)

A TCR-pMHC complex needs five chain entries, each `{"id": ..., "sequence": ...}`:

| id | chain |
|----|-------|
| `A` | TCR alpha |
| `B` | TCR beta |
| `M` | MHC class I heavy chain |
| `N` | beta-2 microglobulin |
| `P` | peptide |

A template lives at `input/tcr_pmhc.template.json` (placeholder sequences, clearly marked - do not run as-is).

## Two probe inputs, in order

1. **Runnability smoke** - run tFold-TCR's own bundled example first (`examples/tcr_pmhc_example.json` in the repo).
   This isolates "does the tool run on the free L4" from "did I build the input correctly".
2. **Real-patient structure** (satisfies AC2's "one real patient structure") - build the 5-chain JSON from *our* pipeline:
   - **peptide (P)**: a real chr22 neoepitope from `results/local/mhcflurry/presentation.tsv` (a `HLA-A*02:01` presenter).
   - **MHC heavy (M)**: HLA-A\*02:01 mature heavy-chain sequence (IMGT/HLA or UniProt `P01892`).
   - **beta-2 microglobulin (N)**: human B2M mature sequence (UniProt `P61769`).
   - **TCR alpha/beta (A/B)**: our DMF5 fallback TCR (the HLA-A\*02:01-restricted default in `run_tcrdock.py` / `gpu_config.yaml`): TRAV12-2 / TRAJ21, CDR3a `CAVNFGGGKLI`; TRBV6-5 / TRBJ2-7, CDR3b `CASSLAGGRPEQYF`. Stitch the gene+CDR3 into full alpha/beta sequences with `stitchr` (already used elsewhere in the pipeline).

Do not fabricate any of these sequences - source each from IMGT/UniProt/stitchr so the probe is a real structure, not a plausible-looking one.

## What to record (feeds AC2 + AC3)

For each run: wall-clock per complex, peak VRAM (`nvidia-smi`), whether the bf16/FP8 path engages, and the observed credit -> L4-hour burn rate (feeds AC1's estimate check).

## Operator handoff

The two live ACs need a Lightning-AI account + an interactive Studio, which the agent cannot create or drive. Operator steps:

1. Sign up for the Lightning-AI free tier (no card); confirm **L4 is selectable at $0** and note the credit->hour rate. (AC1)
2. In a free L4 Studio, run the Environment + Weights + Invocation blocks above on the bundled example, then on the real-patient input. Capture the metrics above. (AC2)
3. Paste the metrics back here (or to the Developer) - the Developer lands the timing/VRAM table in the lab notebook (AC3) and posts the verdict to [Issue #601](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/601): does free L4 clear its hardware park, and does it move the P2/P3 call on the AF3-class re-eval? (AC4)

## Sources

- tFold repo: https://github.com/TencentAI4S/tfold (PolyForm Noncommercial 1.0.0)
- tFold-TCR paper: Wu et al., bioRxiv 2025.01.12.632367 (+30.7% DockQ vs AlphaFold-3, >25x faster)
- Lightning-AI free tier: https://lightning.ai/docs/overview/faq/billing , https://lightning.ai/pricing/
