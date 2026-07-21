# AF3-class structure-backend re-eval - updated verdict (Issue #601 AC2)

**Date:** 2026-07-21
**Author:** Scientist
**Issue:** [#601](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/601) AC2 (re-run the [#316](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/316) survey verdict against the new hardware envelope), arc `scoring-tcr-pmhc`, best-bet 7 of the open-gap program ([#678](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/678)).
**Scope:** desk re-run of the #316 verdict. No GPU runs here; the empirical DockQ benchmark this verdict *recommends* is a downstream deliverable, not this one.

## TL;DR

The park was **hardware-gated**, and that gate is now down (Dev's [#1035](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1035) spike: free Lightning L40S 48 GB / A100 40 GB run tFold-TCR at 8.4 s/complex, 26 GB peak, $0). **But lifting the hardware gate does not by itself settle "migrate TCRdock -> an AF3-class backend, and to which one."** The accuracy evidence is genuinely mixed, and the two variables that actually decide the migration - *which* open AF3-class backend, and whether it beats our AF2 baseline on **our** kind of complex - are still untested by us.

**Verdict: UNPARK to a bounded empirical comparison, not to a migration.** Recommend a small, public-data DockQ benchmark of the runnable open AF3-class candidates against our current TCRdock/AF2 baseline *before* committing any backend swap. Keep TCRdock as the production backend until that benchmark reports.

## What changed since #316

| Axis | #316 verdict (parked) | Now (2026-07-21) |
|---|---|---|
| **Hardware** | AF3-class blocked on P100 (no bf16, 16 GB < 24 GB) | Lifted. Free L40S (48 GB) / A100 (40 GB) run tFold-TCR at $0 ([#1035](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1035)). Real limiter is **VRAM ~40 GB**, not FP8: tFold-TCR's 26 GB peak OOMs a 24 GB L4. |
| **License** | AF3 non-redistributable; Chai-1/Protenix Apache; ESMFold2 MIT | Unchanged and re-confirmed: **AF3 weights non-redistributable** (manual-grant, non-commercial); **Boltz-2 MIT for code AND weights** (commercial-OK); **Chai-1 Apache-2.0 for code + weights** (since Nov 2024); **tFold-TCR PolyForm-Noncommercial** (research/portfolio OK, bars commercial productization). |
| **Accuracy** | AF3 best DockQ (Lu et al.); AF2 lags | Still AF3-first, but the *open* picture is nuanced (see below). |

## Accuracy landscape (from our shelf + verified web)

**Verified against our Zotero shelf (folder 4, TCR-pMHC structure & binding):**

- **Lu et al. unified benchmark** (`3GS2FXXZ`, Briefings in Bioinformatics 2026, DOI `10.1093/bib/bbag289`; 10 methods, 70 unseen complexes, DockQ + interface): **AF3 leads, AF2 lags.** CDR3-pLDDT reranking adds +4.3% Top-1 and captures >80% affinity-impact mutations - but that signal is validated on **AF3 only** and (per #601's own note) does **not** transfer to AF2/TCRdock. tFold-TCR is one of the benchmarked methods.
- **JCIM comparative** (`Z5AP3IT3`, 2025, DOI `10.1021/acs.jcim.5c00298`; 7 isolated-TCR + 6 TCR-pMHC tools, 40 alphabeta TCRs + 27 complexes): tFold-TCR + AF2/AF3 all excel on *overall TCR structure*; CDR3 loops + docking orientations remain hard for all; explicitly **"supports TCRdock retention where overall accuracy matters"** and flags a **"backbone concern for rescoring atop tFold-TCR."**

**From verified web (July 2026), corroborating the license axis and adding the open-model ranking:**

- Among AF3 replications, **Boltz-1/Boltz-2 track AF3 most closely**; Chai-1-MSA is comparable to AF3 on class II but underperforms on class I. MSA-based methods (esp. AF3) outperform PLM-based and docking-based approaches on TCR-pMHC.
- A crystal-structure benchmark signal (via web summary) puts **tFold-TCR and Boltz-2 at acceptable-quality-only (no medium/high-quality) models** on TCR-pMHC. *(Directional only, NOT load-bearing: I could not fetch a primary this session - biorxiv 403'd, the Shi thesis 405'd, the JCIM version is paywalled - and the pre/post-training framing overlaps our already-shelved `Z5AP3IT3`, so this may be a re-description of that study at a different metric/N rather than a distinct source. Not cited as a number; the verdict below stands on the two verified shelf notes alone.)*

## Why the hardware fix does not settle the migration

1. **The strongest model is still license-blocked.** AF3 wins DockQ but its weights are non-redistributable - it cannot ship in our Docker image regardless of hardware. So "AF3 leads" does not translate into a shippable backend.
2. **The runnable-at-$0 backend is not the accuracy winner.** Dev spiked **tFold-TCR** (because it's TCR-specific and packaged), but the shelf evidence has tFold-TCR at parity-or-below AF2/TCRdock for the thing we care about, and one benchmark puts it at acceptable-quality-only. Adopting it because it *runs* free would be optimizing for runnability, not accuracy.
3. **The best-licensed open options are untested by us on both axes.** Boltz-2 (MIT) and Chai-1 (Apache) are the models that actually track AF3 - but (a) their TCR-pMHC-*specific* standing is weaker/mixed in the crystal benchmarks, and (b) their **free-tier runnability is unverified** (Dev: "Chai-1/Boltz-2 may fit 24 GB; untested"). We have a runnability result for the *wrong* model relative to the accuracy evidence.
4. **CDR3-pLDDT reranking (AC3) is AF3-class-coupled and does not transfer to AF2.** So its value only materializes *if* we actually move to an AF3-class backend - which circles back to picking one that both ships (license) and wins (accuracy).

## Verdict and recommended next step

**UNPARK #601; do NOT migrate on desk evidence alone.** The desk verdict is that the migration decision now hinges on an empirical question we can answer cheaply on public data:

**Recommended next deliverable (new AC / issue, the elevated AC1 "optional spike"):** a bounded DockQ benchmark on the **Lu et al. public complex set** (or a subset), comparing:
- our **TCRdock / AF2** (production baseline),
- **tFold-TCR** (already runnable per #1035),
- **Boltz-2** (MIT) and **Chai-1** (Apache) - *if* they clear the free-tier VRAM/precision walls (runnability sub-check first),

on **DockQ + interface metrics**, no patient data. Decision rule: migrate only to an open backend that (a) ships (MIT/Apache preferred; PolyForm-NC acceptable for research/portfolio) **and** (b) beats TCRdock/AF2 on DockQ by a margin worth the integration cost. If none does, TCRdock retention is the correct call and #601 re-parks on the AC3/AC4 legs.

**AC dispositions after this verdict:**
- **AC1** (GPU target): discharged by #1035 (working cards L40S/A100; real floor ~40 GB VRAM, not FP8).
- **AC2** (this verdict): delivered. Verdict = unpark-to-benchmark, above.
- **AC3** (deploy CDR3-pLDDT column + reranking): stays gated on an actual AF3-class backend move - which the recommended benchmark decides.
- **AC4** (patient-specific spike): carved to [#1245](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1245) (needs HLA-matched TCR), blocked, unchanged.

## Prior-art / curation

1. **"A Unified Framework for TCR-pMHC Structural Model Assessment"** (biorxiv `2025.10.09.681411`) - a **reference-free QA classifier** (RF over pLDDT/ipTM/iPAE/pDockQ, trained on 1160 models of 232 PDB class-I complexes; used AF3 to build a 33,820-complex synthetic set). Directly relevant to AC3-style confidence-based reranking/filtering. **ADDED to Zotero this session** (key `E8XR8267`, folder 4 TCR-pMHC structure & binding).
2. **The crystal-structure "acceptable-quality-only" benchmark signal** - primary unreachable via open routes this session (biorxiv 403, Shi thesis 405, JCIM paywalled). It likely overlaps our already-shelved `Z5AP3IT3` (the Shi JCIM study) rather than being a distinct paper; not added, and not cited as a number. Follow-up if it turns out distinct: pin the primary via Europe PMC / author request before shelving and citing.

## Honest scope

This is a **desk verdict**, not an accuracy measurement. It re-runs the #316 decision against the new hardware + license + published-accuracy picture and concludes the migration is now an empirical question worth one bounded benchmark - it does not itself measure DockQ on any backend. All accuracy claims trace to our shelved benchmarks (`3GS2FXXZ`, `Z5AP3IT3`) or verified web; the one un-pinned attribution is flagged inline.
