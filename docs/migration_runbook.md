# Migration Runbook — GCP → RunPod + Cloudflare R2

> **Status:** plan verified & green-lit (2026-06-22), **not yet executed**. Account signup + cutover are operator-gated.
> **Why now:** the GCP free 90-day $300 trial expires **~2026-07-03**. Infra is self-funded (no academic affiliation → no academic-credit/HPC routes). Usage is **bursty** (per-patient runs, VMs idle between), so the right model is **per-second GPU billing + zero-egress storage**, not an always-allocated VM.

Companion bridge work on the existing GCP path (Spot toggle + `pd-balanced` disk) shipped in [#833](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/833) / [PR #834](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/834) — it trims burn on the remaining runway but does **not** replace the migration (provisioning model is fixed at VM create; the real lever is leaving GCP).

---

## 1. Target stack (verified)

| Layer | Choice | Key numbers |
|-------|--------|-------------|
| **Compute** | **RunPod GPU Pod — NVIDIA L4, On-Demand Secure Cloud** | **$0.39/hr**, 24 GB VRAM, **50 GB RAM**, 12 vCPU, Ada arch |
| **Static storage** | RunPod **Network Volume** (references + indices + conda envs + AlphaFold params) | **$0.07/GB/mo** (~$14 for 200 GB), survives stop/terminate/preemption, **DC-pinned** |
| **Bulk / results / cold backup** | **Cloudflare R2** | **$0.015/GB/mo** (~$3 for 200 GB), **zero egress** |

**Why L4 is the pick (not the cheaper RTX 4090):**
- **The #1 migration risk is retired by RAM, not VRAM.** OptiType (HLA typing) peaks at **~36 GB RAM** — the reason the GCP VM is `n1-highmem-8`/52 GB. A naive "cheap GPU pod" often bundles only 24–32 GB RAM and would OOM. **L4 bundles 50 GB RAM** (14 GB headroom) and 12 vCPU. RTX 4090 also clears it (41 GB / 6 vCPU) but with thin headroom and fewer cores, and on Secure Cloud the 4090 is *more* expensive ($0.69) than the L4.
- **24 GB VRAM is safe (not borderline) for AlphaFold/TCRdock.** The ~600–800-residue TCR-pMHC complex is roughly *half* the 16 GB-VRAM ceiling and already runs on the **16 GB P100** — 24 GB is +50% headroom. Ada (CC 8.9) also unlocks bf16 + Pallas (~2× faster than Pascal).
- **Use On-Demand Secure, not Spot/Community,** for patient runs — a ~5-second-warning eviction mid-AlphaFold isn't worth the ~$0.15/hr discount for clinical determinism. Reserve Spot for non-urgent batch reprocessing.

**This migration deletes the P100/Pascal driver-pin saga.** On Ada, current CUDA/torch wheels just work — the CLAUDE.md `python.yaml` `+cu126` pin, the `IMAGE_NAME` DLVM image pin, and the closed-550 driver gymnastics ([#522](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/522)) all become moot on the new stack.

## 2. Cost

| Stack | Compute (20 GPU-hr/mo) | Storage | Egress | **Monthly** |
|-------|------------------------|---------|--------|-------------|
| **RunPod L4 + hybrid (vol + R2)** | 20 × $0.39 = $7.80 | $14 (vol) + $3 (R2) | $0 | **~$25** |
| RunPod L4 + pure-R2-pull | $7.80 | $3 (R2) | $0 | **~$11** (re-downloads 200 GB each run) |
| **GCP on-demand-equivalent** | ~$38 (P100 + n1-highmem-8) | $34 pd-ssd (billed while stopped!) + $4 GCS | $0.12/GB | **~$75–80** |

The hybrid's $14 Network Volume buys instant bursts (no per-run re-download) and preemption safety; the static core stays mounted. Pure-R2 is cheaper but adds a 200 GB pull to every run start. **Recommendation: hybrid.**

> The GCP row shows the **pre-#834** `pd-ssd` steady-state being *replaced* (the historical baseline). Post-bridge, #834's `pd-balanced` cuts that disk line to ~$20/mo for 200 GB — still above RunPod's ~$17 hybrid storage, and the always-allocated-while-stopped problem persists either way; the migration removes it.

## 3. Phased checklist (~11 days)

### Days 1–2 — Data exit (hard-deadline item; do first)
- [ ] Create a Cloudflare R2 bucket + an S3-API token (note the `<ACCOUNT_ID>.r2.cloudflarestorage.com` endpoint + keys).
- [ ] On the **GCP VM** (so GCS→internet egress is the only hop), install `rclone` (**≥ v1.59** — older versions 401 against R2).
- [ ] Configure remotes and copy (see §4). ~115 GB ≈ **~$14** GCS egress.
- [ ] `rclone check` to verify checksums; spot-open a few references + indices.

### Days 3–5 — Compute port
- [ ] RunPod account + payment. Pick **L4 / Secure**. Verify the deploy-time rate (the two RunPod pricing pages disagree on L4: $0.39 vs $0.44 = Secure vs Community).
- [ ] **Slim the Docker image:** move AlphaFold2 params + BLAST DBs *out* of the 25 GB image onto the Network Volume; rebuild image as code + CUDA only (a few GB). RunPod's own example drops a 16 GB-weights cold start from ~3 min → ~20 s this way.
- [ ] Rebuild the image on **current CUDA/torch** for Ada (drop the Pascal pins). Keep the TCRdock container's internal CUDA 11.8 as-is (self-contained).
- [ ] Push to a registry RunPod can pull (GHCR / Docker Hub); create a RunPod **Pod template** (image, `openssh-server`, container-disk sized **above** the unpacked image footprint, Network Volume mounted at `/workspace`).

### Days 6–9 — Smoke test (the integration-run discipline)
- [ ] Launch a pod; `nvidia-smi`; confirm `torch.cuda.is_available()` **and** a real kernel dispatch on Ada.
- [ ] Run **chr22 end-to-end** (`config/test_config.yaml`) in tmux with `--use-conda` — catches the conda-solver / subprocess-CLI / NaN classes that dry-run misses.
- [ ] Run the MHCflurry `Class1PresentationPredictor.predict()` GPU path on the test patient (the `_has_gpu()` smoke test) + one TCRdock/AlphaFold structure run.
- [ ] Time a representative run → confirm the GPU-hr budget estimate.

### Days 10–11 — Cutover + decommission
- [ ] Point the orchestration + docs at RunPod + R2; branch or retire the GCP-specific zone/driver-pin logic in `run_cloud_gpu.sh`.
- [ ] **Delete the `neoepitope-orchestrator` e2-micro too.** It exists only to start/manage the pipeline VM in GCP detached mode — RunPod's per-second billing + API makes that pattern obsolete, so the orchestrator *retires* (not migrates). Easy to leave it quietly billing otherwise.
- [ ] Final `rclone check` GCS↔R2. **Tear down GCP VMs + delete the bucket before the trial suspends** (a suspended-but-not-deleted bucket can still incur charges once billing flips).
- [ ] Keep one cold backup until ≥1 real patient run succeeds on the new stack.

## 4. Exact commands

### rclone GCS → R2
```bash
# On the GCP VM. rclone >= 1.59 required.
rclone config
#  remote 'gcs': type=google cloud storage, service_account_file=<sa.json>
#  remote 'r2' : type=s3, provider=Cloudflare, region=auto,
#                endpoint=https://<ACCOUNT_ID>.r2.cloudflarestorage.com,
#                access_key_id=<R2_KEY>, secret_access_key=<R2_SECRET>

# Copy. --s3-chunk-size 64M is REQUIRED: a 200 GB object at the 5 MiB
# default part size = ~40,000 parts > R2's 10,000-part multipart cap.
# 64 MiB × 10,000 = ~625 GiB ceiling — comfortable.
rclone copy gcs:splice-neoepitope-project r2:splice-neoepitope-project \
  --s3-chunk-size 64M --transfers 16 --checkers 32 --progress

# For object-scoped R2 tokens, add: --s3-no-check-bucket
# Keep checksums ON for reference data (do NOT pass --s3-disable-checksum):
# multipart ETag != MD5, so rclone stores X-Amz-Meta-Md5chksum for integrity.

rclone check gcs:splice-neoepitope-project r2:splice-neoepitope-project --one-way
```

### Reference fetch at run start (replaces `gs://` paths)
```bash
# On the pod, if NOT keeping the static core on the Network Volume:
rclone copy r2:splice-neoepitope-project/references ./references --s3-chunk-size 64M --transfers 16
# With the hybrid design, references live on the mounted Network Volume — no pull needed.
```

### AlphaFold params on the volume (cold-start fix)
```bash
# One-time: stage params + BLAST DBs onto the Network Volume, not the image.
# Pod template mounts the volume at /workspace; point TCRdock/AlphaFold at:
#   --data_dir=/workspace/alphafold_params   (or the equivalent env/flag)
# Image then carries only code + CUDA (a few GB) → fast cold pulls.
```

## 5. Verify-in-console before committing
- **L4 capacity in the datacenter you pin the Network Volume to.** A volume is region/DC-locked, and attaching it constrains the GPU pool to that DC — the direct analog of the GCP P100 zone-exhaustion pain. Check L4 availability in-region *before* creating the volume; provision a volume per region if needed.
- **L4 price** ($0.39 vs $0.44) at your chosen tier/region.
- **25 GB cold-pull time** on real RunPod hosts (docs say "a few minutes," no SLA) — benchmark once, then decide how hard to slim/offload. Note: FlashBoot image caching is **Serverless-only**; Pods land on new hosts each burst, so plan for a cold pull most runs unless the image is slim.
- **Community Cloud + Network Volume** support (volumes are best-documented on Secure) — only relevant if you ever choose Community.

## 6. Excluded / fallback
- **Vast.ai** — cost-floor plan B (4090 interruptible ~$0.30/hr) but host-set egress, **bills on *stopped* instances (must *delete* to stop charges)**, and marketplace reliability tax. Not plan A for clinical work.
- **Kaggle** free P100 (~30 GPU-hr/wk) — notebook-only, no root/Docker/Snakemake → **experiment notebooks only**, never the production pipeline.
- **Lambda / Paperspace / Modal** — no per-second-idle advantage (Lambda), weaker bursty economics (Paperspace), or serverless-function mismatch (Modal).
- **Phantom prices:** RunPod T4/A10 are **not** in the current (2026-06-22) catalog — earlier aggregator figures were dropped; L4 is the verified low tier.

## 7. Provenance
Two research passes (2026-06-22), all core figures from official `runpod.io` / `developers.cloudflare.com`; aggregator-only numbers dropped. Full findings + dated sources are in the Developer project memory (`gcp-trial-expiry-self-funded-infra`). Burn audit + bridge-PR context: [#833](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/833).
