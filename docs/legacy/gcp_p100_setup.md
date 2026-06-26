# Legacy: GCP + P100 compute setup (archived 2026-06-26)

> **Status: ARCHIVED — not current infrastructure.**
> The GCP stack was fully decommissioned on 2026-06-26 (zero-budget cost-out; the free trial was expiring ~2026-07-03). All VMs, disks, and the GCS bucket were deleted; project data was preserved on Cloudflare R2 first (verified - [Issue #854](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/854)).
> This document preserves the GCP/P100 operational recipe so a **funded revival** is a matter of following a checklist, not reconstructing it from git history.
> The forward-looking RunPod+R2 migration plan (if reviving on a paid provider) lives in [`docs/migration_runbook.md`](../migration_runbook.md). Current posture and rationale: project memory `gcp-trial-expiry-self-funded-infra`.

## Why this was archived

With no budget for paid compute, the project moved to a `$0` keep-alive posture: the CPU-only core runs locally (chr22, MHCflurry CPU fallback), and the optional TCRdock structural validation is deferred to opportunistic free GPU (Kaggle/Colab). The GCP/P100 stack below describes infrastructure that **no longer exists**. Do not treat any VM name, zone, or bucket here as live.

## The compute stack (as it was)

- **Production VM:** `neoepitope-pipeline` — `n1-highmem-8` (8 vCPU / 52 GB RAM) + 1× NVIDIA Tesla **P100**, zone `europe-west4-a` (migrated from `europe-west1-b` after the May 2026 sustained P100 outage). The 52 GB RAM was load-bearing: OptiType peaks ~36 GB.
- **Orchestrator VM:** `neoepitope-orchestrator` — `e2-micro`, lightweight companion that started/managed the pipeline VM in detached mode and stayed up cheaply between runs.
- **GCS bucket:** `gs://splice-neoepitope-project` — results at `.../results/<patient_id>/`, logs at `.../logs/`. **Deleted 2026-06-26**; the verified copy now lives on Cloudflare R2.
- **Launch script:** `scripts/run_cloud_gpu.sh` (still in-repo). Defaults to **SPOT** provisioning + `pd-balanced` boot disk ([Issue #833](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/833)), applied only at VM **create** (provisioning model + disk type are immutable post-create). `--on-demand` forces STANDARD. On SPOT, a mid-run preemption stops the VM (`--instance-termination-action=STOP` → TERMINATED, disk + partial results preserved); re-run the same command to resume via `snakemake --rerun-incomplete`.

## NVIDIA driver pin (the P100/Pascal gymnastics)

This is the part most easily lost. P100 is Pascal (SM 6.0), which current `-open` driver variants do **not** support (they need GSP firmware Pascal lacks → `NVRM: ... not supported by open nvidia.ko because it does not include the required GPU System Processor`).

- **Pin the image** `pytorch-latest-cu124-v20250327-ubuntu-2204` (deprecated but READY) **+ `install-nvidia-driver=True` metadata.** The DLVM `install-driver.sh` runs on first boot and installs the proprietary closed driver **550.90.07** via DKMS (Pascal-compatible). It branches on machine type: `machine_type =~ ^n` → `open_kernel_module_arg=""` (proprietary), required for P100.
- **Do NOT** switch to current family aliases (`common-cu129-*`, `pytorch-2-7-cu128-*`): every fresh DLVM image as of 2026-05-28 ships `-open` driver variants (580-open + 570-open both verified to fail on P100 with the GSP error).
- **Do NOT** install drivers via apt (`nvidia-headless-570-server` is a hard-`Depends:` wrapper for 580; bare `nvidia-utils-535-server` churns userspace past the image's kernel-module version and DKMS can't rebuild).
- `run_cloud_gpu.sh` verifies via an `nvidia-smi` smoke test with a 3-min retry for the first-boot install window; subsequent boots persist via DKMS ([Issue #522](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/522)).

## PyTorch / `python.yaml` — P100 (Pascal SM 6.0) torch pin

PyTorch wheels built against CUDA 12.8/12.9 dropped SM 6.0 (Pascal). On a P100, `torch.cuda.is_available()` still returns `True` but kernel dispatch fails. The fix was pinning a Pascal-compatible CUDA 12.6 build:

```yaml
# the historical pin (Linux/x86 + P100 only)
- --extra-index-url https://download.pytorch.org/whl/cu126
- torch ==2.12.0+cu126   # Pascal-compatible; the +cu126 local-version tag forces the cu126 wheel
```

The `+cu126` local-version tag forced resolution to the Pascal-compatible build; PyPI's same-numbered wheel lacks the tag and sorts lower per PEP 440. Empirically verified on `neoepitope-pipeline` 2026-05-27 ([Issue #352](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/352)): `torch==2.12.0+cu126` dispatches ReLU on a P100 cleanly.

**Current state (post-archive):** `workflow/envs/python.yaml` now selects torch by platform via PEP 508 markers — the `+cu126` pin is gated to `platform_system == "Linux"` (preserved for this revival path), and macOS arm64 gets a plain `torch >=2.12` wheel so the env builds for local CPU runs. To revive a P100 host, the Linux marker already resolves to the cu126 pin; no change needed. CUDA 12.6 wheels remain Pascal-compatible through `torch 2.12`; when the cu126 channel falls behind a torch version that's needed, pin to the last cu126 build or migrate off P100.

## Historical P100 unavailability — `europe-west1-b` (pre-migration)

P100 capacity in the original zone `europe-west1-b` went sustained-exhausted multiple times before the VM was migrated to `europe-west4-a`:

- 2026-05-06: 11 launch attempts over ~7h16m, all `ZONE_RESOURCE_POOL_EXHAUSTED` on `n1-highmem-8 + nvidia-tesla-p100` (10:32 → 17:48 BST).
- 2026-05-08: capacity probe still exhausted at ~16:00 UTC; the rolling outage spanned ~46h across 05-06 → 05-08.

The migration to `europe-west4-a` resolved the immediate pain. The longer-term hardware-diversity strategy (T4/L4 hybrid fallback) was tracked in [Issue #310](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/310). If reviving and a zone exhausts, the historic mitigation was to retry overnight — capacity typically recovered after a sustained gap.

## Reviving (checklist)

1. Recreate the GCS bucket (or point storage at R2) and stage references/indices.
2. `scripts/run_cloud_gpu.sh` to create the VM (pick `--on-demand` vs SPOT; both immutable post-create).
3. Confirm the driver pin above lands (`nvidia-smi` smoke test) — the closed 550.90.07 via the pinned DLVM image.
4. `python.yaml`'s Linux marker already resolves to the cu126 Pascal-safe torch.
5. See `docs/migration_runbook.md` if migrating to a modern paid GPU (RunPod L4/Ada) instead — Ada deletes the entire Pascal driver-pin saga.
