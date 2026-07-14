# ADR-0002: Treat cloud storage and compute differently after the GCP decommission

- **Status:** Accepted
- **Date:** 2026-07-08 (decision date; ratified on merge of [PR #1133](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/1133))
- **Deciders:** Developer, PM

## Context

The GCP stack was decommissioned on 2026-06-26 to hold a $0-budget keep-alive posture ([Issue #854](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/854)): the VMs, their disks, and the `gs://splice-neoepitope-project` bucket were deleted, with project data preserved on Cloudflare R2 first.

That single event stranded a wave of open Issues whose designs assume now-deleted infrastructure - `gs://` paths, `gsutil`, `run_cloud_gpu.sh` hooks.
And we have **no committed future cloud provider**: the RunPod migration is shelved, a return to GCP is possible, and nothing is funded.

The trap is re-scoping each stranded Issue ad hoc.
Doing so would bake fresh provider assumptions into a dozen places - assumptions that may be wrong again the next time the funding picture moves.
We need one decision, applied uniformly, rather than N pieces of guesswork.

## Decision

**Do not build a bespoke cloud-provider-agnostic framework.**
Abstracting against one dead provider and zero active ones is premature generality: by the rule of three, no second concrete implementation exists to validate the abstraction.

Instead, **split the two concerns, because they abstract very differently**:

- **Storage - neutralize cheaply.**
  GCS, S3, and R2 are approximately the same shape (a bucket, a key, a blob).
  Route every remote path through a single config key (`storage.remote_prefix`) plus a thin `rclone` wrapper.
  `rclone` already speaks all three through one interface, so **it *is* the agnostic layer** - we configure it rather than build one.

- **Compute - do NOT abstract.**
  VM provisioning, GPU drivers, and batch executors diverge sharply between providers (the P100 / closed-driver / cu126 saga is the standing evidence).
  Keep **per-provider runbooks**, extending the established pattern: [`docs/legacy/gcp_p100_setup.md`](../legacy/gcp_p100_setup.md) (archived GCP recipe) and [`docs/migration_runbook.md`](../migration_runbook.md) (RunPod forward path).
  "A funded revival is a checklist, not a reconstruction" stays the model.

**Revival trigger.** Revisit this decision only when **both** hold: (a) funding exists, **and** (b) a real second provider is in play.
The `rclone` helper itself is **deliberately not implemented now** (YAGNI) - it is a follow-up Issue activated by that trigger, not a prerequisite of this ADR.

## Consequences

- A return to GCP, or a move to any S3-compatible provider, becomes **another `rclone` remote plus the archived recipe** - no rework.
- Storage-coupled Issues re-scope to **one uniform shape** (`storage.remote_prefix` + `rclone`, local-first) instead of N bespoke rewrites.
- Compute-coupled Issues **park** under `arc:cloud-reproducibility`, blocked on funding and a provider choice, rather than being rewritten against a provider we have not chosen.
- **The reader must not hand-roll a provider abstraction.** Adding a storage backend is an `rclone` **config** change, not code. If you find yourself writing a `StorageBackend` interface, this ADR is the thing you are contradicting.
- The helper is not built yet, so a storage-coupled Issue picked up before the follow-up lands should stay local-first rather than reach for a remote path.

## Alternatives considered

- **Build a bespoke cloud-agnostic framework.**
  Rejected: premature generality. One dead provider and zero active ones give nothing to validate the abstraction against, and a wrong abstraction is more expensive than none.

- **Abstract compute as well as storage.**
  Rejected: providers diverge sharply exactly where it hurts (driver pins, GPU architecture compatibility, executor semantics). The archived GCP recipe exists *because* those details resisted abstraction the first time.

- **Do nothing; re-scope each stranded Issue ad hoc.**
  Rejected: it bakes provider assumptions into a dozen Issues and invites building against a bucket that no longer exists - which [Issue #183](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/183) nearly did.
