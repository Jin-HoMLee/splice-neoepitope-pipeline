# Lightning-AI free-GPU runner for tFold-TCR

Hands-free driver for running **tFold-TCR** structure prediction on a free Lightning-AI
Studio, from your laptop. Replaces the manual screenshot loop from the
[Issue #1035](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1035) spike.
Built in [Issue #1249](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1249).

## Setup (one time)

```bash
uv venv --seed tools/lightning/.venv
uv pip install --python tools/lightning/.venv/bin/python -r tools/lightning/requirements.txt
tools/lightning/.venv/bin/lightning login          # browser sign-in; stores creds in ~/.lightning/credentials.json
```

`lightning login` writes a matched `user_id` + `api_key` to `~/.lightning/credentials.json`,
which the SDK reads automatically. No secrets go in the repo or in `.env`. The `.venv` is
gitignored and per-clone.

## Usage

```bash
PY=tools/lightning/.venv/bin/python
$PY tools/lightning/run_tfold.py status     # status + Drive usage, no compute
$PY tools/lightning/run_tfold.py smoke      # $0 CPU dry-run of the whole loop
$PY tools/lightning/run_tfold.py prestage   # download the 2.43 GB weights to persistent cache (free CPU, ~20 min once)
$PY tools/lightning/run_tfold.py predict --input complex.json --gpu   # real run on the L40S (~8 s, billed in seconds)
```

`predict` defaults to **CPU** (free but ~79 min/complex); `--gpu` uses the L40S (~8 s but
billed per second). It always stops the Studio on exit unless `--keep-alive`.

## Why the loop is shaped this way (all verified live)

- **Weights persistence.** The 2.43 GB ESM-650M weights download via `torch.hub.get_dir()`
  (`tfold/model/pretrain.py`), which honors **`TORCH_HOME`**. By default they land in an
  ephemeral path and a Studio **stop wipes them** (a machine *switch* keeps them). The runner
  sets `TORCH_HOME=$HOME/.torch` - `$HOME` is the persistent mount
  (`/teamspace/studios/this_studio`) - so once `prestage` runs, the weights survive stops and
  every later run skips the 20-min re-download. No patch to tfold, just the env var.
- **CPU is free, the L40S bills per second.** Downloads run on CPU; the L40S is entered only
  for the seconds of inference. Never leave a GPU idling (it bills while on, even idle).
- **Teamspace is org-owned.** `Studio("pipeline-devbox", teamspace="pipeline-automation-project", org="jh-m-lee-lab")`
  - the `org=` kwarg is required. The CLI's `owner/teamspace` form mis-resolves as
  user-owned in the SDK and raises "Teamspace ... does not exist".
- **Storage cost.** First **10 GB** of persistent Drive is free, then $0.10/GB/mo. The weight
  cache is ~2.5 GB, well under the cap. With **no payment method on file, Lightning cannot
  charge** - it blocks at the limit instead of billing. The runner prints Drive usage vs the
  10 GB cap on every run so we never drift over it silently.

## Safety rails in the runner

- Refuses to start a paid GPU unless `--gpu` is passed, and hard-aborts (stopping the Studio)
  if a start resolves to a non-CPU machine unexpectedly.
- Stops the Studio on exit (including on error) unless `--keep-alive`.
- Never prints credentials.

## Scope

This drives the **optional structural-validation** step only (a free-GPU opportunistic run),
not any core pipeline stage. TCR selection and structure interpretation are the Scientist's
call ([Issue #1245](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1245)); this
is only the instrument.
