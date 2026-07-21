#!/usr/bin/env python3
"""Drive a free-GPU tFold-TCR run on a Lightning Studio, hands-free (Issue #1249).

Replaces the manual screenshot loop from the Issue #1035 spike. Auth is read from
~/.lightning/credentials.json (run `lightning login` once); this script never prints
secrets.

The loop, and why each step is shaped this way (all verified live, see tools/lightning/README.md):
  - The 2.43 GB ESM-650M weights download to torch.hub's dir, which honors TORCH_HOME.
    Pointing TORCH_HOME under $HOME (the persistent mount) makes them survive a stop.
  - CPU is free; the L40S bills per second. So downloads run on CPU, and the L40S is
    only entered for the seconds of inference.
  - A Studio `stop` wipes ephemeral paths but keeps $HOME. A machine `switch` keeps
    everything. We never leave a GPU running.

Subcommands:
  status    Report Studio status, machine, and Drive usage vs the 10 GB free cap. No compute.
  smoke     $0 CPU dry-run of the whole loop (start CPU, echo, upload/download round-trip, stop).
  prestage  Download the tFold weights to the persistent TORCH_HOME on free CPU, then stop.
  predict   Upload an input complex, run predict.py, download the output PDB, stop.
            Defaults to CPU (free but ~79 min/complex); pass --gpu for the L40S (~8 s, billed).
"""
import argparse
import os
import sys
import tempfile
from pathlib import Path

from lightning_sdk import Machine, Studio

# --- Target Studio (verified live). The teamspace is ORG-owned: org= is required; the
# --- CLI's owner/teamspace form mis-resolves as user-owned in the SDK.
STUDIO_NAME = "pipeline-devbox"
TEAMSPACE = "pipeline-automation-project"
ORG = "jh-m-lee-lab"

# Persistent weight cache: under $HOME (/teamspace/studios/this_studio), which survives a stop.
REMOTE_TORCH_HOME = "$HOME/.torch"
REMOTE_TFOLD = "$HOME/tfold"
FREE_STORAGE_GB = 10  # Lightning free Drive allowance; $0.10/GB/mo beyond.

# tFold weight loaders that must be resident before a TCR-pMHC "Complex" prediction.
WEIGHT_LOADERS = "esm_ppi_650m_tcr, tfold_tcr_trunk, tfold_pmhc_trunk, tfold_tcr_pmhc_trunk"


def _connect() -> Studio:
    return Studio(STUDIO_NAME, teamspace=TEAMSPACE, org=ORG)


def _is_cpu(machine) -> bool:
    # Doubles as the "is this machine free?" test. Sound for the only machines the runner
    # ever requests (Machine.CPU is free; L40S is not). If Lightning ever adds a *billed*
    # tier whose name contains "cpu", this heuristic would misclassify it as free and the
    # --gpu guard would not fire on it (PR #1253 review nit).
    return machine == Machine.CPU or "cpu" in str(machine).lower()


def _same_machine(a, b) -> bool:
    """True if a and b are the same Lightning machine (for the CPU/L40S the runner uses)."""
    if _is_cpu(a) or _is_cpu(b):
        return _is_cpu(a) and _is_cpu(b)
    return str(a).lower() == str(b).lower()


def _ensure_running(studio: Studio, machine, allow_gpu: bool) -> bool:
    """Bring the Studio up on `machine`. Returns True if this call started it.

    Safety: refuses to proceed on a paid GPU unless allow_gpu is set, and hard-aborts
    if the resolved machine is not what was asked for.
    """
    if not _is_cpu(machine) and not allow_gpu:
        raise SystemExit(f"refusing to start paid GPU {machine} without --gpu")
    started = False
    if "running" not in str(studio.status).lower():
        print(f"  starting on {machine} ...", flush=True)
        studio.start(machine)
        started = True
    else:
        # Already running. start() sets the machine only on the stopped path, and we
        # deliberately do NOT switch a studio we didn't start: a switch to a paid GPU
        # would leave started=False, so the caller's finally-stop never fires and the GPU
        # bills idle. So a machine mismatch here is a hard abort with a clear message, not
        # a silent run on the wrong machine (e.g. `predict --gpu` against an already-up CPU
        # studio would otherwise run ~79 min on free CPU instead of the L40S). PR #1253 #1.
        current = studio.machine
        if not _is_cpu(current) and not allow_gpu:
            raise SystemExit(f"studio already running on paid GPU {current}; stop it or pass --gpu")
        if not _same_machine(current, machine):
            raise SystemExit(
                f"studio already running on {current}, but {machine} was requested; "
                f"stop it first and re-run so it comes up on the right machine"
            )
    got = studio.machine
    if not _is_cpu(got) and not allow_gpu:
        studio.stop()
        raise SystemExit(f"ABORT: studio came up on {got}, not CPU - stopped to avoid billing")
    return started


def _drive_usage(studio: Studio) -> None:
    """Print persistent Drive usage vs the free cap. $HOME is the billable mount.

    A usage readout must never fail a real run, so the remote command is forced to
    exit 0 (missing dirs are normal before prestage) and any error here is swallowed.
    """
    try:
        out = studio.run(
            f"echo \"USED_MB=$(du -sm $HOME 2>/dev/null | cut -f1)\"; "
            f"echo \"TORCH=$(du -sh {REMOTE_TORCH_HOME} 2>/dev/null | cut -f1)\"; "
            f"echo \"TFOLD=$(du -sh {REMOTE_TFOLD} 2>/dev/null | cut -f1)\"; "
            f"true"
        )
    except Exception as e:  # readout is best-effort, never fatal
        print(f"  Drive: (usage readout unavailable: {e})")
        return
    fields = dict(
        l.split("=", 1) for l in out.splitlines() if "=" in l
    )
    try:
        used_mb = int(fields.get("USED_MB", "").strip())
        pct = 100.0 * used_mb / (FREE_STORAGE_GB * 1024)
        flag = "  <-- OVER FREE CAP" if used_mb > FREE_STORAGE_GB * 1024 else ""
        print(f"  Drive: {used_mb} MB / {FREE_STORAGE_GB} GB free ({pct:.0f}%){flag}")
    except (ValueError, TypeError):
        print(f"  Drive: (could not parse usage: {fields})")
    for label in ("TORCH", "TFOLD"):
        val = fields.get(label, "").strip()
        if val:
            print(f"    {label.lower()} cache: {val}")


def cmd_status(args) -> None:
    studio = _connect()
    print(f"studio {STUDIO_NAME!r} ({ORG}/{TEAMSPACE}): status={studio.status}")
    if "running" in str(studio.status).lower():
        print(f"  machine={studio.machine}")
        _drive_usage(studio)
    else:
        print("  (stopped - start it to read Drive usage)")


def cmd_smoke(args) -> None:
    studio = _connect()
    print(f"[smoke] status={studio.status}")
    started = False
    try:
        started = _ensure_running(studio, Machine.CPU, allow_gpu=False)
        print(f"[smoke] machine={studio.machine}")
        print("[smoke] run():", studio.run("echo SMOKE_OK && whoami").strip())
        with tempfile.TemporaryDirectory() as d:
            up, down = os.path.join(d, "u.txt"), os.path.join(d, "d.txt")
            token = "run-tfold-smoke-roundtrip"
            Path(up).write_text(token + "\n")
            studio.upload_file(up, remote_path="run_tfold_smoke.txt", progress_bar=False)
            studio.download_file("run_tfold_smoke.txt", down)
            ok = os.path.exists(down) and Path(down).read_text().strip() == token
            print(f"[smoke] upload/download round-trip byte-match: {ok}")
        studio.run("rm -f $HOME/run_tfold_smoke.txt")
        _drive_usage(studio)
    finally:
        if started:
            studio.stop()
            print(f"[smoke] stopped; status={studio.status}")


def cmd_prestage(args) -> None:
    studio = _connect()
    print(f"[prestage] status={studio.status}; downloading weights to persistent {REMOTE_TORCH_HOME}")
    started = False
    try:
        started = _ensure_running(studio, Machine.CPU, allow_gpu=False)
        cmd = (
            f"export TORCH_HOME={REMOTE_TORCH_HOME} && cd {REMOTE_TFOLD} && "
            f"python -c \"from tfold.model.pretrain import {WEIGHT_LOADERS}; "
            f"[f() for f in ({WEIGHT_LOADERS},)]\" && "
            f"echo '--- cached weights ---' && ls -lh {REMOTE_TORCH_HOME}/hub/checkpoints/"
        )
        print("[prestage] running weight download on free CPU (this takes ~20 min the first time)...")
        print(studio.run(cmd))
        _drive_usage(studio)
    finally:
        if started:
            studio.stop()
            print(f"[prestage] stopped; status={studio.status}")


def cmd_predict(args) -> None:
    if not os.path.exists(args.input):
        raise SystemExit(f"input not found: {args.input}")
    os.makedirs(args.output_dir, exist_ok=True)
    machine = Machine.L40S if args.gpu else Machine.CPU
    studio = _connect()
    print(f"[predict] input={args.input} machine={machine} status={studio.status}")
    started = False
    try:
        started = _ensure_running(studio, machine, allow_gpu=args.gpu)
        remote_in = "run_tfold_input.json"
        remote_out = "run_tfold_predictions"
        studio.upload_file(args.input, remote_path=remote_in, progress_bar=False)
        cmd = (
            f"export TORCH_HOME={REMOTE_TORCH_HOME} && cd {REMOTE_TFOLD} && "
            f"rm -rf $HOME/{remote_out} && mkdir -p $HOME/{remote_out} && "
            f"python projects/tfold_tcr/predict.py --json $HOME/{remote_in} "
            f"--output $HOME/{remote_out}/ --model_version Complex"
        )
        print("[predict] running predict.py ...")
        print(studio.run(cmd))
        # List the output dir in a SEPARATE call so a chatty predict.py log line ending in
        # ".pdb" can't be misread as an output filename (PR #1253 review #2).
        listing = studio.run(f"ls $HOME/{remote_out}/")
        for name in [l.strip() for l in listing.splitlines() if l.strip().endswith(".pdb")]:
            local = os.path.join(args.output_dir, name)
            studio.download_file(f"{remote_out}/{name}", local)
            print(f"[predict] downloaded {local}")
        _drive_usage(studio)
    finally:
        if started and not args.keep_alive:
            studio.stop()
            print(f"[predict] stopped; status={studio.status}")
        elif args.keep_alive:
            print(f"[predict] --keep-alive: studio left on {studio.machine} (REMEMBER TO STOP IT)")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("status", help="report status + Drive usage").set_defaults(func=cmd_status)
    sub.add_parser("smoke", help="$0 CPU dry-run of the loop").set_defaults(func=cmd_smoke)
    sub.add_parser("prestage", help="download weights to persistent TORCH_HOME on CPU").set_defaults(func=cmd_prestage)
    pp = sub.add_parser("predict", help="run a complex end to end")
    pp.add_argument("--input", required=True, help="path to a tFold input complex JSON")
    pp.add_argument("--output-dir", default="results/tfold", help="local dir for downloaded PDBs")
    pp.add_argument("--gpu", action="store_true", help="use the L40S (billed) instead of free CPU")
    pp.add_argument("--keep-alive", action="store_true", help="do NOT stop the studio on exit")
    pp.set_defaults(func=cmd_predict)
    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
