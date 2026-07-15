# Issue #1140 - hook-liveness contract specs (workflow provenance)

Provenance record for the per-hook contract specs that back `tools/ci/test_hook_liveness_contract.py`.
Related: [Issue #1140](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1140), epic [#1135](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1135), [PR #1195](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/pull/1195), design [`docs/superpowers/specs/2026-07-15-hook-liveness-contract-tests-design.md`](../../../docs/superpowers/specs/2026-07-15-hook-liveness-contract-tests-design.md).

## How these were produced

A 22-agent workflow (one **analyze** + one adversarial **verify** agent per hook, over the 11 hooks registered in `.agents/settings.json`). Each analyze agent read its hook's real internal matcher and *drove* the hook through `hook_contract.drive()` to produce an empirically-validated `HookContract` spec (real fire trigger, real close counterexample, observable kind, `gh`-stub need). Each verify agent independently re-drove all three directions and ran a per-hook red-on-break check. The main session assembled the verified specs into the single `CONTRACTS` registry in the suite.

The assembled, in-use form of this data lives in the committed suite; these files are the raw record, kept because the workflow scratchpad is ephemeral.

## Files

- `workflow_agent_results.json` - the full per-hook `{spec, verify}` records (analyze spec + adversarial-verify verdict, incl. each agent's evidence and rationale). 11 hooks.
- `verified_specs.json` - the merged effective spec per hook (verify corrections applied over the analyze spec), the intermediate that fed registry assembly.

## What the adversarial pass caught (the value of the verify stage)

Two real harness gaps that the analyze specs alone would have shipped:

1. **`check_memory_path_cwd_drift`** keys on the session `cwd` in the PreToolUse envelope, which the harness's `pretooluse_bash` did not carry - so the contract could never fire through the harness as first built. Fixed by adding an optional `cwd` to `pretooluse_bash`.
2. **`write_session_watermark`** has no discriminating matcher (it writes unconditionally on Stop), and the first observable helper always forced a writable `CLAUDE_PROJECT_DIR`, so the hook could never be made silent - the counterexample was structurally impossible. Fixed by making the matched pair writable-root (writes) vs unwritable-root (fails open, silent).

One hook (`write_session_watermark`) came back `confirmed=false` from verify for exactly reason 2 above, and one (`check_memory_path_cwd_drift`) carried corrections; both were reconciled before the suite's first run. All 11 red-on-break checks passed.
