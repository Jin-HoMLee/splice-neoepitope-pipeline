# Open splice-caller benchmark harness (epic #679, leaf B / #966)

The reusable substrate that runs the open splice-neoantigen callers on a shared input, collects their junction + neoepitope-candidate calls into one common schema, and hands a unified table to the downstream concordance leaves (detection #968, burden #969).

This is **leaf B** of the [#679](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/679) decomposition.
Leaf A ([#965](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/965), closed) established which callers run open-only; this leaf turns their outputs into a comparable schema.

## Layout

```
inputs/                   # committed small caller inputs (e.g. the chr22 SJ.out.tab)
harness/
  common_schema.py        # CommonRecord + record_level + junction_id normalization + validate()
  collect.py              # the collect CLI: run and/or ingest -> merge -> validate -> unified TSV
  adapters/
    splice2neo.py         # ingest adapters (caller output -> list[CommonRecord])
    asneo.py
  runners/
    asneo.py              # run adapters (SJ.out.tab -> invoke the caller -> output path)
    run_asneo.sh          # the proven clone/patch/genome/run sequence (subprocess)
  tests/                  # pytest suite for the schema, adapters, runners, and collect CLI
  outputs/                # unified TSVs produced by the CLI (small; demo evidence)
```

## Invocation

`collect.py` takes one or more `caller:path` inputs and writes a unified, plainly tab-delimited TSV:

```bash
cd harness
python collect.py \
  --input splice2neo:../../issue_965_open_caller_audit/outputs/splice2neo_smoke_out.tsv \
  --out outputs/unified_smoke.tsv
```

`--input` is repeatable, so several callers merge into one table:

```bash
python collect.py \
  --input splice2neo:s2n_out.tsv \
  --input asneo:asneo_out.tsv \
  --out outputs/unified.tsv
```

The output columns are the `CommonRecord` fields in declaration order; `provenance` is a JSON string per row.
The format is greppable and `pandas.read_csv(sep="\t")`-friendly (no csv quoting).

### Running a caller end-to-end (AC-1)

`--run CALLER:SJ_TAB` invokes a caller before ingesting, rather than ingesting a pre-made file.
ASNEO (open-only, option-B patched) runs locally on the committed chr22 `SJ.out.tab`:

```bash
cd harness
python collect.py --run asneo:../inputs/chr22_SRR2660032.SJ.out.tab \
  --out outputs/asneo_chr22_unified.tsv
```

The runner clones + patches ASNEO and fetches the hg19 chr22 FASTA into an out-of-tree scratch dir (no >100 MB artifact committed), runs the caller through `conda run -n asneo`, and ingests its peptide-only output.
The committed `outputs/asneo_chr22_unified.tsv` is that run: 6194 chr22 junctions -> 110 pass filters -> 60 novel isoforms -> **800 candidate peptides**.
ASNEO records are `record_level=peptide` with null junction fields (its candidate-peptide set discards the junction linkage; recovery tracked in [#1258](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1258)).

### chr22-scale smoke run (splice2neo ingest)

The committed `outputs/unified_smoke.tsv` is the harness ingesting leaf A's splice2neo toy-fixture smoke output (`issue_965_open_caller_audit/outputs/splice2neo_smoke_out.tsv`): 17 caller rows collapse to 14 candidate records (3 rows carry no in-frame peptide and are dropped by the adapter).
Reproduce it with the first invocation above.

## Adding a caller (the extension point)

1. Write `adapters/<caller>.py` with a function that parses that caller's output into `list[CommonRecord]` (set `record_level` to `"junction"` or `"peptide"` per what the caller emits).
2. Register the adapter with one line in `collect.py`'s `ADAPTERS` dict: `"<caller>": parse_<caller>`.
3. To let the harness *run* the caller (not just ingest it), add `runners/<caller>.py` with `run_<caller>(input, workdir) -> output_path` and register it in the `RUNNERS` dict: `"<caller>": run_<caller>`.
4. Add `tests/test_<caller>_adapter.py` (and a runner test for the pure output-contract helpers) mirroring the ASNEO ones.

No change to `collect.py`'s dispatch/merge/write logic is needed - `get_adapter` raises a clear error for an unregistered caller, and every record is validated against the `record_level` legal combinations at merge.

## Caller status (from leaf A #965)

| Caller | Open-only | Local (arm64 CPU) | Adapter | Runner | Notes |
|---|---|---|---|---|---|
| splice2neo | GO (MIT) | yes | **built** | not yet | junction -> peptide library; `record_level=junction`; hg19 (needs liftOver to hg38) |
| ASNEO | GO (Apache-2.0, #566 patch) | yes (chr22 e2e) | **built** | **built** | `putative_peptide.txt` is peptide-only, so records are `record_level=peptide` with null junction fields; per-peptide junction linkage recovery tracked in [#1258](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/1258) |
| SNAF | GO (MIT) | no (Linux/Docker only) | not yet | not yet | AltAnalyze ships amd64-container-only |
| NeoSplice | NO-GO as-is | no | - | - | Py2.7 + bundled non-redistributable NetMHCpan binaries |

## Test suite

```bash
# from repo root, using the pytest venv
workflow/tests/.venv/bin/python -m pytest research/experiments/issue_679_caller_benchmark/harness
```

`conftest.py` puts the harness dir on `sys.path` so `common_schema` / `adapters` import without packaging.
