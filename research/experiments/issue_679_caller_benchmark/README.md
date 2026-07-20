# Open splice-caller benchmark harness (epic #679, leaf B / #966)

The reusable substrate that runs the open splice-neoantigen callers on a shared input, collects their junction + neoepitope-candidate calls into one common schema, and hands a unified table to the downstream concordance leaves (detection #968, burden #969).

This is **leaf B** of the [#679](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/679) decomposition.
Leaf A ([#965](https://github.com/Jin-HoMLee/splice-neoepitope-pipeline/issues/965), closed) established which callers run open-only; this leaf turns their outputs into a comparable schema.

## Layout

```
harness/
  common_schema.py        # CommonRecord + junction_id normalization (the shared contract)
  collect.py              # the collect CLI: dispatch adapters -> merge -> unified TSV
  adapters/
    splice2neo.py         # one adapter per caller (caller output -> list[CommonRecord])
  tests/                  # pytest suite for the schema, adapters, and collect CLI
  outputs/                # unified TSVs produced by the CLI (small; demo evidence)
```

## Invocation

`collect.py` takes one or more `caller:path` inputs and writes a unified, plainly tab-delimited TSV:

```bash
python collect.py \
  --input splice2neo:../issue_965_open_caller_audit/outputs/splice2neo_smoke_out.tsv \
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

### chr22-scale smoke run

The committed `outputs/unified_smoke.tsv` is the harness run on leaf A's real chr22-scale splice2neo smoke output (`issue_965_open_caller_audit/outputs/splice2neo_smoke_out.tsv`): 17 caller rows collapse to 14 candidate records (3 rows carry no in-frame peptide and are dropped by the adapter).
Reproduce it with the first invocation above.

## Adding a caller (the extension point)

1. Write `adapters/<caller>.py` with a function that parses that caller's output into `list[CommonRecord]`.
2. Register it with one line in `collect.py`'s `ADAPTERS` dict: `"<caller>": parse_<caller>`.
3. Add a `tests/test_<caller>_adapter.py` mirroring `tests/test_splice2neo_adapter.py`.

No change to `collect.py`'s dispatch/merge/write logic is needed - `get_adapter` raises a clear error for an unregistered caller.

## Caller status (from leaf A #965)

| Caller | Open-only | Local (arm64 CPU) | Adapter | Notes |
|---|---|---|---|---|
| splice2neo | GO (MIT) | yes | **built** | junction -> peptide library; hg19 (needs liftOver to hg38) |
| ASNEO | GO (Apache-2.0, #566 patch) | yes (chr22 smoke) | not yet | `putative_peptide.txt` is peptide-only - **no per-peptide junction coords**, so it does not fill the required `junction_id`. A schema-identity decision precedes this adapter (tracked separately). |
| SNAF | GO (MIT) | no (Linux/Docker only) | not yet | AltAnalyze ships amd64-container-only |
| NeoSplice | NO-GO as-is | no | - | Py2.7 + bundled non-redistributable NetMHCpan binaries |

## Test suite

```bash
# from repo root, using the pytest venv
workflow/tests/.venv/bin/python -m pytest research/experiments/issue_679_caller_benchmark/harness
```

`conftest.py` puts the harness dir on `sys.path` so `common_schema` / `adapters` import without packaging.
