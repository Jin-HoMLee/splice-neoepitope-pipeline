"""Pure helpers for the NH-uniqueness BAM prefilter (Issue #919).

`regtools junctions extract` exposes no alignment-quality filter, and by the
time we hold its BED12 output the per-alignment `NH`/`MAPQ` tags are gone. So
suppressing spurious junction calls from multimapped reads has to happen as a
BAM prefilter feeding regtools.

Why `NH`, not a `MAPQ` floor
---------------------------
HISAT2's `MAPQ` conflates two populations (see `unique.h`): `MAPQ=0` is
reachable both from the equal-best-tie branch *and* from the no-second-best
branch, i.e. a uniquely-placed but low-confidence alignment. A `-q 2` floor
therefore discards genuine unique alignments as collateral. The `NH` tag
separates the two exactly, so `[NH]==1` is the primary lever.

Why a preflight, and why it runs *first*
----------------------------------------
The filter-expression engine (`samtools view -e/--expr`) landed in samtools
1.12. Production provisions samtools via `apt-get` in `scripts/setup_vm.sh`,
which on Ubuntu 22.04 yields 1.13 - new enough, but only by one release. An
older binary must fail loudly rather than silently emit an unfiltered BED12
that looks like a successful filtered run.

The gate is emitted at the **top** of the rule's shell, before HISAT2 runs. The
capability probe is pure and costs milliseconds, so there is no reason to burn
an entire alignment before discovering that the binary cannot honor the knob.

The `NH` tag is assumed present
-------------------------------
`[NH]==1` evaluates false for a record with no `NH` tag, so an aligner that omits
it would have every read silently dropped. HISAT2 always emits `NH` (verified on
the chr22 run: all 41,634 mapped records carry it), which is why this is safe
today - but the assumption is load-bearing, and whoever wires the STAR lever must
re-check it rather than inherit it.

These helpers build shell text rather than run anything, so the rendered
`hisat2_align` command is unit-testable without a BAM. Snakemake placeholders
(`{output.bam}`, `{log}`, `{threads}`) are left literal for Snakemake to expand
at job time; nothing here may be passed through `str.format`.
"""

# Filter expressions (`samtools view -e`) were introduced in samtools 1.12.
#
# MIN_SAMTOOLS_VERSION and the two version helpers below are a *spec* of that
# floor, not the live gate: the rendered preflight capability-probes the binary
# (`samtools view --help | grep -- --expr`) instead, so a vendored or backported
# build is judged on what it can actually do rather than on how it numbers itself.
# They are kept, and tested, because the floor is the thing the error message
# cites and the thing a future reader will want to check against upstream.
MIN_SAMTOOLS_VERSION = (1, 12)

# Keep brace-free: this string is embedded in a Snakemake shell template, and a
# `{` would be read as a placeholder.
NH_UNIQUE_EXPR = "[NH]==1"

# Suffix for the prefiltered BAM. Only written when the knob is on.
FILTERED_BAM_SUFFIX = ".nh_unique.bam"


def is_enabled(config):
    """True when `alignment.uniqueness_filter.enabled` is set truthy in config.

    Defaults to False so an untouched config renders exactly today's command.
    """
    alignment = config.get("alignment", {}) or {}
    knob = alignment.get("uniqueness_filter", {}) or {}
    return bool(knob.get("enabled", False))


def filtered_bam_path(bam_path):
    """Derive the prefiltered BAM path from the sorted BAM path.

    `a/b/s.bam` -> `a/b/s.nh_unique.bam`. A path not ending in `.bam` simply
    gets the suffix appended, so a Snakemake placeholder is never mangled.
    """
    if bam_path.endswith(".bam"):
        return bam_path[: -len(".bam")] + FILTERED_BAM_SUFFIX
    return bam_path + FILTERED_BAM_SUFFIX


def parse_samtools_version(version_output):
    """Extract `(major, minor)` from `samtools --version` output.

    Returns None when no version can be parsed. Handles the Debian/Ubuntu
    `1.13-4` packaging suffix and plain `1.23.1`.
    """
    for line in (version_output or "").splitlines():
        line = line.strip()
        if not line.lower().startswith("samtools"):
            continue
        parts = line.split()
        if len(parts) < 2:
            continue
        token = parts[1].split("-")[0]
        bits = token.split(".")
        try:
            return (int(bits[0]), int(bits[1]) if len(bits) > 1 else 0)
        except (ValueError, IndexError):
            return None
    return None


def supports_filter_expressions(version):
    """True when `version` (a `(major, minor)` tuple) has `samtools view -e`."""
    if version is None:
        return False
    return version >= MIN_SAMTOOLS_VERSION


def build_preflight_block():
    """Shell that aborts the rule when samtools lacks `-e`/`--expr`.

    Capability-probes the actual binary (`samtools view --help`) rather than
    parsing a version string, so a vendored or backported build is judged on
    what it can do. `--help` exits non-zero on some builds, hence the pipe into
    grep rather than a bare exit-code test.

    The `|| true` is load-bearing, not defensive noise. This block runs inside a
    shell that opens with `set -euo pipefail`, and under `pipefail` a pipeline
    takes the exit status of its rightmost *failing* component - so a samtools
    whose `--help` exits non-zero would propagate that code even when grep
    matched, `if !` would fire, and the preflight would abort claiming "no
    filter-expression support" about a binary that has it. `pipefail` silently
    re-couples the probe to the very exit code the pipe-into-grep was chosen to
    ignore. Neutralizing samtools' status keeps grep authoritative, which is the
    whole point of the design. (Caught in review on PR #1113; behaviour locked in
    by test_preflight_survives_a_samtools_whose_help_exits_nonzero.)

    The status is neutralized by **capturing into a variable**, not by the more
    obvious `if ! { samtools ...|| true; } | grep`: a brace group would be read by
    Snakemake as a placeholder and fail to expand at job time (the same trap the
    brace-free `NH_UNIQUE_EXPR` avoids). `$(...)` and `"$SAMTOOLS_HELP"` are
    brace-free, and capturing also sidesteps any SIGPIPE edge from `grep -q`
    closing the pipe early.
    """
    return (
        "\n"
        "# Preflight: `samtools view -e` (filter expressions) landed in samtools 1.12.\n"
        "# Fail loudly - a silent fallback would emit an unfiltered BED12 that is\n"
        "# indistinguishable from a successful filtered run.\n"
        "#\n"
        "# Capture first, THEN grep. Under `set -o pipefail` a pipeline takes the exit\n"
        "# status of its rightmost failing component, so piping a `--help` that exits\n"
        "# non-zero straight into grep would falsely abort on a samtools that DOES\n"
        "# support --expr. `|| true` keeps grep authoritative. Captured into a var\n"
        "# rather than a brace group, which Snakemake would read as a placeholder.\n"
        'SAMTOOLS_HELP="$(samtools view --help 2>&1 || true)"\n'
        'if ! printf \'%s\\n\' "$SAMTOOLS_HELP" | grep -q -- "--expr"; then\n'
        '    echo "ERROR: alignment.uniqueness_filter.enabled is true, but this samtools"'
        " | tee -a {log}\n"
        '    echo "       has no filter-expression support (samtools view -e/--expr)."'
        " | tee -a {log}\n"
        '    echo "       Found: $(samtools --version 2>&1 | head -1)" | tee -a {log}\n'
        '    echo "       Need:  samtools >= 1.12 (Ubuntu 22.04 apt ships 1.13)."'
        " | tee -a {log}\n"
        '    echo "       Fix:   upgrade samtools, or set"'
        " | tee -a {log}\n"
        '    echo "              alignment.uniqueness_filter.enabled: false" | tee -a {log}\n'
        "    exit 1\n"
        "fi\n"
    )


def build_prefilter_block(expr=NH_UNIQUE_EXPR):
    """Shell that writes and indexes the NH-unique BAM regtools will read.

    regtools requires an indexed BAM, so the filtered alignments must be
    materialized rather than streamed.
    """
    return (
        "\n"
        "# Uniqueness prefilter: drop multimapped alignments before junction\n"
        "# extraction. regtools takes an indexed BAM, so this is materialized.\n"
        "samtools view \\\n"
        "    -b \\\n"
        "    -@ {threads} \\\n"
        "    -e '" + expr + "' \\\n"
        "    -o {output.filtered_bam} \\\n"
        "    {output.bam} \\\n"
        "    2>> {log}\n"
        "\n"
        "samtools index {output.filtered_bam} 2>> {log}\n"
        "\n"
        'echo "Uniqueness filter (' + expr + '): '
        "$(samtools view -c {output.bam}) -> "
        '$(samtools view -c {output.filtered_bam}) alignments" >> {log}\n'
    )


def build_preflight_gate(enabled):
    """Shell emitted at the TOP of the rule, before HISAT2 runs.

    Empty when the knob is off. When on, aborts before any expensive work if the
    binary cannot honor the filter - failing after align+sort+index would waste
    the whole alignment to learn something a millisecond probe can tell us.
    """
    return build_preflight_block() if enabled else ""


def build_prefilter_gate(enabled, expr=NH_UNIQUE_EXPR):
    """Shell emitted between `samtools index` and the BED12 conversion.

    Empty when the knob is off, so the rendered command is byte-identical to the
    pre-#919 one and no filtered BAM is written.
    """
    return build_prefilter_block(expr) if enabled else ""


def regtools_input_bam(enabled):
    """The Snakemake placeholder regtools should read alignments from."""
    return "{output.filtered_bam}" if enabled else "{output.bam}"


def filtered_bam_outputs(enabled, bam_template):
    """Extra `output:` entries for the align rule, keyed for `**` unpacking.

    Empty when the knob is off, so Snakemake declares (and later cleans) no
    filtered BAM at all. `bam_template` is the rule's `bam` output path, which
    may contain `{patient_id}`/`{sample}` wildcards.
    """
    if not enabled:
        return {}
    filtered = filtered_bam_path(bam_template)
    return {"filtered_bam": filtered, "filtered_bai": filtered + ".bai"}
