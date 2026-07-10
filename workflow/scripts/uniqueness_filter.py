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

Why a preflight
---------------
The filter-expression engine (`samtools view -e/--expr`) landed in samtools
1.12. Production provisions samtools via `apt-get` in `scripts/setup_vm.sh`,
which on Ubuntu 22.04 yields 1.13 - new enough, but only by one release. An
older binary must fail loudly rather than silently emit an unfiltered BED12
that looks like a successful filtered run.

These helpers build shell text rather than run anything, so the rendered
`hisat2_align` command is unit-testable without a BAM. Snakemake placeholders
(`{output.bam}`, `{log}`, `{threads}`) are left literal for Snakemake to expand
at job time; nothing here may be passed through `str.format`.
"""

# Filter expressions (`samtools view -e`) were introduced in samtools 1.12.
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
    """
    return (
        "\n"
        "# Preflight: `samtools view -e` (filter expressions) landed in samtools 1.12.\n"
        "# Fail loudly - a silent fallback would emit an unfiltered BED12 that is\n"
        "# indistinguishable from a successful filtered run.\n"
        'if ! samtools view --help 2>&1 | grep -q -- "--expr"; then\n'
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


def build_junction_extraction_block(enabled, expr=NH_UNIQUE_EXPR):
    """Full shell between `samtools index` and the BED12 conversion.

    When `enabled` is False this returns the empty string, so the rendered
    command is byte-identical to the pre-#919 one and no filtered BAM is
    written. When True it returns preflight + prefilter.
    """
    if not enabled:
        return ""
    return build_preflight_block() + build_prefilter_block(expr)


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
