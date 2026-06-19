#!/usr/bin/env python3
"""Deep project structure extractor — fully granular Pipeline Core layer.

Extracts per-rule, per-function, per-output-schema, per-dependency-edge detail
from the entire project tree, with special depth inside workflow/.

Usage:
    python tools/project_map/extract_graph.py [--output tools/project_map/graph.json]
"""

import json
import os
import re
import subprocess
import sys
from pathlib import Path

# This file lives at tools/project_map/extract_graph.py → repo root is 3 levels up.
PROJECT_ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent

# Our own build artifacts / vendored deps are not "project structure" — skip them
# so the map doesn't list a 280 KB minified D3 blob or its own generated outputs.
SELF_EXCLUDE_FILES = {"tools/project_map/graph.json", "tools/project_map/project_map.html"}
SELF_EXCLUDE_DIRS = {"tools/project_map/vendor"}


def slug(name):
    return name.replace("/", "__").replace(".", "_").replace("-", "_")


def git_tracked_paths():
    """Return ``(tracked_files, tracked_dirs)`` from ``git ls-files``.

    Discovery is restricted to the git-tracked set so gitignored runtime/data
    dirs (``references/``, ``results/``, ``logs/``, ``data/``, ``indices/``)
    never enter the atlas — making it reproducible across working-tree states
    rather than dependent on whether the pipeline has been run locally (#780).
    git naturally respects ``.gitignore``, so this is a single source of truth.

    Both sets hold repo-root-relative, forward-slash paths. ``tracked_dirs`` is
    every ancestor directory of a tracked file (so the walk can keep a dir that
    only contains tracked files deeper down).
    """
    try:
        out = subprocess.run(
            ["git", "-C", str(PROJECT_ROOT), "ls-files", "-z"],
            capture_output=True, text=True, check=True).stdout
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        stderr_hint = (getattr(exc, "stderr", "") or "").strip()
        raise RuntimeError(
            "project_map's extractor requires a git working tree "
            "(`git ls-files` failed); run it from inside the repo."
            + (f"\ngit stderr: {stderr_hint}" if stderr_hint else "")) from exc
    tracked_files = {p for p in out.split("\0") if p}
    tracked_dirs = set()
    for path in tracked_files:
        parts = path.split("/")
        for i in range(1, len(parts)):
            tracked_dirs.add("/".join(parts[:i]))
    return tracked_files, tracked_dirs


ISSUE_RE = re.compile(r'(?:^|/)issue_(\d+)')


def issue_of(rel_path):
    """The GitHub issue number an issue_NNN/ path belongs to, or None."""
    m = ISSUE_RE.search(rel_path)
    return int(m.group(1)) if m else None


def classify_path(rel_path, is_dir=False):
    p = rel_path.lower()
    if is_dir:
        if p == "workflow/rules": return "rules_dir"
        if p == "workflow/scripts": return "scripts_dir"
        if p == "workflow/envs": return "envs_dir"
        if p.startswith("workflow/tests"): return "tests_dir"
        if p == "config/samples": return "samples_dir"
        if p.startswith("research/evals"): return "evals_dir"
        if p.startswith("research/experiments"): return "experiments_dir"
        if p.startswith("research/decisions"): return "decisions_dir"
        if p.startswith("research/manuscript"): return "manuscript_dir"
        if p.startswith("research/notebooks"): return "notebooks_dir"
        if p.startswith("research/lab_notebook"): return "labnotebook_dir"
        if p.startswith("research/slides"): return "slides_dir"
        if p.startswith("tools/ci"): return "ci_dir"
        if p.startswith("tools/news"): return "news_dir"
        if p == "docker": return "docker_dir"
        if p == "resources": return "resources_dir"
        if p.startswith("scripts/pm"): return "pm_dir"
        if p.startswith("docs"): return "docs_dir"
        if p == "config": return "config_dir"
        return "dir"

    # ── file classification ──
    # os.walk passes relative paths ("research/x.md"); the directory tests below
    # are written for absolute-style ("/research/"), so normalise with a leading
    # slash. basename checks use `base`.
    q = "/" + p
    base = p.rsplit("/", 1)[-1]

    if q.endswith(".smk"): return "rule"
    if q.endswith((".py", ".sh")) and "/workflow/scripts/" in q: return "script"
    if q.endswith(".py") and "/scripts/" in q and "/pm/" not in q: return "script"
    if q.endswith(".yaml") and "/workflow/envs/" in q: return "env"
    if q.endswith(".yaml") and "/config/" in q: return "config"
    if q.endswith(".tsv") and "/config/samples/" in q: return "sample"
    if q.endswith(".tsv"): return "data"
    if q.endswith(".md") and "/docs/" in q: return "doc"
    if q.endswith(".md") and "/research/manuscript/" in q: return "manuscript"
    if q.endswith(".md") and "/research/lab_notebook/" in q: return "labnotebook"
    if q.endswith("/lab_notebook.md"): return "labnotebook"
    if q.endswith(".md") and "/research/" in q: return "doc"  # research-root notes
    if q.endswith(".ipynb") and "/research/" in q: return "notebook"
    if q.endswith(".qmd") and ("/research/evals/" in q or "/research/experiments/" in q
                               or "/research/decisions/" in q or "/research/slides/" in q
                               or "/docs/features/" in q): return "slide"
    if q.endswith(".sh") and "/scripts/" in q: return "script"
    if q.endswith(".py") and "/tools/ci/" in q: return "ci"
    if q.endswith(".py") and "/tools/news/" in q: return "news"
    if q.endswith(".py") and "/workflow/tests/" in q: return "test"
    if q.endswith(".py") and "/scripts/tests/" in q: return "test"
    if q.endswith(".ini"): return "config"
    if base == "snakefile": return "rule"
    if base == "dockerfile.pipeline": return "docker"
    if q.endswith(".json") and "/scripts/" in q: return "config"
    if q.endswith((".csl", ".scss")): return "resource"
    if q.endswith(".bib"): return "reference"
    if base == "readme.md": return "doc"
    if q.endswith(".yaml") and "/tools/news/" in q: return "config"
    return "resource"


def assign_group(rel_path):
    if rel_path.startswith("workflow"): return "pipeline"
    if rel_path.startswith("config"): return "config"
    if rel_path.startswith("research"): return "research"
    if rel_path.startswith("docs"): return "docs"
    if rel_path.startswith("scripts") or rel_path.startswith("tools"): return "infrastructure"
    if rel_path.startswith("docker") or rel_path.startswith("resources"): return "infrastructure"
    return "project"


# Top-level dirs nest under their group's root node, building a complete
# containment tree (group root → top dir → subdir → file).
GROUP_ROOT_ID = {
    "pipeline": "pipeline_root",
    "research": "research_root",
    "infrastructure": "infrastructure_root",
    "docs": "docs_root",
    "config": "config_root",
}


def describe_file(rel_path, ftype):
    labels = {
        "rule": "Snakemake rule module",
        "script": "Python/shell script",
        "env": "Conda environment definition",
        "config": "Configuration file",
        "sample": "Sample manifest",
        "doc": "Documentation",
        "manuscript": "Manuscript section",
        "notebook": "Jupyter analysis notebook",
        "labnotebook": "Lab notebook entry",
        "slide": "Quarto slide deck",
        "test": "Test module",
        "ci": "CI automation tool",
        "news": "Release polling tool",
        "docker": "Docker build configuration",
        "reference": "Bibliography / reference",
        "data": "Tabular data",
    }
    return labels.get(ftype, "Project file")


# ── Deep pipeline analysis ────────────────────────────────────────────────

def get_script_metadata(script_path):
    """Extract functions, imports, argparse args, and complexity info from a pipeline script."""
    path = PROJECT_ROOT / script_path
    if not path.exists():
        return {}, [], [], 0
    content = path.read_text()
    functions = list(re.findall(r'^def\s+(\w+)\(', content, re.MULTILINE))
    # sorted(), not list(set(...)): bare set ordering is PYTHONHASHSEED-dependent,
    # so an unsorted serialization makes graph.json regenerate non-identically
    # across runs/machines — breaking the committed-vs-regenerated invariant.
    imports = sorted(set(re.findall(r'^import\s+(\S+)', content, re.MULTILINE) +
                         re.findall(r'^from\s+(\S+)\s+import', content, re.MULTILINE)))
    argparse_args = [m.group(1) for m in re.finditer(r"add_argument\(['\"]-{1,2}(\w+)['\"]", content)]
    lines = content.count("\n") + 1
    return functions, imports, argparse_args, lines


def get_pipeline_env_deps(env_path):
    """Return package deps from a conda yaml."""
    path = PROJECT_ROOT / env_path
    if not path.exists():
        return []
    content = path.read_text()
    deps = re.findall(r'^\s*-\s+(\S+)', content, re.MULTILINE)
    return [d for d in deps if not d.startswith('#') and not d.startswith('-')]


def parse_rule_resources(content):
    """Extract threads, mem_mb, and conda env from a rule block.

    Returns a dict with keys: threads (int|None), mem_mb (int|None),
    gpu (bool), docker (str|None), conditional (bool).
    """
    resources = {"threads": None, "mem_mb": None, "gpu": False,
                 "docker": None, "conditional": False, "config_keys": []}

    # Match threads: <number> or threads: config.get(...)
    tm = re.search(r'threads:\s*(\d+)', content)
    if tm:
        resources["threads"] = int(tm.group(1))
    else:
        # threads: config.get("alignment", {}).get("threads", 8)
        # The matcher regex has no capture group, so the default must be pulled
        # with a separate capturing search — never call tm.group(1) on it.
        tm_cfg = re.search(
            r'threads:\s*.*\.get\(["\']\w+["\'],\s*["\']?\w+["\']?\s*\)', content)
        if tm_cfg:
            tm2 = re.search(r'\.get\(["\']\w+["\'],\s*(\d+)\)', content)
            if tm2:
                resources["threads"] = int(tm2.group(1))

    # Match resources: block with mem_mb=N
    mm = re.search(r'resources:\s*\n\s+mem_mb=(\d+)', content)
    if mm:
        resources["mem_mb"] = int(mm.group(1))

    # GPU detection: check for enabled gating via gpu_config or Docker GPU image
    if re.search(r'config\[["\']\w+["\']\]\[["\']enabled["\']\]', content):
        resources["conditional"] = True
    if "gpu" in content.lower() or "nvidia" in content.lower() or "cuda" in content.lower():
        resources["gpu"] = True
    if "tcrdock" in content.lower():
        resources["gpu"] = True

    dm = re.search(r'docker_image\s*=\s*config\[["\'](\w+)["\']\]\[["\'](\w+)["\']\]', content)
    if dm:
        resources["docker"] = dm.group(0)

    # Config keys accessed (snakemake config["..."]["..."])
    for cm in re.finditer(r'config\[["\'](\w+)["\']\]\[["\'](\w+)["\']\]', content):
        key = f"{cm.group(1)}.{cm.group(2)}"
        if key not in resources["config_keys"]:
            resources["config_keys"].append(key)
    for cm in re.finditer(r'config\[["\'](\w+)["\']\]', content):
        key = cm.group(1)
        if key not in resources["config_keys"]:
            resources["config_keys"].append(key)

    return resources


MD_LINK_RE = re.compile(r'\]\(\s*([^)\s]+\.md)(?:#[^)]*)?\s*\)')


def resolve_link(src_rel, link):
    """Resolve a markdown link target to a repo-relative path, or None if external/anchor."""
    link = link.strip()
    if not link or link.startswith(("http://", "https://", "#", "mailto:")):
        return None
    link = link.split("#", 1)[0]
    if not link:
        return None
    joined = os.path.normpath(os.path.join(os.path.dirname(src_rel), link))
    return joined.replace(os.sep, "/")


def parse_markdown_links(node_ids):
    """Scan docs/ + research/ markdown for links to other repo .md files → references edges."""
    edges = []
    for base in ("docs", "research"):
        base_dir = PROJECT_ROOT / base
        if not base_dir.exists():
            continue
        for md in base_dir.rglob("*.md"):
            src_rel = str(md.relative_to(PROJECT_ROOT))
            src_id = slug(src_rel)
            if src_id not in node_ids:
                continue
            try:
                text = md.read_text()
            except Exception:
                continue
            for m in MD_LINK_RE.finditer(text):
                tgt_rel = resolve_link(src_rel, m.group(1))
                if not tgt_rel:
                    continue
                tgt_id = slug(tgt_rel)
                if tgt_id in node_ids and tgt_id != src_id:
                    edges.append({"source": src_id, "target": tgt_id,
                                  "type": "references", "label": "links to"})
    return edges


PIPELINE_SCRIPTS_META = {
    "workflow/scripts/filter_junctions.py": {
        "desc": "Classify splice junctions by origin (annotated/normal_shared/tumor_exclusive)",
        "inputs": "raw_junctions.tsv × N, manifest.tsv, reference_junctions.bed, gencode.gtf",
        "outputs": "novel_junctions.tsv (junction_id, chrom, start, end, strand, mapped_reads, sample_id, sample_type, junction_origin, reading_frame)",
    },
    "workflow/scripts/build_reference_junctions.py": {
        "desc": "Parse GENCODE GTF → reference junction BED (chrom, start, end, name, score, strand)",
        "inputs": "gencode.v47.annotation.gtf.gz",
        "outputs": "reference_junctions.bed",
    },
    "workflow/scripts/assemble_contigs.py": {
        "desc": "Extract junction-spanning contigs via bedtools getfasta, exclude soft-clipped",
        "inputs": "novel_junctions.tsv, GRCh38.fa",
        "outputs": "contigs.fa (FASTA, 54 nt contigs = 27+27 symmetric flanks)",
    },
    "workflow/scripts/translate_peptides.py": {
        "desc": "Translate junction contigs → spanning 8/9/10-mers in 3 reading frames, stop-codon filter",
        "inputs": "contigs.fa",
        "outputs": "peptides.tsv (contig_key, start_nt, peptide)",
    },
    "workflow/scripts/proteome_filter.py": {
        "desc": "Remove self-peptides via k-mer index against UniProt Swiss-Prot",
        "inputs": "peptides.tsv, human_proteome.fasta",
        "outputs": "peptides_novel.tsv, peptides_excluded.tsv (with UniProt accessions)",
    },
    "workflow/scripts/run_mhcflurry.py": {
        "desc": "MHCflurry 2.x Class1PresentationPredictor — per-allele + genotype presentation scores",
        "inputs": "peptides_novel.tsv or peptides.tsv, alleles.tsv or fallback_alleles",
        "outputs": "mhc_presentation.tsv (contig_key, peptide, best_allele, ic50_nM, processing_score, presentation_score, presentation_percentile, presentation_class, genotype_presentation_score)",
    },
    "workflow/scripts/aggregate_hla_alleles.py": {
        "desc": "Aggregate per-sample OptiType → per-patient alleles.tsv with normal-first policy + QC",
        "inputs": "{sample}_result.tsv × N, samples.tsv",
        "outputs": "alleles.tsv (locus, allele1, allele2), hla_qc.tsv (source, reads, discrepancy)",
    },
    "workflow/scripts/run_tcrdock.py": {
        "desc": "TCRdock AlphaFold-based TCR-pMHC ternary structure prediction (GPU only)",
        "inputs": "mhc_presentation.tsv, vdjdb_panel.tsv (optional)",
        "outputs": "top_candidate.pdb, docking_scores.tsv (pLDDT, PAE, pTM, TCR provenance)",
    },
    "workflow/scripts/fetch_vdjdb_panel.py": {
        "desc": "Filter VDJdb → per-patient HLA-matched paired α/β TCR panel via stitchr",
        "inputs": "vdjdb_full.txt, IMGT germlines, alleles.tsv",
        "outputs": "panel.tsv (TCRs with full sequences), panel_qc.tsv (per-allele coverage)",
    },
    "workflow/scripts/generate_report.py": {
        "desc": "Self-contained HTML report: junction summary, HLA QC, top presenters, Mol* 3D viewer",
        "inputs": "novel_junctions.tsv, filtering_stats.tsv, mhc_presentation.tsv, hla_qc.tsv, vdjdb_panel.tsv, contigs.fa, tcrdock outputs",
        "outputs": "report.html, report.tsv, report_top_candidates.tsv",
    },
    "workflow/scripts/aggregate_filtering_stats.py": {
        "desc": "Combine per-step funnel stats into unified filtering_stats.tsv",
        "inputs": "junction_filter_stats.tsv, contig_assemble_stats.tsv, translate_stats.tsv, mhc_stats.tsv, proteome_stats.tsv",
        "outputs": "filtering_stats.tsv (patient_id, sample_id, sample_type, step, category, count)",
    },
    "workflow/scripts/strandness.py": {
        "desc": "Map biological strandness → HISAT2 --rna-strandness flag (F/R/FR/RF)",
        "inputs": "samples.tsv strandness column",
        "outputs": "strandness flag string",
    },
    "workflow/scripts/bed12_to_junctions.py": {
        "desc": "Convert regtools BED12 → 2-column junctions TSV (HISAT2 path)",
        "inputs": "regtools BED12",
        "outputs": "raw_junctions.tsv",
    },
    "workflow/scripts/star_sj_to_junctions.py": {
        "desc": "Convert STAR SJ.out.tab → 2-column junctions TSV with strand rescue from intron motif",
        "inputs": "SJ.out.tab",
        "outputs": "raw_junctions.tsv",
    },
    "workflow/scripts/build_gtex_pan_tissue_ref.py": {
        "desc": "Build GTEx pan-tissue novel-junction blacklist BED from Snaptron gtexv2",
        "inputs": "Snaptron bulk bgz (remote tabix)",
        "outputs": "pan-tissue junction BED",
    },
    "workflow/scripts/run_star_alignment.py": {
        "desc": "STAR genome index build + alignment on FASTQ samples",
        "inputs": "FASTQ files, STAR index",
        "outputs": "SJ.out.tab → raw_junctions.tsv (via star_sj_to_junctions.py)",
    },
}

PIPELINE_RULES_DETAIL = {
    "build_reference_junctions": {
        "smk": "filter_junctions.smk",
        "desc": "Parse GENCODE GTF → reference junction BED",
        "type": "reference build",
        "outputs": ["references/reference_junctions.bed"],
    },
    "filter_junctions": {
        "smk": "filter_junctions.smk",
        "desc": "Classify junctions by origin: annotated→discard, normal_shared, tumor_exclusive",
        "type": "core filter",
        "outputs": ["junctions/novel_junctions.tsv", "junctions/junction_filter_stats.tsv"],
    },
    "create_alignment_manifest": {
        "smk": "alignment.smk",
        "desc": "Create manifest TSV from samples config",
        "type": "metadata",
        "outputs": ["alignment/manifest.tsv"],
    },
    "hisat2_download_index": {
        "smk": "alignment.smk",
        "desc": "Download pre-built HISAT2 hg38 index",
        "type": "index",
        "outputs": ["indices/hisat2/"],
    },
    "hisat2_index": {
        "smk": "alignment.smk",
        "desc": "Build HISAT2 genome index (~8 GB RAM)",
        "type": "index",
        "outputs": ["indices/hisat2/"],
    },
    "hisat2_align": {
        "smk": "alignment.smk",
        "desc": "HISAT2 alignment + regtools junction extraction",
        "type": "alignment",
        "outputs": ["alignment/{sample}/raw_junctions.tsv", "alignment/{sample}/{sample}.bam", "alignment/{sample}/{sample}.bam.bai", "alignment/{sample}/{sample}_junctions.bed"],
    },
    "star_index": {
        "smk": "alignment.smk",
        "desc": "Build STAR genome index (~32 GB RAM)",
        "type": "index",
        "outputs": ["indices/star/"],
    },
    "star_align": {
        "smk": "alignment.smk",
        "desc": "STAR alignment with 2-pass mode",
        "type": "alignment",
        "outputs": ["alignment/{sample}/raw_junctions.tsv"],
    },
    "run_optitype": {
        "smk": "hla_typing.smk",
        "desc": "OptiType HLA-A/B/C typing on sample FASTQ",
        "type": "hla",
        "outputs": ["hla_typing/{sample}/{sample}_result.tsv", "hla_typing/{sample}/{sample}_coverage_plot.pdf"],
    },
    "aggregate_hla_alleles": {
        "smk": "hla_typing.smk",
        "desc": "Aggregate OptiType → per-patient alleles with normal-first policy",
        "type": "hla",
        "outputs": ["hla_typing/alleles.tsv", "hla_typing/hla_qc.tsv"],
    },
    "assemble_contigs": {
        "smk": "assemble_contigs.smk",
        "desc": "Extract junction-spanning 54 nt contigs via bedtools getfasta",
        "type": "core assembly",
        "outputs": ["contigs/contigs.fa", "contigs/contig_assemble_stats.tsv"],
    },
    "translate_peptides": {
        "smk": "translate_peptides.smk",
        "desc": "Translate contigs → 8/9/10-mers in 3 reading frames",
        "type": "core translation",
        "outputs": ["peptides/peptides.tsv", "peptides/translate_stats.tsv"],
    },
    "download_human_proteome": {
        "smk": "proteome_filter.smk",
        "desc": "Download UniProt Swiss-Prot human proteome FASTA",
        "type": "download",
        "outputs": ["references/human_proteome.fasta"],
    },
    "proteome_filter_peptides": {
        "smk": "proteome_filter.smk",
        "desc": "Remove self-peptides via k-mer index against Swiss-Prot",
        "type": "core filter",
        "outputs": ["peptides/peptides_novel.tsv", "peptides/peptides_excluded.tsv", "peptides/proteome_stats.tsv"],
    },
    "download_mhcflurry_models": {
        "smk": "mhc_affinity.smk",
        "desc": "Download MHCflurry trained models (~1 GB)",
        "type": "download",
        "outputs": ["~/.mhcflurry/.download_done"],
    },
    "run_mhcflurry": {
        "smk": "mhc_affinity.smk",
        "desc": "MHCflurry 2.x prediction — per-allele + genotype presentation scores",
        "type": "core prediction",
        "outputs": ["predictions/mhc_presentation.tsv", "predictions/mhc_stats.tsv"],
    },
    "fetch_vdjdb_panel": {
        "smk": "tcr_panel.smk",
        "desc": "Build per-patient VDJdb TCR panel (HLA-matched paired α/β)",
        "type": "tcr reference",
        "outputs": ["tcr_panel/vdjdb/panel.tsv", "tcr_panel/vdjdb/panel_qc.tsv"],
    },
    "download_vdjdb_release": {
        "smk": "download.smk",
        "desc": "Download pinned VDJdb release + SHA256 verify",
        "type": "download",
        "outputs": ["references/vdjdb/{release}/vdjdb_full.txt"],
    },
    "download_imgt_germlines": {
        "smk": "download.smk",
        "desc": "Download IMGT germline reference via stitchrdl",
        "type": "download",
        "outputs": ["references/imgt_germlines/.download.done"],
    },
    "download_fastq": {
        "smk": "download.smk",
        "desc": "Download FASTQ from GCS or public HTTPS",
        "type": "download",
        "outputs": ["data/{patient_id}/{sample}/{filename}"],
    },
    "run_tcrdock": {
        "smk": "structure.smk",
        "desc": "TCRdock AlphaFold TCR-pMHC ternary structure (GPU/Docker)",
        "type": "structural",
        "outputs": ["predictions/tcrdock/top_candidate.pdb", "predictions/tcrdock/docking_scores.tsv"],
    },
    "aggregate_filtering_stats": {
        "smk": "analysis.smk",
        "desc": "Concatenate per-step funnel stats into unified filtering_stats.tsv",
        "type": "analysis",
        "outputs": ["reports/filtering_stats.tsv"],
    },
    "generate_report": {
        "smk": "analysis.smk",
        "desc": "Self-contained HTML report with Mol* 3D viewer",
        "type": "analysis",
        "outputs": ["reports/report.html", "reports/report.tsv", "reports/report_top_candidates.tsv"],
    },
}

PIPELINE_ENVS_DETAIL = {
    "workflow/envs/python.yaml": {
        "packages": ["python >=3.11,<3.13", "biopython >=1.81", "pandas >=2.0", "numpy >=1.25",
                     "scipy >=1.11", "matplotlib >=3.8", "seaborn >=0.13", "requests >=2.31",
                     "jinja2 >=3.1", "pytest >=8.0", "torch ==2.12.0+cu126", "mhcflurry >=2.0"],
        "channel": "conda-forge + bioconda + pip (pytorch cu126)",
    },
    "workflow/envs/biotools.yaml": {
        "packages": ["bedtools >=2.31", "pandas >=2.0", "python >=3.11,<3.13"],
        "channel": "bioconda + conda-forge",
    },
    "workflow/envs/hisat2.yaml": {
        "packages": ["hisat2 >=2.2.1", "python >=3.11,<3.13", "pandas", "regtools >=1.0.0"],
        "channel": "bioconda + conda-forge",
        "notes": "samtools intentionally omitted (libdeflate ABI conflict, tracked in #237)",
    },
    "workflow/envs/star.yaml": {
        "packages": ["star=2.7.10b", "python >=3.11", "pandas"],
        "channel": "bioconda + conda-forge",
        "notes": "Pinned to 2.7.10b (not 2.7.11+) due to libdeflate ABI conflict (#629)",
    },
    "workflow/envs/optitype.yaml": {
        "packages": ["optitype =1.3.5", "coincbc"],
        "channel": "bioconda + conda-forge",
    },
    "workflow/envs/vdjdb.yaml": {
        "packages": ["python=3.11", "pandas", "requests", "stitchr>=1.1.0", "IMGTgeneDL>=0.6.1"],
        "channel": "conda-forge + bioconda + pip",
    },
    "workflow/envs/alphagenome.yaml": {
        "packages": ["python >=3.11,<3.13", "pandas >=2.0", "numpy >=1.25", "matplotlib >=3.8",
                     "matplotlib-venn >=1.1", "scipy >=1.11", "seaborn >=0.13", "scikit-learn >=1.3",
                     "pyarrow >=15.0", "ipykernel", "alphagenome >=0.6,<0.7"],
        "channel": "conda-forge + pip",
    },
}


def parse_smk_rules():
    """Parse Snakemake .smk files for rule definitions and dependencies."""
    edges = []
    rules_dir = PROJECT_ROOT / "workflow" / "rules"
    if not rules_dir.exists():
        return edges

    for smk_file in sorted(rules_dir.glob("*.smk")):
        rel_path = str(smk_file.relative_to(PROJECT_ROOT))
        content = smk_file.read_text()
        smk_file_id = slug(rel_path)

        # Find rule names inside this file (including indented rules inside if/elif blocks)
        _entries = []
        for m in re.finditer(r'(?:^|\n)\s*rule\s+(\w+):', content):
            rule_name = m.group(1)
            rule_id = f"rule__{rule_name}"

            # Create a sub-node for the individual rule
            edges.append({
                "source": rule_id, "target": smk_file_id,
                "type": "defined_in", "label": "defined in rule file",
            })

            # Collect position info for block extraction
            _entries.append((rule_name, rule_id, m.end()))

        # Extract script/conda/dep edges using position-based block extraction
        for i, (rule_name, rule_id, block_start) in enumerate(_entries):
            next_pos = _entries[i + 1][2] + 1 if i + 1 < len(_entries) else len(content)
            rule_block = content[block_start:next_pos]
            for sm in re.finditer(r'script:\s*"\.\./scripts/(.+?)"', rule_block):
                script_name = sm.group(1)
                script_id = slug(f"workflow/scripts/{script_name}")
                edges.append({
                    "source": rule_id, "target": script_id,
                    "type": "calls", "label": "runs script",
                })

            # Extract conda: directives
            for cm in re.finditer(r'conda:\s*"\.\./envs/(.+?)"', rule_block):
                env_name = cm.group(1)
                env_id = slug(f"workflow/envs/{env_name}")
                edges.append({
                    "source": rule_id, "target": env_id,
                    "type": "runs_in", "label": "uses env",
                })

            # Extract input: references to other rule outputs
            for im in re.finditer(r'input:\s*(.+?)(?=\n\s*(?:output|log|params|conda|script|run))',
                                  rule_block, re.DOTALL):
                input_block = im.group(1)
                for ref in re.finditer(r'rules\.(\w+)\.output', input_block):
                    dep_rule = ref.group(1)
                    dep_id = f"rule__{dep_rule}"
                    edges.append({
                        "source": rule_id, "target": dep_id,
                        "type": "depends_on", "label": f"requires {dep_rule}",
                    })

    return edges


def parse_snakefile():
    edges = []
    snakefile = PROJECT_ROOT / "Snakefile"
    if not snakefile.exists():
        return edges
    content = snakefile.read_text()
    for m in re.finditer(r'include:\s*"workflow/rules/(.+?)"', content):
        rule_file = m.group(1)
        rule_id = f"workflow__rules__{slug(rule_file)}"
        edges.append({
            "source": "snakefile_main", "target": rule_id,
            "type": "includes", "label": "includes rule",
        })
    return edges


def build_graph():
    nodes = []
    edges = []

    # ── Root nodes ──
    root_nodes = [
        ("pipeline_root", "Pipeline Core", "root", "pipeline",
         "Snakemake workflow: alignment → filtering → HLA typing → contigs → translation → proteome filter → MHC affinity → TCRdock → report"),
        ("research_root", "Research Program", "root", "research",
         "Manuscript, tool evaluations, experiments, decisions, notebooks, slides"),
        ("infrastructure_root", "Infrastructure & Ops", "root", "infrastructure",
         "GCP orchestration, CI toolchain, release polling, Docker, project management scripts"),
        ("docs_root", "Documentation", "root", "docs",
         "User guides, feature decks, architecture notes, superpower specs"),
        ("config_root", "Configuration", "root", "config",
         "Pipeline configs, sample manifests, test configs, GPU overlay"),
    ]
    for nid, label, ntype, group, desc in root_nodes:
        nodes.append({
            "id": nid, "label": label, "type": ntype, "group": group,
            "path": "", "description": desc, "children": [], "size": 5,
        })

    # ── Root containment ──
    # Top-level dirs are nested under their group root by the directory walk
    # (each dir → its parent dir, top-level dir → GROUP_ROOT_ID[group]), so no
    # hardcoded root→subdir list is needed here — it produced an inconsistent,
    # partial tree (some groups pointed at subdirs, skipping the top-level dir).

    # ── Snakefile ──
    nodes.append({
        "id": "snakefile_main", "label": "Snakefile", "type": "rule",
        "group": "pipeline", "path": "Snakefile",
        "description": "Main workflow entry point — includes all rule modules",
        "size": 4,
    })
    edges.append({"source": "pipeline_root", "target": "snakefile_main", "type": "contained_in", "label": ""})

    # ── Pipeline individual rule nodes ──
    # Pre-parse .smk files for resource metadata so we can attach it to rule nodes.
    _rule_resources = {}
    for smk_path in sorted((PROJECT_ROOT / "workflow" / "rules").glob("*.smk")):
        content = smk_path.read_text()
        # Collect all rule positions first (including indented rules inside if/elif blocks)
        _entries = []
        for m in re.finditer(r'(?:^|\n)\s*rule\s+(\w+):', content):
            _entries.append((m.group(1), m.end()))
        for i, (rn, rule_end) in enumerate(_entries):
            next_pos = _entries[i + 1][1] + 1 if i + 1 < len(_entries) else len(content)
            rule_block = content[rule_end:next_pos]
            _rule_resources[rn] = parse_rule_resources(rule_block)

    for rule_name, meta in PIPELINE_RULES_DETAIL.items():
        node_id = f"rule__{rule_name}"
        smk_file = f"workflow/rules/{meta['smk']}"
        res = _rule_resources.get(rule_name, {})
        node_data = {
            "id": node_id, "label": rule_name, "type": "rule",
            "group": "pipeline", "path": smk_file,
            "description": meta["desc"],
            "detail_type": meta["type"],
            "outputs": meta["outputs"],
            "size": 3,
        }
        if res.get("threads"):
            node_data["threads"] = res["threads"]
        if res.get("mem_mb"):
            node_data["mem_mb"] = res["mem_mb"]
        if res.get("gpu"):
            node_data["gpu"] = True
        if res.get("docker"):
            node_data["docker"] = res["docker"]
        if res.get("conditional"):
            node_data["conditional"] = True
        if res.get("config_keys"):
            node_data["config_keys"] = res["config_keys"]
        nodes.append(node_data)
        edges.append({"source": node_id, "target": slug(smk_file), "type": "contained_in", "label": ""})

    # ── Pipeline script nodes with function-level metadata ──
    for script_path, meta in PIPELINE_SCRIPTS_META.items():
        script_id = slug(script_path)
        functions, imports, argparse_args, loc = get_script_metadata(script_path)
        nodes.append({
            "id": script_id, "label": os.path.basename(script_path), "type": "script",
            "group": "pipeline", "path": script_path,
            "description": meta["desc"],
            "inputs_schema": meta["inputs"],
            "outputs_schema": meta["outputs"],
            "functions": functions[:8],  # top functions
            "function_count": len(functions),
            "imports": imports[:12],
            "import_count": len(imports),
            "argparse_args": argparse_args[:6],
            "lines_of_code": loc,
            "size": 2,
        })
        # Nest under the workflow/scripts dir node (created by the walk).
        edges.append({"source": script_id, "target": "workflow__scripts",
                      "type": "contained_in", "label": ""})

    # ── Pipeline env nodes with package detail ──
    for env_path, meta in PIPELINE_ENVS_DETAIL.items():
        env_id = slug(env_path)
        nodes.append({
            "id": env_id, "label": os.path.basename(env_path).replace(".yaml", ""),
            "type": "env", "group": "pipeline", "path": env_path,
            "description": f"Conda env: {meta['channel']}",
            "packages": meta["packages"],
            "notes": meta.get("notes", ""),
            "size": 2,
        })
        # Nest under the workflow/envs dir node (created by the walk).
        edges.append({"source": env_id, "target": "workflow__envs",
                      "type": "contained_in", "label": ""})

    # ── Pipeline rule-file nodes (.smk) ──
    # The rule sub-nodes (`contained_in`), the Snakefile `include:` edges, and the
    # per-rule `defined_in` edges all point at these per-file nodes. The directory
    # walk skips workflow/rules/*.smk (the is_already_added rule guard), so without
    # this loop those edges dangle and D3's forceLink throws `missing: <id>`,
    # blanking the entire graph in the default "full" view.
    rules_dir = PROJECT_ROOT / "workflow" / "rules"
    if rules_dir.exists():
        for smk_file in sorted(rules_dir.glob("*.smk")):
            rel = str(smk_file.relative_to(PROJECT_ROOT))
            smk_id = slug(rel)
            nodes.append({
                "id": smk_id, "label": smk_file.name, "type": "rule",
                "group": "pipeline", "path": rel,
                "description": "Snakemake rule module", "size": 3,
            })
            edges.append({"source": smk_id, "target": "workflow__rules",
                          "type": "contained_in", "label": ""})

    # ── Walk remaining files ──
    # Seed the dedup map with every explicitly-added node id so the walk never
    # re-adds an already-detailed node under a different classification. (The env
    # files were being duplicated because classify_path mis-types them as
    # "resource" — its substring checks expect a leading slash the relative path
    # lacks — so the `ftype == "env"` dedup guard never fired.)
    file_nodes = {n["id"]: True for n in nodes}
    # Restrict discovery to the git-tracked tree so gitignored runtime/data dirs
    # never enter the atlas — the map is then identical between a clean checkout
    # and a clone that has run the pipeline (#780).
    tracked_files, tracked_dirs = git_tracked_paths()
    for root, dirs, files in os.walk(PROJECT_ROOT):
        # Prune build artifacts: __pycache__, VCS/env dirs, and Quarto's rendered
        # slide output (`*_files/` holds a vendored reveal.js copy — pure noise; the
        # canonical deck is slides.qmd).
        dirs[:] = [d for d in dirs if not d.startswith(".") and d not in (
            "__pycache__", ".snakemake", ".git", ".venv", "node_modules")
            and not d.endswith("_files")]
        rel_root = os.path.relpath(root, PROJECT_ROOT)
        # Descend only into dirs that contain a tracked file (drops references/,
        # results/, logs/, data/, indices/ and any other gitignored dir).
        dirs[:] = [d for d in dirs
                   if (d if rel_root == "." else f"{rel_root}/{d}") in tracked_dirs]
        if rel_root == ".":
            continue
        if rel_root in SELF_EXCLUDE_DIRS:
            dirs[:] = []  # don't descend into vendored deps
            continue
        # Keep only tracked files, so the dir node's `size` and the file loop are
        # both reproducible w.r.t. git rather than the local working tree.
        files = [f for f in files
                 if os.path.join(rel_root, f).replace(os.sep, "/") in tracked_files]
        dir_id = slug(rel_root)
        if dir_id in file_nodes:
            continue
        file_nodes[dir_id] = True
        ftype = classify_path(rel_root, is_dir=True)
        group = assign_group(rel_root)
        if ftype != "dir" or group != "project":
            _issue = issue_of(rel_root)
            nodes.append({
                "id": dir_id, "label": os.path.basename(rel_root),
                "type": ftype, "group": group, "path": rel_root,
                "size": len(files),
                **({"issue": _issue} if _issue else {}),
            })
            # Containment: nest this dir under its parent dir, or — for a
            # top-level dir — under its group root, completing the tree.
            parent_rel = os.path.dirname(rel_root)
            parent_id = slug(parent_rel) if parent_rel else GROUP_ROOT_ID.get(group)
            if parent_id:
                edges.append({"source": dir_id, "target": parent_id,
                              "type": "contained_in", "label": ""})
        for f in files:
            if f.endswith((".pyc", ".pyo")):
                continue
            rel_path = os.path.join(rel_root, f)
            if rel_path in SELF_EXCLUDE_FILES:
                continue
            file_id = slug(rel_path)
            if file_id in file_nodes:
                continue
            file_nodes[file_id] = True
            ftype = classify_path(rel_path, is_dir=False)
            if ftype == "ignore":
                continue
            group = assign_group(rel_path)
            # Skip if we already added it as a deep pipeline node
            is_already_added = (
                (ftype == "script" and rel_path in PIPELINE_SCRIPTS_META) or
                (ftype == "env" and rel_path in PIPELINE_ENVS_DETAIL) or
                (ftype == "rule" and rel_path.startswith("workflow/rules/"))
            )
            if is_already_added:
                continue

            sz = 1
            if ftype in ("rule", "env"): sz = 3
            elif ftype in ("script", "slide", "manuscript", "notebook"): sz = 2
            elif ftype == "ci": sz = 1

            _issue = issue_of(rel_path)
            nodes.append({
                "id": file_id, "label": f, "type": ftype,
                "group": group, "path": rel_path,
                "description": describe_file(rel_path, ftype), "size": sz,
                **({"issue": _issue} if _issue else {}),
            })
            edges.append({"source": file_id, "target": dir_id, "type": "contained_in", "label": ""})

    # ── Parse Snakemake relationships ──
    smk_edges = parse_smk_rules()
    edges.extend(smk_edges)

    snakefile_edges = parse_snakefile()
    edges.extend(snakefile_edges)

    # ── Pipeline dataflow edges ──
    pipeline_dataflow = [
        ("rule__build_reference_junctions", "rule__filter_junctions", "produces_reference"),
        ("rule__filter_junctions", "rule__assemble_contigs", "produces_novel_junctions"),
        ("rule__assemble_contigs", "rule__translate_peptides", "produces_contigs"),
        ("rule__translate_peptides", "rule__proteome_filter_peptides", "produces_peptides"),
        ("rule__proteome_filter_peptides", "rule__run_mhcflurry", "produces_novel_peptides"),
        ("rule__run_mhcflurry", "rule__run_tcrdock", "produces_predictions"),
        ("rule__run_mhcflurry", "rule__generate_report", "produces_predictions"),
        ("rule__aggregate_filtering_stats", "rule__generate_report", "produces_filtering_stats"),
        ("rule__run_optitype", "rule__aggregate_hla_alleles", "produces_optitype_calls"),
        ("rule__aggregate_hla_alleles", "rule__run_mhcflurry", "produces_alleles"),
        ("rule__aggregate_hla_alleles", "rule__fetch_vdjdb_panel", "produces_alleles_for_tcr"),
        ("rule__fetch_vdjdb_panel", "rule__run_tcrdock", "produces_tcr_panel"),
        ("rule__download_vdjdb_release", "rule__fetch_vdjdb_panel", "produces_vdjdb_data"),
        ("rule__download_imgt_germlines", "rule__fetch_vdjdb_panel", "produces_imgt_data"),
        ("rule__download_human_proteome", "rule__proteome_filter_peptides", "produces_proteome_fasta"),
    ]
    for src, tgt, lbl in pipeline_dataflow:
        edges.append({"source": src, "target": tgt, "type": "produces", "label": lbl})

    # ── Research reference edges ──
    evals_dir = PROJECT_ROOT / "research" / "evals"
    if evals_dir.exists():
        for eval_dir in sorted(evals_dir.iterdir()):
            if eval_dir.is_dir():
                slide_file = eval_dir / "slides.qmd"
                if slide_file.exists():
                    eval_id = slug(f"research__evals__{eval_dir.name}")
                    slide_id = slug(f"research__evals__{eval_dir.name}__slides_qmd")
                    edges.append({"source": eval_id, "target": slide_id, "type": "references", "label": "slide deck"})
    exp_dir = PROJECT_ROOT / "research" / "experiments"
    if exp_dir.exists():
        for exp_d in sorted(exp_dir.iterdir()):
            if exp_d.is_dir():
                slide_file = exp_d / "slides.qmd"
                if slide_file.exists():
                    exp_id = slug(f"research__experiments__{exp_d.name}")
                    slide_id = slug(f"research__experiments__{exp_d.name}__slides_qmd")
                    edges.append({"source": exp_id, "target": slide_id, "type": "references", "label": "slide deck"})

    # ── Test coverage edges ──
    test_dir = PROJECT_ROOT / "workflow" / "tests"
    if test_dir.exists():
        for tf in sorted(test_dir.glob("test_*.py")):
            test_stem = tf.stem.replace("test_", "", 1)
            test_id = slug(f"workflow__tests__{tf.name}")
            for script_file in (PROJECT_ROOT / "workflow" / "scripts").glob(f"{test_stem}.py"):
                script_id = slug(f"workflow/scripts/{script_file.name}")
                edges.append({"source": test_id, "target": script_id, "type": "test_of", "label": "tests"})

    # ── Research/docs cross-links (markdown → markdown) ──
    # Turns ~157 prose links into real `references` edges so the research half
    # stops being a disconnected file dump. The dangling-edge filter below drops
    # any link whose target isn't a node, so this is safe.
    edges.extend(parse_markdown_links({n["id"] for n in nodes}))

    # ── Deduplicate nodes by id (keep first) ──
    # Belt-and-suspenders against any id collision the explicit/walk node sources
    # might reintroduce; duplicate ids inflate the count and confuse D3's keying.
    seen_node_ids = set()
    unique_nodes = []
    for n in nodes:
        if n["id"] in seen_node_ids:
            continue
        seen_node_ids.add(n["id"])
        unique_nodes.append(n)
    dropped_nodes = len(nodes) - len(unique_nodes)
    nodes = unique_nodes

    # ── Drop edges with unknown endpoints ──
    # D3's forceLink throws `missing: <id>` on the first edge whose source or
    # target has no node, which blanks the whole graph. Guarantee a render-safe
    # graph by dropping dangling edges here — and report the count, never silently.
    valid_edges = []
    dropped_edges = []
    for e in edges:
        if e["source"] in seen_node_ids and e["target"] in seen_node_ids:
            valid_edges.append(e)
        else:
            dropped_edges.append(e)
    edges = valid_edges

    # ── Deduplicate edges ──
    seen_edges = set()
    unique_edges = []
    for e in edges:
        key = f"{e['source']}|{e['target']}|{e['type']}"
        if key not in seen_edges:
            seen_edges.add(key)
            unique_edges.append(e)

    # ── Degrees ──
    # Tally from unique_edges, NOT edges: a duplicated edge would otherwise
    # double-count both endpoints, inflating degree → sizeOf() renders those
    # nodes artificially large in the HTML.
    degree_count = {}
    for e in unique_edges:
        degree_count[e["source"]] = degree_count.get(e["source"], 0) + 1
        degree_count[e["target"]] = degree_count.get(e["target"], 0) + 1
    for n in nodes:
        n["degree"] = degree_count.get(n["id"], 0)

    if dropped_nodes:
        print(f"  ⚠ dropped {dropped_nodes} duplicate node id(s)")
    if dropped_edges:
        sample = ", ".join(sorted({f"{e['source']}->{e['target']}" for e in dropped_edges})[:5])
        print(f"  ⚠ dropped {len(dropped_edges)} dangling edge(s): {sample}")

    return {"nodes": nodes, "edges": unique_edges}


def main():
    import argparse
    parser = argparse.ArgumentParser(description="Extract project structure as JSON graph")
    parser.add_argument("--output", default=str(HERE / "graph.json"))
    args = parser.parse_args()

    graph = build_graph()
    output_path = Path(args.output)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(graph, indent=2))
    print(f"Wrote {len(graph['nodes'])} nodes and {len(graph['edges'])} edges to {output_path}")

    # Stats
    pipeline_nodes = [n for n in graph['nodes'] if n['group'] == 'pipeline']
    pipeline_edges = [e for e in graph['edges'] if any(
        n['id'] in (e['source'], e['target']) for n in pipeline_nodes)]
    print(f"  Pipeline: {len(pipeline_nodes)} nodes, {len(pipeline_edges)} edges")
    rules = [n for n in pipeline_nodes if n.get('detail_type')]
    print(f"  Individual rules with metadata: {len(rules)}")
    scripts = [n for n in pipeline_nodes if n['type'] == 'script' and 'functions' in n]
    print(f"  Scripts with function/import detail: {len(scripts)}")
    res = sum(1 for n in graph["nodes"] if n["type"] == "resource")
    refs = sum(1 for e in graph["edges"] if e["type"] == "references")
    issued = sum(1 for n in graph["nodes"] if n.get("issue"))
    print(f"  Unclassified 'resource' nodes: {res}/{len(graph['nodes'])} "
          f"({100*res//len(graph['nodes'])}%)  |  references edges: {refs}  |  issue-tagged nodes: {issued}")

if __name__ == "__main__":
    main()