#!/usr/bin/env Rscript
# Issue #965 leaf A - splice2neo open-only local smoke (macOS arm64, CPU).
# A "true smoke, not full scale": exercise the junction -> context-sequence ->
# peptide path on the package's OWN bundled toy fixtures, so no external caller
# or patient data is needed. Proves the MIT/arm64 install runs end to end and
# emits the junction-derived peptide sequences that would feed our MHCflurry stage.
# API verified against splice2neo 0.6.14:
#   add_context_seq(df, transcripts, size=400, bsg)  -> adds `cts_seq`
#   add_peptide(df, cds, flanking_size=14, bsg)      -> adds `protein`, `peptide_context`, `frame_shift`
# Run from the repo root (the outdir path below is repo-root-relative).

suppressMessages({
  library(splice2neo)
  library(BSgenome.Hsapiens.UCSC.hg19)
})
bsg <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19

outdir <- "research/experiments/issue_965_open_caller_audit/outputs"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

cat("splice2neo", as.character(packageVersion("splice2neo")),
    "| R", as.character(getRversion()),
    "| hg19", as.character(packageVersion("BSgenome.Hsapiens.UCSC.hg19")), "\n")

# Bundled toy fixtures: toy_junc_df (junc_id + tx_id), toy_transcripts, toy_cds.
stopifnot(exists("toy_junc_df"), exists("toy_transcripts"), exists("toy_cds"))
cat("toy_junc_df:", nrow(toy_junc_df), "junctions x", ncol(toy_junc_df), "cols\n")

# Count only real, non-empty strings. `nzchar(NA)` returns TRUE by default
# (keepNA = FALSE), so an NA cell would be miscounted as present - guard with !is.na.
non_empty <- function(x) sum(!is.na(x) & nzchar(as.character(x)))

# 1) Junction -> transcript context sequence.
ctx <- add_context_seq(toy_junc_df, transcripts = toy_transcripts, size = 400, bsg = bsg)
n_ctx <- non_empty(ctx$cts_seq)
cat("[1] add_context_seq ->", nrow(ctx), "rows;", n_ctx, "with a context sequence (cts_seq)\n")

# 2) Context -> mutated CDS -> peptide (the artifact that would feed MHCflurry).
pep <- add_peptide(ctx, cds = toy_cds, flanking_size = 14, bsg = bsg)
n_prot <- non_empty(pep$protein)
n_pctx <- non_empty(pep$peptide_context)
cat("[2] add_peptide ->", nrow(pep), "rows;", n_prot, "with a mutated protein;",
    n_pctx, "with a junction peptide_context;", sum(pep$frame_shift, na.rm = TRUE),
    "frame-shift junctions\n")

# Persist a small, human-readable artifact (drop the huge full-protein column).
out <- data.frame(
  junc_id          = pep$junc_id,
  tx_id            = pep$tx_id,
  frame_shift      = pep$frame_shift,
  cts_seq_len      = nchar(as.character(pep$cts_seq)),
  peptide_context  = as.character(pep$peptide_context),
  stringsAsFactors = FALSE
)
utils::write.table(head(out, 20), file.path(outdir, "splice2neo_smoke_out.tsv"),
                   sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nWrote", file.path(outdir, "splice2neo_smoke_out.tsv"),
    "(", min(20, nrow(out)), "rows )\n")
stopifnot(n_ctx > 0, n_prot > 0)
cat("SMOKE_OK\n")
