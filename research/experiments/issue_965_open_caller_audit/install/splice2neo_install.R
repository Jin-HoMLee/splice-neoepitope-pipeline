#!/usr/bin/env Rscript
# Issue #965 leaf A - splice2neo open-only install recipe (macOS arm64, CPU).
# splice2neo = TRON-Bioinformatics/splice2neo, MIT, R>=4.4, junction->peptide library.
# Installs the package + its Bioconductor deps, then the hg19 BSgenome needed by the
# toy context-sequence smoke. Idempotent: skips anything already installed.

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(timeout = 3600)  # BSgenome hg19 is ~700 MB

# macOS toolchain fix (homebrew R 4.6.0 / arm64): the default compile line omits
# -isysroot, so Apple clang can't find libc++ headers and every C++ source package
# dies with `fatal error: 'cmath' file not found`. Point the C++ include path at
# the active SDK's headers. (No CRAN/Bioconductor arm64 binaries exist yet for R
# 4.6.0 / Bioc 3.23, so the whole stack builds from source and needs this.)
sdk <- tryCatch(system2("xcrun", "--show-sdk-path", stdout = TRUE), error = function(e) "")
if (length(sdk) && nzchar(sdk) && dir.exists(file.path(sdk, "usr/include/c++/v1")))
  Sys.setenv(CPLUS_INCLUDE_PATH = paste(file.path(sdk, "usr/include/c++/v1"),
                                        file.path(sdk, "usr/include"), sep = ":"))

need <- function(pkg) !requireNamespace(pkg, quietly = TRUE)

if (need("BiocManager")) install.packages("BiocManager")
if (need("remotes"))     install.packages("remotes")

# Bioconductor genome dep used by add_context_seq() in the toy smoke.
if (need("BSgenome.Hsapiens.UCSC.hg19"))
  BiocManager::install("BSgenome.Hsapiens.UCSC.hg19", update = FALSE, ask = FALSE)

# The package itself (pulls its Bioconductor Imports as binaries on arm64).
if (need("splice2neo"))
  remotes::install_github("TRON-Bioinformatics/splice2neo", upgrade = "never")

cat("\n=== versions ===\n")
cat("R:", as.character(getRversion()), "\n")
cat("splice2neo:", as.character(packageVersion("splice2neo")), "\n")
cat("BSgenome.hg19:", as.character(packageVersion("BSgenome.Hsapiens.UCSC.hg19")), "\n")
cat("INSTALL_OK\n")
