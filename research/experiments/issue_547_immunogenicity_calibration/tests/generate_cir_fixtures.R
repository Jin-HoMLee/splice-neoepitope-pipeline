# Generates frozen centered-isotonic-regression oracle outputs from the R `cir`
# package for the Python port test. Run once: Rscript generate_cir_fixtures.R
# Requires: install.packages("cir"); jsonlite. Records package versions.
library(cir); library(jsonlite)
cases <- list(
  # Oron & Flournoy 2017 style: monotone-with-plateau dose-response
  onf_example = list(x=c(1,2,3,4,5,6), y=c(0.1,0.2,0.2,0.2,0.5,0.9), w=c(10,10,10,10,10,10)),
  multi_levelset = list(x=c(1,2,3,4,5,6,7,8), y=c(0.1,0.1,0.1,0.4,0.4,0.7,0.7,0.7), w=c(5,8,3,10,6,4,9,2)),
  boundary_plateau = list(x=c(1,2,3,4,5), y=c(0.3,0.3,0.3,0.6,0.9), w=c(7,7,7,7,7)),
  no_collapse = list(x=c(1,2,3,4), y=c(0.1,0.3,0.6,0.9), w=c(4,4,4,4))
)
out <- list(r_version=R.version.string, cir_version=as.character(packageVersion("cir")), cases=list())
for (nm in names(cases)) {
  cs <- cases[[nm]]
  fit <- cirPAVA(y=cs$y, x=cs$x, wt=cs$w, full=TRUE)  # centered isotonic (CIR)
  out$cases[[nm]] <- list(x=cs$x, y=cs$y, w=cs$w, cir_x=as.numeric(fit$x), cir_y=as.numeric(fit$y))
}
write_json(out, "fixtures/cir_fixtures.json", auto_unbox=TRUE, digits=10, pretty=TRUE)
