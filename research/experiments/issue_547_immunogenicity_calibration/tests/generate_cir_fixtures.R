# Generates frozen centered-isotonic-regression oracle outputs from the R `cir`
# package for the Python port test. Run once: Rscript generate_cir_fixtures.R
# Requires: install.packages("cir"); jsonlite. Records package versions.
library(cir); library(jsonlite)
cases <- list(
  # Oron & Flournoy 2017 style: monotone-with-plateau dose-response
  onf_example = list(x=c(1,2,3,4,5,6), y=c(0.1,0.2,0.2,0.2,0.5,0.9), w=c(10,10,10,10,10,10)),
  multi_levelset = list(x=c(1,2,3,4,5,6,7,8), y=c(0.1,0.1,0.1,0.4,0.4,0.7,0.7,0.7), w=c(5,8,3,10,6,4,9,2)),
  boundary_plateau = list(x=c(1,2,3,4,5), y=c(0.3,0.3,0.3,0.6,0.9), w=c(7,7,7,7,7)),
  no_collapse = list(x=c(1,2,3,4), y=c(0.1,0.3,0.6,0.9), w=c(4,4,4,4)),
  # PAVA-pooling probe: raw y is NON-monotone (0.5 then 0.3), so PAVA must pool
  # the violator block before centering. Equal weights -> the x=2,3 violation
  # pools with x=4 into a y=0.4 plateau. Exercises the integrated pool+center
  # path that the Python port splits into sklearn-isotonic + centered_isotonic.
  pooling_violation = list(x=c(1,2,3,4,5), y=c(0.1,0.5,0.3,0.4,0.9), w=c(6,6,6,6,6)),
  # Asymmetric-weight pooling: violator block with unequal weights, so the
  # pooled plateau value is a weighted mean and the centroid x is weight-skewed.
  pooling_asym = list(x=c(1,2,3,4,5), y=c(0.2,0.6,0.3,0.5,0.8), w=c(3,9,2,4,5))
)
out <- list(r_version=R.version.string, cir_version=as.character(packageVersion("cir")), cases=list())
for (nm in names(cases)) {
  cs <- cases[[nm]]
  fit <- cirPAVA(y=cs$y, x=cs$x, wt=cs$w, full=TRUE)  # centered isotonic (CIR)
  # cir 2.5.1: full=TRUE returns list(output, input, shrinkage). `shrinkage`
  # holds the collapsed/centered knots (matching the Python port's
  # centered_isotonic() return); `output` is that fit interpolated back to all
  # input x. We compare against the centered knots -> use `shrinkage`.
  knots <- fit$shrinkage
  out$cases[[nm]] <- list(x=cs$x, y=cs$y, w=cs$w,
                          cir_x=as.numeric(knots$x), cir_y=as.numeric(knots$y))
}
write_json(out, "fixtures/cir_fixtures.json", auto_unbox=TRUE, digits=10, pretty=TRUE)
