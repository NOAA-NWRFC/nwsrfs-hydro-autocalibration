#!/usr/bin/env Rscript

# Checks R/metrics.R against hydroGOF, the reference implementation of the same
# published metrics. hydroGOF is GPL licensed and is NOT a dependency of this
# project; it is only needed to run this check. Install it first with:
#
#   Rscript -e "install.packages('hydroGOF')"
#   ./R/test-metrics.R
#
# Written by Cameron Bracken and Geoffrey Walters (2025)
# Please see the LICENSE file for license information

box::use(
  ./metrics
)

if (!requireNamespace("hydroGOF", quietly = TRUE)) {
  stop("hydroGOF is needed to run this check, install it with install.packages('hydroGOF')")
}
hg <- asNamespace("hydroGOF")

tol <- 1e-9
failures <- 0
checks <- 0

check <- function(label, mine, theirs) {
  checks <<- checks + 1
  ok <- isTRUE(all.equal(mine, theirs, tolerance = tol))
  if (!ok) failures <<- failures + 1
  cat(sprintf(
    "  %-28s %-20.12g %-20.12g %s\n", label, mine, theirs,
    if (ok) "ok" else "MISMATCH"
  ))
}

# Several regimes, so that the checks cover skew, near-zero flows and negative
# correlation rather than one comfortable case.
datasets <- list(
  gaussian = function() {
    o <- abs(rnorm(500, 100, 40)) + 1
    list(obs = o, sim = pmax(o * runif(500, 0.6, 1.4) + rnorm(500, 2, 5), 0.1))
  },
  skewed = function() {
    o <- rlnorm(400, 4, 1.2) + 0.5
    list(obs = o, sim = pmax(o^1.05 * runif(400, 0.5, 1.5), 0.05))
  },
  low_flow = function() {
    o <- rlnorm(300, 0, 1.5) + 0.01
    list(obs = o, sim = pmax(o * runif(300, 0.3, 2.0), 0.001))
  },
  poor_fit = function() {
    o <- abs(rnorm(250, 50, 20)) + 1
    list(obs = o, sim = pmax(rev(o) * runif(250, 0.5, 1.5), 0.1))
  },
  near_perfect = function() {
    o <- abs(rnorm(200, 80, 30)) + 1
    list(obs = o, sim = o * runif(200, 0.98, 1.02))
  }
)

set.seed(2025)
for (nm in names(datasets)) {
  dat <- datasets[[nm]]()
  s <- dat$sim
  o <- dat$obs
  cat("\n=== dataset:", nm, "===\n")

  check("me", metrics$me(s, o), hg$me(s, o))
  check("mae", metrics$mae(s, o), hg$mae(s, o))
  check("mse", metrics$mse(s, o), hg$mse(s, o))
  check("rmse", metrics$rmse(s, o), hg$rmse(s, o))
  check("ubRMSE", metrics$ubRMSE(s, o), hg$ubRMSE(s, o))
  check("nrmse sd", metrics$nrmse(s, o), hg$nrmse(s, o, norm = "sd"))
  check("nrmse maxmin", metrics$nrmse(s, o, norm = "maxmin"), hg$nrmse(s, o, norm = "maxmin"))
  check("pbias", metrics$pbias(s, o), hg$pbias(s, o))
  check("rsr", metrics$rsr(s, o), hg$rsr(s, o))
  check("rSD", metrics$rSD(s, o), hg$rSD(s, o))
  check("NSE", metrics$NSE(s, o), hg$NSE(s, o))
  check("mNSE j=1", metrics$mNSE(s, o, j = 1), hg$mNSE(s, o, j = 1))
  check("mNSE j=2", metrics$mNSE(s, o, j = 2), hg$mNSE(s, o, j = 2))
  check("rNSE", metrics$rNSE(s, o), hg$rNSE(s, o))
  check("wNSE", metrics$wNSE(s, o), hg$wNSE(s, o))
  check("wsNSE default", metrics$wsNSE(s, o), hg$wsNSE(s, o))
  check("wsNSE j=1", metrics$wsNSE(s, o, j = 1), hg$wsNSE(s, o, j = 1))
  check("d", metrics$d(s, o), hg$d(s, o))
  check("dr", metrics$dr(s, o), hg$dr(s, o))
  check("md j=1", metrics$md(s, o, j = 1), hg$md(s, o, j = 1))
  check("rd", metrics$rd(s, o), hg$rd(s, o))
  check("cp", metrics$cp(s, o), hg$cp(s, o))
  check("rPearson", metrics$rPearson(s, o), hg$rPearson(s, o))
  check("rSpearman", metrics$rSpearman(s, o), hg$rSpearman(s, o))
  check("R2", metrics$R2(s, o), hg$R2(s, o))
  check("br2", metrics$br2(s, o), hg$br2(s, o))
  check("VE", metrics$VE(s, o), hg$VE(s, o))
  check("KGE 2009", metrics$KGE(s, o, method = "2009"), hg$KGE(s, o, method = "2009"))
  check("KGE 2012", metrics$KGE(s, o, method = "2012"), hg$KGE(s, o, method = "2012"))
  check("KGE 2021", metrics$KGE(s, o, method = "2021"), hg$KGE(s, o, method = "2021"))
  check(
    "KGE weighted", metrics$KGE(s, o, s = c(0.5, 0.3, 0.2)),
    hg$KGE(s, o, s = c(0.5, 0.3, 0.2))
  )
  check("KGElf 2009", metrics$KGElf(s, o, method = "2009"), hg$KGElf(s, o, method = "2009"))
  check("KGElf 2012", metrics$KGElf(s, o, method = "2012"), hg$KGElf(s, o, method = "2012"))
  check(
    "KGElf eps none", metrics$KGElf(s, o, epsilon.type = "none"),
    hg$KGElf(s, o, epsilon.type = "none")
  )
  check("KGEnp", metrics$KGEnp(s, o), hg$KGEnp(s, o))
  check("KGEkm 2009", metrics$KGEkm(s, o, method = "2009"), hg$KGEkm(s, o, method = "2009"))
  check("KGEkm 2012", metrics$KGEkm(s, o, method = "2012"), hg$KGEkm(s, o, method = "2012"))
  check("LME", metrics$LME(s, o), hg$LME(s, o))
  check("LCE", metrics$LCE(s, o), hg$LCE(s, o))

  # The whole summary table, row by row.
  mine <- metrics$gof(s, o)
  theirs <- hg$gof(s, o)
  shared <- intersect(rownames(mine), rownames(theirs))
  cat("  -- gof table --\n")
  for (rn in shared) check(paste("gof", rn), mine[rn, 1], theirs[rn, 1])
  extra <- setdiff(rownames(theirs), rownames(mine))
  if (length(extra)) {
    cat("  not reproduced:", paste(extra, collapse = ", "), "\n")
  }
}

# NA handling: pairs where either series is missing must drop out.
cat("\n=== NA handling ===\n")
set.seed(11)
o <- abs(rnorm(200, 100, 30)) + 1
s <- pmax(o * runif(200, 0.7, 1.3), 0.1)
o[c(3, 17, 45)] <- NA
s[c(9, 17, 88)] <- NA
check("NSE with NA", metrics$NSE(s, o), hg$NSE(s, o))
check("KGE with NA", metrics$KGE(s, o, method = "2012"), hg$KGE(s, o, method = "2012"))
check("pbias with NA", metrics$pbias(s, o), hg$pbias(s, o))
check("cp with NA", metrics$cp(s, o), hg$cp(s, o))
check("gof NSE with NA", metrics$gof(s, o)["NSE", 1], hg$gof(s, o)["NSE", 1])

cat(sprintf("\n%d checks, %d mismatches\n", checks, failures))
if (failures > 0) quit(status = 1)
