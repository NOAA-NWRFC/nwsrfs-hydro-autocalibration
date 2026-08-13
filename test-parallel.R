#!/usr/bin/env Rscript

# Checks the parallel backend that ep_dds() in R/wrappers.R and the bootstrap in
# cv-plots.R rely on. Both pick FORK on Linux and macOS and PSOCK on Windows,
# then bootstrap each worker by hand: set the working directory, load the local
# box modules, and export the wrapper functions. A PSOCK worker starts from an
# empty session, so anything the master picked up implicitly has to be loaded
# again on the worker. That is the part that breaks on Windows, and this script
# exercises it without waiting for a full calibration.
#
# Run it from the repository root:
#
#   Rscript test-parallel.R
#
# PSOCK is the strict case, so it is worth running on Linux and macOS too rather
# than waiting for a Windows machine to disagree. Pass the type to force it:
#
#   Rscript test-parallel.R PSOCK
#
# It lives beside run-controller.R rather than under R/ on purpose. box resolves
# a relative module against the directory of the script that asks for it, so
# `./R/metrics` only resolves for a caller at the repository root, which is where
# the worker bootstrap being tested runs from.
#
# Written by Cameron Bracken and Geoffrey Walters (2025)
# Please see the LICENSE file for license information

box::use(
  parallel[clusterCall, clusterEvalQ, clusterExport, clusterSetRNGStream,
           detectCores, makeCluster, stopCluster],
  ./R/metrics[NSE]
)

failures <- 0
checks <- 0

check <- function(label, ok, detail = "") {
  checks <<- checks + 1
  if (!isTRUE(ok)) failures <<- failures + 1
  cat(sprintf("  %-46s %s%s\n", label, if (isTRUE(ok)) "ok" else "FAIL",
              if (nzchar(detail)) paste0("  (", detail, ")") else ""))
}

cat("platform: ", R.version$platform, " (", .Platform$OS.type, ")\n", sep = "")
cat("cores detected: ", detectCores(), "\n\n", sep = "")

# The same choice ep_dds() and cv-plots.R make.
expected <- if (.Platform$OS.type == "windows") "PSOCK" else "FORK"

forced <- commandArgs(trailingOnly = TRUE)[1]
if (!is.na(forced)) {
  forced <- toupper(forced)
  if (!forced %in% c("FORK", "PSOCK")) stop("cluster type must be FORK or PSOCK")
  if (forced == "FORK" && .Platform$OS.type == "windows") {
    stop("Windows has no FORK cluster, nothing to force")
  }
  os_type <- forced
  cat("forcing ", os_type, " (platform default is ", expected, ")\n\n", sep = "")
} else {
  os_type <- expected
  check(paste0("cluster type is ", expected), identical(os_type, expected), os_type)
}

# Two workers is enough to tell a per-worker problem from a shared one, and it
# fits the smallest runner.
n_cores <- max(2, min(2, detectCores()))

cat("\n=== worker startup ===\n")
my_cluster <- makeCluster(n_cores, type = os_type)
on.exit(try(stopCluster(my_cluster), silent = TRUE), add = TRUE)
check("cluster started", inherits(my_cluster, "cluster"),
      paste(n_cores, os_type, "workers"))

RNGkind("L'Ecuyer-CMRG")
clusterSetRNGStream(cl = my_cluster, iseed = NULL)

# Byte for byte the bootstrap ep_dds() performs, so a failure here is a failure
# there. A FORK worker inherits all of this and a PSOCK worker does not.
master_wd <- getwd()
clusterExport(my_cluster, "master_wd", envir = environment())
boot <- try(
  clusterEvalQ(my_cluster, {
    setwd(master_wd)
    box::use(
      stats[runif, rnorm, setNames],
      dplyr[filter, select, summarise, group_by, ungroup, mutate,
            lead, right_join, bind_rows],
      data.table[as.data.table, data.table, merge.data.table, copy,
                 rbindlist, nafill],
      tibble[as_tibble],
      tidyr[fill],
      ./R/metrics[NSE, pbias, rPearson, KGE],
      nwsrfsr[sac_snow_uh, sac_snow_uh_lagk, lagk, chanloss,
              consuse, fa_nwrfc],
      obj_funs = ./R/obj_fun
    )
    TRUE
  }),
  silent = TRUE
)
check("workers loaded local and package modules", !inherits(boot, "try-error"),
      if (inherits(boot, "try-error")) trimws(as.character(boot)) else "")

# The working directory the modules were resolved against.
wds <- unlist(clusterCall(my_cluster, getwd))
check("workers share the master working directory",
      all(normalizePath(wds) == normalizePath(master_wd)))

cat("\n=== worker execution ===\n")
# A function defined here and shipped to the workers. Under PSOCK it arrives
# without its enclosing environment and resolves `NSE` against the worker's
# global environment, which is what the box::use above populated.
worker_task <- function() {
  set.seed(NULL)
  obs <- c(10, 12, 9, 14, 11, 13, 8, 15)
  sim <- obs * 1.1
  list(pid = Sys.getpid(), seed = globalenv()$.Random.seed, nse = NSE(sim, obs))
}
res <- try(clusterCall(my_cluster, worker_task), silent = TRUE)
check("workers ran a local module function", !inherits(res, "try-error"),
      if (inherits(res, "try-error")) trimws(as.character(res)) else "")

if (!inherits(res, "try-error")) {
  worker_nse <- sapply(res, function(x) x$nse)
  master_nse <- NSE(c(10, 12, 9, 14, 11, 13, 8, 15) * 1.1,
                    c(10, 12, 9, 14, 11, 13, 8, 15))
  check("worker results match the master",
        isTRUE(all.equal(unname(worker_nse), rep(master_nse, length(worker_nse)))),
        sprintf("%.6f", master_nse))

  # run-controller.R aborts the whole calibration when two workers report the
  # same stream, so check it here where the message can be useful.
  seeds <- sapply(res, function(x) paste(x$seed, collapse = ","))
  check("each worker has a distinct RNG stream",
        length(unique(seeds)) == length(seeds),
        sprintf("%d streams across %d workers", length(unique(seeds)), length(seeds)))

  kinds <- unique(sapply(res, function(x) x$seed[1]))
  check("workers use L'Ecuyer-CMRG", all(kinds %% 100 == 7),
        paste(kinds, collapse = ", "))
}

cat("\n=== teardown ===\n")
stopped <- try(stopCluster(my_cluster), silent = TRUE)
check("cluster stopped", !inherits(stopped, "try-error"))

cat(sprintf("\n%d checks, %d failures\n", checks, failures))
if (failures > 0) quit(status = 1)
