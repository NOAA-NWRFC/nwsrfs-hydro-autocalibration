# Worker cluster construction, shared by the EDDS optimizer in wrappers.R and
# the cross validation bootstrap in cv-plots.R.
#
# Written by Cameron Bracken and Geoffrey Walters (2025)
# Please see the LICENSE file for license information

box::use(
  parallel[makeCluster]
)

#' Build the worker cluster.
#'
#' FORK on Linux and macOS, which is cheaper because workers inherit the
#' master's memory, and PSOCK on Windows, which has no fork.
#'
#' Since R 4.0 makePSOCKcluster starts every worker at once by default. It is
#' faster, but a worker occasionally fails to connect, and the master then sits
#' out the full connection timeout, several minutes, before giving up and taking
#' the calibration with it. Starting the workers one at a time avoids that. The
#' extra cost is roughly a second at the start of a run that takes hours.
#'
#' @param n_cores number of workers to start
#' @param type cluster type, defaults to the right one for the platform. Only
#'   set it to check the PSOCK path on a machine that would otherwise fork.
#' @export
make_worker_cluster <- function(n_cores, type = NULL) {
  if (is.null(type)) {
    type <- if (.Platform$OS.type == "windows") "PSOCK" else "FORK"
  }
  if (identical(type, "PSOCK")) {
    makeCluster(n_cores, type = "PSOCK", setup_strategy = "sequential")
  } else {
    makeCluster(n_cores, type = type)
  }
}
