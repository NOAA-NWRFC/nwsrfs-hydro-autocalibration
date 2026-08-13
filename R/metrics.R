# Goodness-of-fit metrics for streamflow simulation.
#
# Written by Cameron Bracken and Geoffrey Walters (2025)
# Please see the LICENSE file for license information
#
# These functions are implemented from the published formulas cited against each
# metric below. They are not derived from, and contain no code from, any existing
# goodness-of-fit package. The R package hydroGOF (Zambrano-Bigiarini, GPL >= 2)
# implements the same set of published metrics and was used as an external
# reference to check these implementations numerically; every function here
# agrees with hydroGOF 0.7-0 to within 1e-9 (see the checks in `test-metrics.R`).
# hydroGOF is GPL licensed and this project is Apache 2.0, so no hydroGOF code
# could be vendored here.
#
# Argument order follows the hydrology convention throughout: metric(sim, obs).
#
# References
#
#   Nash, J.E. and Sutcliffe, J.V. (1970). River flow forecasting through
#     conceptual models part I: A discussion of principles. Journal of Hydrology
#     10(3), 282-290. doi:10.1016/0022-1694(70)90255-6
#   Willmott, C.J. (1981). On the validation of models. Physical Geography 2(2),
#     184-194. doi:10.1080/02723646.1981.10642213
#   Kitanidis, P.K. and Bras, R.L. (1980). Real-time forecasting with a conceptual
#     hydrologic model. Water Resources Research 16(6), 1025-1033.
#     doi:10.1029/WR016i006p01025
#   Legates, D.R. and McCabe, G.J. (1999). Evaluating the use of goodness-of-fit
#     measures in hydrologic and hydroclimatic model validation. Water Resources
#     Research 35(1), 233-241. doi:10.1029/1998WR900018
#   Krause, P., Boyle, D.P. and Base, F. (2005). Comparison of different efficiency
#     criteria for hydrological model assessment. Advances in Geosciences 5, 89-97.
#     doi:10.5194/adgeo-5-89-2005
#   Moriasi, D.N. et al. (2007). Model evaluation guidelines for systematic
#     quantification of accuracy in watershed simulations. Transactions of the
#     ASABE 50(3), 885-900. doi:10.13031/2013.23153
#   Criss, R.E. and Winston, W.E. (2008). Do Nash values have value? Discussion and
#     alternate proposals. Hydrological Processes 22(14), 2723-2725.
#     doi:10.1002/hyp.7072
#   Gupta, H.V., Kling, H., Yilmaz, K.K. and Martinez, G.F. (2009). Decomposition of
#     the mean squared error and NSE performance criteria: Implications for improving
#     hydrological modelling. Journal of Hydrology 377(1-2), 80-91.
#     doi:10.1016/j.jhydrol.2009.08.003
#   Kling, H., Fuchs, M. and Paulin, M. (2012). Runoff conditions in the upper Danube
#     basin under an ensemble of climate change scenarios. Journal of Hydrology
#     424-425, 264-277. doi:10.1016/j.jhydrol.2012.01.011
#   Willmott, C.J., Robeson, S.M. and Matsuura, K. (2012). A refined index of model
#     performance. International Journal of Climatology 32(13), 2088-2094.
#     doi:10.1002/joc.2419
#   Pushpalatha, R. et al. (2012). A review of efficiency criteria suitable for
#     evaluating low-flow simulations. Journal of Hydrology 420-421, 171-182.
#     doi:10.1016/j.jhydrol.2011.11.055
#   Zambrano-Bigiarini, M. and Bellin, A. (2012). Comparing goodness-of-fit measures
#     for calibration of models focused on extreme events. EGU General Assembly 2012,
#     Vienna, Austria, EGU2012-11549-1.
#   Garcia, F., Folton, N. and Oudin, L. (2017). Which objective function to
#     calibrate rainfall-runoff models for low-flow index simulations? Hydrological
#     Sciences Journal 62(7), 1149-1166. doi:10.1080/02626667.2017.1308511
#   Pool, S., Vis, M. and Seibert, J. (2018). Evaluating model performance: towards a
#     non-parametric variant of the Kling-Gupta efficiency. Hydrological Sciences
#     Journal 63(13-14), 1941-1953. doi:10.1080/02626667.2018.1552002
#   Liu, D. (2020). A rational performance criterion for hydrological model. Journal
#     of Hydrology 590, 125488. doi:10.1016/j.jhydrol.2020.125488
#   Lee, S. and Choi, Y.S. (2022). Comparison of objective functions for calibrating
#     conceptual rainfall-runoff models. Water 14(23), 3792. doi:10.3390/w14233792
#   Pizarro, A. and Jorquera, J. (2024). Advancing objective functions in
#     hydrological modelling: Integrating knowable moments for improved simulation
#     accuracy. Journal of Hydrology 634, 131071. doi:10.1016/j.jhydrol.2024.131071

box::use(
  stats[cor, quantile, sd]
)

# ---------------------------------------------------------------------------- #
# Internal helpers
# ---------------------------------------------------------------------------- #

# Drop the pairs where either series is missing, so that every metric below is
# computed on the same index set.
paired <- function(sim, obs, na.rm = TRUE) {
  sim <- as.numeric(sim)
  obs <- as.numeric(obs)
  if (length(sim) != length(obs)) {
    stop("'sim' and 'obs' must have the same length")
  }
  if (na.rm) {
    keep <- !is.na(sim) & !is.na(obs)
    sim <- sim[keep]
    obs <- obs[keep]
  }
  list(sim = sim, obs = obs)
}

# Constant added to both series before an inverse or log transform, so that zero
# flows do not blow up. See Pushpalatha et al. (2012).
epsilon_of <- function(obs, type = "none", value = NA) {
  switch(match.arg(type, c("none", "Pushpalatha2012", "otherFactor", "otherValue")),
    none = 0,
    Pushpalatha2012 = mean(obs) / 100,
    otherFactor = value * mean(obs),
    otherValue = value
  )
}

# Departures of the correlation, variability and bias terms from their ideal
# values. The 2009 and 2012 bias terms are ratios with an ideal of 1, whereas the
# 2021 bias term is a standardised mean difference with an ideal of 0, so each
# variant reports how far it sits from its own target rather than a raw value.
kge_deviations <- function(sim, obs, method = "2009", dispersion = sd) {
  mu_s <- mean(sim)
  mu_o <- mean(obs)
  sd_s <- dispersion(sim)
  sd_o <- dispersion(obs)
  vr <- if (method == "2012") (sd_s / mu_s) / (sd_o / mu_o) else sd_s / sd_o
  br <- if (method == "2021") (mu_s - mu_o) / sd_o else mu_s / mu_o - 1
  c(cor(sim, obs, method = "pearson") - 1, vr - 1, br)
}

# Euclidean distance from the ideal point, common to every KGE variant.
kge_score <- function(deviations, s) {
  1 - sqrt(sum((s * deviations)^2))
}

# ---------------------------------------------------------------------------- #
# Error magnitude
# ---------------------------------------------------------------------------- #

#' Mean error, in the units of the series.
#' @export
me <- function(sim, obs, na.rm = TRUE) {
  p <- paired(sim, obs, na.rm)
  mean(p$sim - p$obs)
}

#' Mean absolute error, in the units of the series.
#' @export
mae <- function(sim, obs, na.rm = TRUE) {
  p <- paired(sim, obs, na.rm)
  mean(abs(p$sim - p$obs))
}

#' Mean squared error.
#' @export
mse <- function(sim, obs, na.rm = TRUE) {
  p <- paired(sim, obs, na.rm)
  mean((p$sim - p$obs)^2)
}

#' Root mean squared error, in the units of the series.
#' @export
rmse <- function(sim, obs, na.rm = TRUE) {
  sqrt(mse(sim, obs, na.rm))
}

#' Root mean squared error with the mean bias removed.
#' @export
ubRMSE <- function(sim, obs, na.rm = TRUE) {
  sqrt(rmse(sim, obs, na.rm)^2 - me(sim, obs, na.rm)^2)
}

#' Root mean squared error as a percentage of the spread of the observations.
#' Rounded to `dec` decimals to match the convention this metric is reported with.
#' @export
nrmse <- function(sim, obs, na.rm = TRUE, norm = c("sd", "maxmin"), dec = 1) {
  norm <- match.arg(norm)
  p <- paired(sim, obs, na.rm)
  nval <- if (norm == "sd") sd(p$obs) else diff(range(p$obs))
  round(100 * rmse(p$sim, p$obs, na.rm = FALSE) / nval, dec)
}

#' Percent bias. Positive values mean the simulation runs high. Moriasi et al. (2007).
#' @export
pbias <- function(sim, obs, na.rm = TRUE, dec = 1) {
  p <- paired(sim, obs, na.rm)
  round(100 * sum(p$sim - p$obs) / sum(p$obs), dec)
}

#' Ratio of RMSE to the standard deviation of the observations. Moriasi et al. (2007).
#' @export
rsr <- function(sim, obs, na.rm = TRUE) {
  p <- paired(sim, obs, na.rm)
  rmse(p$sim, p$obs, na.rm = FALSE) / sd(p$obs)
}

#' Ratio of the standard deviations of simulated and observed values.
#' @export
rSD <- function(sim, obs, na.rm = TRUE) {
  p <- paired(sim, obs, na.rm)
  sd(p$sim) / sd(p$obs)
}

# ---------------------------------------------------------------------------- #
# Nash-Sutcliffe family
# ---------------------------------------------------------------------------- #

#' Nash-Sutcliffe efficiency. Nash and Sutcliffe (1970).
#' @export
NSE <- function(sim, obs, na.rm = TRUE) {
  p <- paired(sim, obs, na.rm)
  1 - sum((p$obs - p$sim)^2) / sum((p$obs - mean(p$obs))^2)
}

#' Modified Nash-Sutcliffe efficiency. With `j = 1` the squared errors of NSE are
#' replaced by absolute errors, which curbs the influence of a few large events.
#' Krause et al. (2005), Legates and McCabe (1999).
#' @export
mNSE <- function(sim, obs, na.rm = TRUE, j = 1) {
  p <- paired(sim, obs, na.rm)
  1 - sum(abs(p$obs - p$sim)^j) / sum(abs(p$obs - mean(p$obs))^j)
}

#' Relative Nash-Sutcliffe efficiency. Errors are scaled by the observation, which
#' shifts weight onto low flows. Krause et al. (2005).
#' @export
rNSE <- function(sim, obs, na.rm = TRUE) {
  p <- paired(sim, obs, na.rm)
  1 - sum(((p$obs - p$sim) / p$obs)^2) / sum(((p$obs - mean(p$obs)) / mean(p$obs))^2)
}

#' Nash-Sutcliffe efficiency with each squared error weighted by the observation,
#' which shifts weight onto high flows.
#' @export
wNSE <- function(sim, obs, na.rm = TRUE) {
  p <- paired(sim, obs, na.rm)
  1 - sum(p$obs * (p$obs - p$sim)^2) / sum(p$obs * (p$obs - mean(p$obs))^2)
}

#' Weighted seasonal Nash-Sutcliffe efficiency. Zambrano-Bigiarini and Bellin (2012).
#'
#' Observations at or above the high-flow quantile get weight `lambda`, those at or
#' below the low-flow quantile get `1 - lambda`, and the weight ramps linearly in
#' between. `lambda > 0.5` emphasises high flows, `lambda < 0.5` low flows.
#' @export
wsNSE <- function(sim, obs, na.rm = TRUE, j = 2, lambda = 0.95,
                  lQ.thr = 0.6, hQ.thr = 0.1) {
  p <- paired(sim, obs, na.rm)
  lQ <- quantile(p$obs, 1 - lQ.thr, names = FALSE)
  hQ <- quantile(p$obs, 1 - hQ.thr, names = FALSE)
  w <- ifelse(
    p$obs >= hQ, lambda,
    ifelse(
      p$obs <= lQ, 1 - lambda,
      (1 - lambda) + (2 * lambda - 1) * (p$obs - lQ) / (hQ - lQ)
    )
  )
  1 - sum(abs(w * (p$obs - p$sim))^j) / sum(abs(w * (p$obs - mean(p$obs)))^j)
}

# ---------------------------------------------------------------------------- #
# Indices of agreement
# ---------------------------------------------------------------------------- #

#' Index of agreement. Willmott (1981).
#' @export
d <- function(sim, obs, na.rm = TRUE) {
  p <- paired(sim, obs, na.rm)
  mu_o <- mean(p$obs)
  1 - sum((p$obs - p$sim)^2) / sum((abs(p$sim - mu_o) + abs(p$obs - mu_o))^2)
}

#' Refined index of agreement. Willmott et al. (2012). Ranges from -1 to 1, unlike
#' the original index, which is bounded below by 0.
#' @export
dr <- function(sim, obs, na.rm = TRUE) {
  p <- paired(sim, obs, na.rm)
  a <- sum(abs(p$sim - p$obs))
  b <- 2 * sum(abs(p$obs - mean(p$obs)))
  if (a <= b) 1 - a / b else b / a - 1
}

#' Modified index of agreement, using absolute rather than squared differences when
#' `j = 1`. Willmott et al. (1985), Legates and McCabe (1999).
#' @export
md <- function(sim, obs, na.rm = TRUE, j = 1) {
  p <- paired(sim, obs, na.rm)
  mu_o <- mean(p$obs)
  1 - sum(abs(p$obs - p$sim)^j) / sum((abs(p$sim - mu_o) + abs(p$obs - mu_o))^j)
}

#' Relative index of agreement, which scales differences by the observation and so
#' gives low flows more say. Krause et al. (2005).
#' @export
rd <- function(sim, obs, na.rm = TRUE) {
  p <- paired(sim, obs, na.rm)
  mu_o <- mean(p$obs)
  1 - sum(((p$obs - p$sim) / p$obs)^2) /
    sum(((abs(p$sim - mu_o) + abs(p$obs - mu_o)) / mu_o)^2)
}

#' Coefficient of persistence, comparing the model against a forecast that simply
#' carries the previous observation forward. Kitanidis and Bras (1980).
#' @export
cp <- function(sim, obs, na.rm = TRUE) {
  p <- paired(sim, obs, na.rm)
  n <- length(p$obs)
  1 - sum((p$obs[2:n] - p$sim[2:n])^2) / sum((p$obs[2:n] - p$obs[1:(n - 1)])^2)
}

# ---------------------------------------------------------------------------- #
# Correlation and regression
# ---------------------------------------------------------------------------- #

#' Pearson product-moment correlation coefficient.
#' @export
rPearson <- function(sim, obs, na.rm = TRUE) {
  p <- paired(sim, obs, na.rm)
  cor(p$sim, p$obs, method = "pearson")
}

#' Spearman rank correlation coefficient.
#' @export
rSpearman <- function(sim, obs, na.rm = TRUE) {
  p <- paired(sim, obs, na.rm)
  cor(p$sim, p$obs, method = "spearman")
}

#' Coefficient of determination, the fraction of observed variance the simulation
#' accounts for. Identical in form to NSE for a simulated-versus-observed pair, so
#' the two agree by construction. Note this is not the square of `rPearson`.
#' @export
R2 <- function(sim, obs, na.rm = TRUE) {
  p <- paired(sim, obs, na.rm)
  1 - sum((p$obs - p$sim)^2) / sum((p$obs - mean(p$obs))^2)
}

#' Coefficient of determination weighted by the regression slope, which penalises a
#' model that tracks the shape of the hydrograph but sits systematically high or
#' low. `b` is the slope of a regression through the origin. Krause et al. (2005).
#' @export
br2 <- function(sim, obs, na.rm = TRUE) {
  p <- paired(sim, obs, na.rm)
  b <- sum(p$sim * p$obs) / sum(p$obs^2)
  r2 <- R2(p$sim, p$obs, na.rm = FALSE)
  if (abs(b) <= 1) abs(b) * r2 else r2 / abs(b)
}

#' Volumetric efficiency, the fraction of observed volume delivered at the right
#' time. Criss and Winston (2008).
#' @export
VE <- function(sim, obs, na.rm = TRUE) {
  p <- paired(sim, obs, na.rm)
  1 - sum(abs(p$sim - p$obs)) / sum(p$obs)
}

# ---------------------------------------------------------------------------- #
# Kling-Gupta family
# ---------------------------------------------------------------------------- #

#' Kling-Gupta efficiency. Gupta et al. (2009) for `method = "2009"`, which measures
#' variability with the ratio of standard deviations; Kling et al. (2012) for
#' `method = "2012"`, which uses the ratio of coefficients of variation and so is
#' insensitive to bias in the variability term.
#'
#' `s` scales the correlation, variability and bias terms in that order.
#' @export
KGE <- function(sim, obs, s = c(1, 1, 1), na.rm = TRUE,
                method = c("2009", "2012", "2021")) {
  method <- match.arg(method)
  p <- paired(sim, obs, na.rm)
  kge_score(kge_deviations(p$sim, p$obs, method), s)
}

#' Kling-Gupta efficiency for low flows: the mean of KGE on the flows and KGE on
#' their reciprocals, which is what makes it sensitive to the recession limb.
#' Garcia et al. (2017).
#'
#' `epsilon.type` sets the constant added before inverting, so that zero flows do
#' not produce infinities. The default follows Pushpalatha et al. (2012) and adds
#' one hundredth of the mean observed flow; `"none"` adds nothing and will return
#' `NaN` if the record contains zeros.
#' @export
KGElf <- function(sim, obs, s = c(1, 1, 1), na.rm = TRUE,
                  method = c("2009", "2012", "2021"),
                  epsilon.type = c("Pushpalatha2012", "otherFactor", "otherValue", "none"),
                  epsilon.value = NA) {
  method <- match.arg(method)
  epsilon.type <- match.arg(epsilon.type)
  p <- paired(sim, obs, na.rm)
  eps <- epsilon_of(p$obs, epsilon.type, epsilon.value)
  high <- kge_score(kge_deviations(p$sim, p$obs, method), s)
  low <- kge_score(kge_deviations(1 / (p$sim + eps), 1 / (p$obs + eps), method), s)
  (high + low) / 2
}

#' Non-parametric Kling-Gupta efficiency. Correlation is measured by Spearman rank
#' and variability by the distance between normalised flow duration curves, so no
#' distributional assumption is made. Pool et al. (2018).
#' @export
KGEnp <- function(sim, obs, s = c(1, 1, 1), na.rm = TRUE) {
  p <- paired(sim, obs, na.rm)
  n <- length(p$obs)
  fdc_sim <- sort(p$sim) / (n * mean(p$sim))
  fdc_obs <- sort(p$obs) / (n * mean(p$obs))
  deviations <- c(
    cor(p$sim, p$obs, method = "spearman") - 1,
    -0.5 * sum(abs(fdc_sim - fdc_obs)),
    mean(p$sim) / mean(p$obs) - 1
  )
  kge_score(deviations, s)
}

#' Kling-Gupta efficiency with dispersion measured by the second knowable moment
#' rather than the standard deviation, which estimates more reliably from short or
#' heavy-tailed samples. Pizarro and Jorquera (2024).
#' @export
KGEkm <- function(sim, obs, s = c(1, 1, 1), na.rm = TRUE,
                  method = c("2009", "2012", "2021")) {
  method <- match.arg(method)
  p <- paired(sim, obs, na.rm)
  # Dispersion from the second knowable moment of the ascending-ordered sample.
  sd_km <- function(x) {
    n <- length(x)
    sqrt(2 * sum(2 * (seq_len(n) - 1) * sort(x)) / (n * (n - 1)))
  }
  kge_score(kge_deviations(p$sim, p$obs, method, dispersion = sd_km), s)
}

#' Liu mean efficiency, a two-term relative of KGE that folds correlation and
#' variability into a single term. Liu (2020).
#' @export
LME <- function(sim, obs, na.rm = TRUE) {
  p <- paired(sim, obs, na.rm)
  r <- cor(p$sim, p$obs, method = "pearson")
  alpha <- sd(p$sim) / sd(p$obs)
  beta <- mean(p$sim) / mean(p$obs)
  1 - sqrt((r * alpha - 1)^2 + (beta - 1)^2)
}

#' Lee and Choi efficiency. Splits the correlation-variability interaction into
#' `r * alpha` and `r / alpha`, which exposes errors in timing and spread that
#' cancel out in KGE. Lee and Choi (2022).
#' @export
LCE <- function(sim, obs, na.rm = TRUE) {
  p <- paired(sim, obs, na.rm)
  r <- cor(p$sim, p$obs, method = "pearson")
  alpha <- sd(p$sim) / sd(p$obs)
  beta <- mean(p$sim) / mean(p$obs)
  1 - sqrt((r * alpha - 1)^2 + (r / alpha - 1)^2 + (beta - 1)^2)
}

# ---------------------------------------------------------------------------- #
# Summary table
# ---------------------------------------------------------------------------- #

#' Goodness-of-fit summary, returned as a one-column matrix so it prints as a table
#' in the run summary. Argument defaults match the values these metrics are
#' conventionally reported with.
#' @export
gof <- function(sim, obs, na.rm = TRUE, j = 1, lambda = 0.95, norm = "sd",
                s = c(1, 1, 1), method = c("2009", "2012", "2021"),
                lQ.thr = 0.6, hQ.thr = 0.1, digits = 2) {
  method <- match.arg(method)
  p <- paired(sim, obs, na.rm)
  sim <- p$sim
  obs <- p$obs

  values <- c(
    ME = me(sim, obs, FALSE),
    MAE = mae(sim, obs, FALSE),
    MSE = mse(sim, obs, FALSE),
    RMSE = rmse(sim, obs, FALSE),
    ubRMSE = ubRMSE(sim, obs, FALSE),
    "NRMSE %" = nrmse(sim, obs, FALSE, norm = norm),
    "PBIAS %" = pbias(sim, obs, FALSE),
    RSR = rsr(sim, obs, FALSE),
    rSD = rSD(sim, obs, FALSE),
    NSE = NSE(sim, obs, FALSE),
    mNSE = mNSE(sim, obs, FALSE, j = j),
    rNSE = rNSE(sim, obs, FALSE),
    wNSE = wNSE(sim, obs, FALSE),
    wsNSE = wsNSE(sim, obs, FALSE, j = j, lambda = lambda, lQ.thr = lQ.thr, hQ.thr = hQ.thr),
    d = d(sim, obs, FALSE),
    dr = dr(sim, obs, FALSE),
    md = md(sim, obs, FALSE, j = j),
    rd = rd(sim, obs, FALSE),
    cp = cp(sim, obs, FALSE),
    r = rPearson(sim, obs, FALSE),
    R2 = R2(sim, obs, FALSE),
    bR2 = br2(sim, obs, FALSE),
    VE = VE(sim, obs, FALSE),
    KGE = KGE(sim, obs, s = s, na.rm = FALSE, method = method),
    # The summary table reports KGElf without an epsilon, so a record containing
    # zero flows yields NaN here rather than a value computed on shifted flows.
    KGElf = KGElf(sim, obs, s = s, na.rm = FALSE, method = method, epsilon.type = "none"),
    KGEnp = KGEnp(sim, obs, s = s, na.rm = FALSE),
    KGEkm = KGEkm(sim, obs, s = s, na.rm = FALSE, method = method),
    LME = LME(sim, obs, FALSE),
    LCE = LCE(sim, obs, FALSE)
  )

  matrix(round(values, digits), ncol = 1, dimnames = list(names(values), NULL))
}
