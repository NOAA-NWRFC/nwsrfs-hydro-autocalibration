# Written by Cameron Bracken and Geoffrey Walters (2025)
# Please see the LICENSE file for license information

# Instructions for Creating Objective Functions in the Calibration Program
#
# Below are the core guidelines for defining objective functions used by the
# calibration (`run-controller.R`) program:
#
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# DEVIATION FROM THESE NAMING CONVENTION MAY RESULT IN UNEXPECTED RESULTS
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#
# objective function naming convention is as follows:
#           <daily_metric>_<inst_metric>_obj
#
# 1. Clear Naming:
#    - Name your objective function explicitly so that other calibrators can
#      easily understand the general construct of the metric being used for
#      daily or subdaily data.
#
#    - The name must always end with "_obj" (e.g., nse_obj or custom_metric_obj).
#
# 2. Function Inputs:
#    - Your function must accept the following inputs: results_daily, results_inst.
#    - These inputs represent daily and instantaneous (subdaily) results,
#      respectively.
#
# 3. Optimization Direction:
#    - The optimizer minimizes the objective function score.
#    - For metrics like NSE or KGE, reverse their sign (e.g., multiply by -100)
#      to ensure the optimizer works correctly.
#
# 4. Weighting Considerations:
#    - Experiment with adjusting the weights of individual metrics within your
#      objective function.
#    - Consider the relative importance of daily versus instantaneous (subdaily)
#      data when assigning weights.
#
# 5. Filtering Considerations:
#    - Experiment with filtering out data based on threshold/quantiles
#    - Filters can also be applied to focus on certain times of the year
#
# Available Metrics (sign has already been reversed for all of these functions)
# nse_fun = Nash-Sutcliffe Efficiency
# lognse_fun = Log Nash-Sutcliffe Efficiency
# kge_fun = Kling-Gupta Efficiency (Kling et al. 2012)
# kgelf_fun = Kling-Gupta Efficiency (Kling et al. 2012) - low flow focus
# npbias_fun = 1-abs(%bias)/100
# r2_fun =  correlation correlation coefficient (r^2)
#
# Examples of Objective Function Names:
#
# 1. NSE for daily data and KGE for instantaneous data
#
#    nse_kge_obj
#
# 2. NSE + logNSE for daily data and no instantaneous data
#    "NULL" specifies that instantaneous data is not used.
#
#    nselognse_NULL_obj
#
# 3. No daily data and KGE for instantaneous data
#    "NULL" specifies that daily data is not used.
#
#    NULL_kge_obj
#
# 4. NSE for daily data and KGE with a 95th percentile focus for instantaneous data
#    "th" specifies the use of percentiles.
#
#    nse_kge95th_obj
#
# 4. KGE for daily data and NPBIAS with a 1,000 cfs focus for instantaneous data
#    "q" specifies the use of a flow threshold.
#
#    kge_npbias1000q_obj
#
# 5. NPBias for daily data (May–July) and KGE with a 95th percentile focus for instantaneous data
#    Months are specified with leading zeros and followed by "m".
#    "th" specifies the use of percentiles.
#
#    npbias050607m_kge95th_obj
#
# 6. NSE for daily data (weighted 25%) and KGE for instantaneous data (weighted 75%)
#    "W" specifies weighting as a fraction. Ensure that Daily W + Inst W = 1.
#
#    nse.25W_kge.75W_obj
#
# 7. Combined daily metrics: 0.2 NSE + 0.8 logNSE (weighted 40% overall)
#    and KGE for instantaneous data (weighted 60%)
#    "w" specifies the weighting between daily metrics (as fractions).
#    "W" specifies the weighting between daily and instantaneous data (as fractions).
#    Ensure that Daily W + Inst W = 1 and Daily metric#1 w + Daily metric#2 w = 1.
#
#    nse.2wlognse.8w.4W_kge.6W_obj
#
# Note:
# If no "W" or "w" appears in the function name, it is assumed that no weighting is used.

###############################################################################
#################### custom objective function space ##########################
###############################################################################





###############################################################################
#################### standard/example objective functions #####################
###############################################################################

# Daily NSE + daily logNSE
nselognse_NULL_obj <- function(results_daily, results_inst) {
  return(nse_fun(results_daily) + lognse_fun(results_daily))
}

# Daily NSE + daily NPBias
nsenpbias_NULL_obj <- function(results_daily, results_inst) {
  return(nse_fun(results_daily) + npbias_fun(results_daily))
}

# Daily KGE
kge_NULL_obj <- function(results_daily, results_inst) {
  return(kge_fun(results_daily))
}

# Instantaneous NSE + instantaneous logNSE
NULL_nselognse_obj <- function(results_daily, results_inst) {
  return(nse_fun(results_inst) + lognse_fun(results_inst))
}

# Instantaneous KGE
NULL_kge_obj <- function(results_daily, results_inst) {
  return(kge_fun(results_inst))
}

# Daily logNSE + instantaneous NSE
lognse_nse_obj <- function(results_daily, results_inst) {
  return(lognse_fun(results_daily) + nse_fun(results_inst))
}

# Daily logNSE + instantaneous KGE
lognse_kge_obj <- function(results_daily, results_inst) {
  return(lognse_fun(results_daily) + kge_fun(results_inst))
}

# Daily logNSE + instantaneous R2
lognse_r2_obj <- function(results_daily, results_inst) {
  return(lognse_fun(results_daily) + r2_fun(results_inst))
}

# Daily KGE low flow + instantaneous kge
kgelf_kge_obj <- function(results_daily, results_inst) {
  return(kgelf_fun(results_daily) + kge_fun(results_inst))
}

# Daily logNSE + instantaneous NPBias for the top 5% flows
lognse_npbias95th_obj <- function(results_daily, results_inst) {
  daily_obj <- lognse_fun(results_daily)

  filter_th <- 0.95
  inst_obj <- npbias_fun(results_inst |>
    filter(quantile(flow_cfs, filter_th, na.rm = TRUE) < flow_cfs))

  return(daily_obj + inst_obj)
}

# Daily kge + instantaneous NPBias for the flows larger than 2,000cfs
kge_npbias2000q_obj <- function(results_daily, results_inst) {
  daily_obj <- kge_fun(results_daily)

  filter_th <- 0.95
  inst_obj <- npbias_fun(results_inst |>
    filter(2000 < flow_cfs))

  return(daily_obj + inst_obj)
}

# Daily NPBias for May–July + instantaneous KGE
npbias050607m_kge_obj <- function(results_daily, results_inst) {
  filter_mths <- c(5, 6, 7)
  daily_obj <- npbias_fun(results_daily |>
    filter(month %in% filter_mths))

  inst_obj <- kge_fun(results_inst)

  return(daily_obj + inst_obj)
}

# Daily 0.75 logNSE + 0.25 NSE + instantaneous NPBias for the top 1% flows
nse.25wlognse.75w_npbias99th_obj <- function(results_daily, results_inst) {
  daily_obj <- 0.25 * nse_fun(results_daily) +
    0.75 * lognse_fun(results_daily)

  filter_th <- 0.99
  inst_obj <- npbias_fun(results_inst |>
    filter(quantile(flow_cfs, filter_th, na.rm = TRUE) < flow_cfs))

  return(daily_obj + inst_obj)
}

# Daily 40% logNSE + instantaneous 60% NSE for November–March
lognse.4W_nse1112010203m.6W_obj <- function(results_daily, results_inst) {
  daily_obj <- 0.4 * lognse_fun(results_daily)

  filter_mths <- c(11, 12, 1, 2, 3)
  inst_obj <- 0.6 * nse_fun(results_inst |>
    filter(month %in% filter_mths))

  return(daily_obj + inst_obj)
}

###############################################################################
######## These are metric functions for use in the objective function #########
###############################################################################

nse_fun <- function(result) {
  with(result, NSE(sim_flow_cfs, flow_cfs)) * -100
}

lognse_fun <- function(result) {
  with(result, NSE(log(sim_flow_cfs + 1), log(flow_cfs + 1))) * -100
}

kge_fun <- function(result, s = c(1, 1, 1)) {
  with(result, KGE(sim_flow_cfs, flow_cfs, method = "2012"), s = s) * -100
}

kgelf_fun <- function(result, s = c(1, 1, 1)) {
  with(result, KGElf(sim_flow_cfs, flow_cfs, method = "2012"), s = s) * -100
}

kgenp_fun <- function(result) {
  with(result, KGEnp(sim_flow_cfs, flow_cfs)) * -100
}

npbias_fun <- function(result) {
  pbias_obj <- (1 - abs(with(result, pbias(sim_flow_cfs, flow_cfs))) / 100) * -100
}

r2_fun <- function(result) {
  with(result, rPearson(sim_flow_cfs, flow_cfs))^2 * -100
}
