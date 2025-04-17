# Written by Cameron Bracken and Geoffrey Walters (2025)
# Please see the LICENSE file for license information

update_params <- function(p, p_names, default_pars) {
  # nl = nml::read_nml('namelist.HHWM8')
  # n_hrus = nl$INIT_CONTROL$n_hrus

  # SAC
  # hru_area - area (basin square miles)
  # uztwm - Upper zone tension water capacity, Units of MM
  # uzfwm - Upper zone free water capacity, Units of MM
  # lztwm - Lower zone tension water capacity, Units of MM
  # lzfpm - Lower zone primary free water capacity; Units of MM;
  #         LZFSM and LZFPM are input as total values and not just as the
  #         visible (channel component) portion.
  # lzfsm - Lower zone supplemental free water capacity; Units of MM
  # adimp - Additional impervious area; Units of percent/100
  # uzk -   Fractional daily upper zone free water withdrawal rate
  # lzpk -  Fractional daily primary withdrawal rate
  # lzsk -  Fractional daily supplemental withdrawal rate
  # zperc - Maximum percolation rate
  # rexp -  Exponent for the percolation equation
  # pctim - Minimum impervious area; Unit of percent/100
  # pfree - Percent/100 of percolated water which always goes directly
  #         to lower zone free water storages
  # riva -  Riparian vegetation area; Units of percent/100
  # side -  Ratio of non-channel baseflow (deep recharge) to
  #         channel (visible) baseflow
  # rserv - Percent/100 of lower zone free water which cannot be
  #         transferred to lower zone tension water

  # SNOW17
  # latitude - in degrees N
  # elev -   elevation (average, I am guessing)
  # scf -    Snowfall correction factor
  # mfmax -  Maximum non-rain melt factor; Units of MM/DEGC/6HR
  # mfmin -  Minimum non-rain melt factor; Units of MM/DEGC/6HR
  # uadj -   Average value of the wind function during rain-on-snow
  #          events Units of MM/MB
  # si -     Areal water-equivalent above which there is always 100%
  # pxtemp - Temperature that separates rain from snow; Units of DEGC
  # nmf -    Maximum negative melt factor; Units of MM/DEGC/6HR
  # tipm -   Antecedent snow temperature index parameter; Range is 0.1 to 1.0
  # mbase -  Base temperature for non-rain melt factor; Units of DEGC
  # plwhc -  Maximum amount of liquid water held against gravity
  #          drainage; decimal fraction
  # daygm -  Daily melt at the snow-soil interface Units of MM
  # adc1-11 - Areal snow cover at WE/Ai ratios of 0.0, 0.1, 0.2, 0.3, 0.4,
  #           0.5, 0.6, 0.7, 0.8, 0.9 and 1.0; Total 11 numbers.
  #           Decimal fraction

  # UH
  # unit_shape - gamma shape of the unit hydrograph
  # unit_scale - horizontal scale of the unit hydrograph
  #

  # zones = default_pars$zone |> unique |> as.character
  updated_pars <- copy(default_pars)
  # updated_pars$value = NULL
  updated_pars[, p_name := paste0(name, "_", zone)]
  updated_pars2 <- merge.data.table(updated_pars,
    data.table(updated = p, p_name = p_names),
    by = "p_name", all.x = TRUE
  )
  updated_pars2[!is.na(updated), value := updated]
  updated_pars2$updated <- NULL
  updated_pars2[order(zone, name)]
}

update_cu_params <- function(pars, zone1_name, cu_zones) {
  zone1_pars <- pars[zone == zone1_name]

  for (cu_zone_name in cu_zones) {
    # copy zone 1 pars, save the area and
    # dont copy the monthly peadj values, those are specially modified
    cu_zone_pars <- copy(zone1_pars[!(substr(name, 1, 6) == "peadj_")])
    cu_zone_area <- pars[zone == cu_zone_name & name == "zone_area"]$value
    cu_zone_peadj <- pars[substr(name, 1, 6) == "peadj_" & zone == cu_zone_name & type == "sac"]

    # update the zone and parameter names
    cu_zone_pars[, zone := cu_zone_name]
    cu_zone_pars[, p_name := paste0(name, "_", zone)]

    # stick all the params back together
    pars <- rbind(pars[zone != cu_zone_name | type == "consuse"], cu_zone_pars, cu_zone_peadj)
    # set back the correct cu zone area, for the uh calc
    pars[zone == cu_zone_name & name == "zone_area", value := cu_zone_area]
  }
  pars[order(zone, name)]
}

inst_to_ave <- function(forcing, sim_flow_cfs, agg_to_daily = FALSE) {
  sim <- as.data.table(forcing[[1]][, c("year", "month", "day", "hour")])

  # compute 6 hour average from instantaneous values
  sim[, next_sim := dplyr::lead(sim_flow_cfs)]
  sim[, sim_flow_cfs := (sim_flow_cfs + next_sim) / 2]
  sim$next_sim <- NULL

  # duplicate last value so no NA's
  sim[, sim_flow_cfs := nafill(sim_flow_cfs, "locf")]

  if (agg_to_daily) {
    simd <- sim[, list(sim_flow_cfs = mean(sim_flow_cfs)), by = .(year, month, day)]
    return(simd)
  } else {
    return(sim)
  }
}

model_wrapper <- function(p, p_names, dt_hours, default_pars, obs_daily, obs_inst,
                          forcing_raw, upflow, obj_fun, n_zones, cu_zones,
                          return_flow = FALSE) {
  # browser()

  # update params with current iteration values
  pars <- update_params(p, p_names, default_pars)

  # write out a file specific to the process
  # pid = Sys.getpid()
  # nodename = tolower(Sys.info()['nodename'])
  # fn = sprintf('%s_%s.txt',nodename,pid)
  # cat(p,'\n',file = fn,append=TRUE)

  # special case for consuse zones
  # duplicate zone 1 params and replace the cu zone params with those
  cu <- ifelse(length(cu_zones) > 0, TRUE, FALSE)
  if (cu) {
    # duplicate zone1 params
    zone1_name <- names(forcing_raw)[1]

    pars <- update_cu_params(pars, zone1_name, cu_zones)
  }

  # adjusted forcings
  if (n_zones > 0) {
    forcing_adj <- fa_nwrfc(dt_hours, forcing_raw, pars)
  } else {
    forcing_adj <- forcing_raw
  }

  # Run the model

  if (is.null(upflow)) {
    # !!includes chanloss but not consuse!!
    sim <- sac_snow_uh(dt_hours, forcing_adj, pars)
  } else if (n_zones == 0) {
    sim <- lagk(dt_hours, upflow, pars)
    sim <- chanloss(sim, forcing_adj, dt_hours, pars)
  } else {
    # !!includes chanloss but not consuse!!
    sim <- sac_snow_uh_lagk(dt_hours, forcing_adj, upflow, pars)
  }

  # format instant sim as datatable
  sim_inst <- as.data.table(forcing_adj[[1]][, c("year", "month", "day", "hour")])
  sim_inst[, sim_flow_cfs := sim]

  # compute 6 hour average from instantaneous values
  sim_per_avg <- inst_to_ave(forcing_adj, sim)

  # If consuse model is being used, need to recalculate the simulation to
  # consider diversion/return flow
  # !!This code will need to be updated to handle basins with multiple CU zones!!
  if (cu) {
    ### Daily calc

    # to match CHPS PET has to be shifted back so it is included in previous day
    sim_per_avg$pet <- forcing_adj[[1]] |>
      mutate(pet_mm = lead(pet_mm)) |>
      as_tibble() |>
      fill(pet_mm) |>
      select(pet_mm)

    # consuse needs a different flow name and a pet timeseries
    sim_daily <- sim_per_avg[, list(flow = mean(sim_flow_cfs), pet = sum(pet)), by = .(year, month, day)]
    # run consuse, replace flow with adjusted flow from cu model
    cu_out <- consuse(sim_daily, pars)

    # If there are multiple CU zone then sum each zone into a single datatable
    if (length(cu_zones) > 1) {
      cu_out <- rbindlist(cu_out, .id = "zone") |>
        group_by(month, day, year) |>
        summarise(
          qadj = sum(qadj),
          qdiv = sum(qdiv),
          qrfout = sum(qrfout),
          .groups = "drop"
        ) |>
        as.data.table()
    }

    # Dont need flow or pet columns anymore, so remove it
    # assign consuse qadj output
    sim_daily$sim_flow_cfs <- cu_out$qadj
    sim_daily <- subset(sim_daily, select = -c(flow, pet))

    ### Inst calc
    # apply daily consuse to 6 hour instantaneous data
    sim_inst <- cu_out |>
      as_tibble() |>
      mutate(qnetdiv = qdiv - qrfout) |>
      select(year, month, day, qnetdiv) |>
      right_join(sim_inst, by = c("year", "month", "day")) |>
      # to match CHPS the 00:00 value needs to be from the previous day
      mutate(qnetdiv = lag(qnetdiv)) |>
      fill(qnetdiv, .direction = "up") |>
      mutate(sim_flow_cfs = sim_flow_cfs - qnetdiv) |>
      select(-qnetdiv) |>
      as.data.table()
  } else {
    sim_daily <- sim_per_avg[, list(sim_flow_cfs = mean(sim_flow_cfs)), by = .(year, month, day)]
  }

  # merge in the obs and cut off the first year
  if (!is.null(obs_inst)) {
    results_inst <- merge(sim_inst, as.data.table(obs_inst), by = c("year", "month", "day", "hour"))[-(1:365 * 4)]
  }
  results_daily <- merge(sim_daily, as.data.table(obs_daily), by = c("year", "month", "day"))[-(1:365)]

  # browser()

  obj <- get(paste0(obj_fun, "_obj"))(results_daily, results_inst)

  if (return_flow) results_daily else c(obj_fun = obj)
}

run_controller_edds <- function(lower, upper, basin, dt_hours, default_pars,
                                obs_daily, obs_inst, forcing, upflow = NULL,
                                obj_fun = "rmse", n_zones, cu_zones = character(0),
                                n_cores, lite = FALSE) {
  # browser()
  # ptm = proc.time()

  t_iter <- ifelse(lite, 2500, 5000)

  out <- ep_dds(
    fn = model_wrapper,
    p_bounds = data.frame(name = names(lower), min = unname(lower), max = unname(upper)),
    t_iter = t_iter,
    n_cores = n_cores,
    r = 0.2,
    p_names = names(lower),
    dt_hours = dt_hours,
    default_pars = default_pars,
    obs_daily = obs_daily,
    obs_inst = obs_inst,
    forcing_raw = forcing,
    upflow = upflow,
    obj_fun = obj_fun,
    n_zones = n_zones,
    cu_zones = cu_zones
  )


  # end = (proc.time() - ptm)['elapsed']
  # hr = floor(end/60/60)
  # min = floor((end-hr*60*60)/60)
  # sec = round(end-hr*60*60-min*60)
  # cat(bgCyan('Optimization took',hr,'hr',min,'min',sec,'sec'),'\n')
  out
}


# Evolving Dynamically Dimensioned Search (EDDS)
ep_dds <- function(fn, p_bounds, t_iter = 1000, n_cores = 4, r = 0.2, ...) {
  # Parallel Registration
  my_cluster <- makeCluster(n_cores, type = "FORK") #' PSOCK'
  # Seeding using L'Ecuyer-CMRG
  RNGkind("L'Ecuyer-CMRG")
  clusterSetRNGStream(cl = my_cluster, iseed = NULL)
  # Set a fresh random seed for the master process
  set.seed(sample.int(.Machine$integer.max, 1))
  (stream <- .Random.seed)

  # Step 1: Check inputs
  if (!is.function(fn)) {
    stop(paste0("'fn' must be a function"))
  }
  if (!is.data.frame(p_bounds)) {
    stop("parameters must be given as a data frame")
  }
  required <- c("name", "min", "max")
  if (is.null(names(p_bounds)) || (!all(required %in% colnames(p_bounds)))) {
    stop(paste0(
      "missing column names in data frame of parameters,",
      " expecting '", paste(required, collapse = "', '"), "'"
    ))
  }
  # p_bounds does not have a default column, so this check is not being used
  if (any(p_bounds$max < p_bounds$min) || any(p_bounds$max < p_bounds$default) ||
    any(p_bounds$min > p_bounds$default)) {
    stop("parameter ranges not reasonable")
  }

  # Step 2: Initial evaluation of objective function
  i <- 0
  # Randomly select initial parameters
  p_bounds$initial <- apply(
    data.frame(p_bounds$min, p_bounds$max), 1,
    function(x) runif(1, x[1], x[2])
  )
  p_best <- setNames(p_bounds$initial, p_bounds$name)
  f_best <- fn(p_best, ...)
  if (any(is.na(f_best))) {
    stop("'fn' returns NA for initial parameter set")
  }
  if (is.null(names(f_best)) || (any(names(f_best) == ""))) {
    stop("'fn' must return a named vector")
  }

  # Save initial state
  p_trace <- data.frame(matrix(p_best, nrow = 1))
  names(p_trace) <- names(p_best)
  f_trace <- data.frame(matrix(f_best, nrow = 1))
  names(f_trace) <- names(f_best)

  message("                                                                                ")
  message("===============================================================================================")
  message("                                     Running optimization                     ")
  message("===============================================================================================")
  message("* Daily Flow Metrics Displayed                                                   ")

  sim <- model_wrapper(p = p_best, ..., return_flow = TRUE)[-(1:365 * 4)]
  nse_best <- NSE(sim$sim_flow_cfs, sim$flow_cfs)
  pbias_best <- pbias(sim$sim_flow_cfs, sim$flow_cfs)
  R2_best <- rPearson(sim$sim_flow_cfs, sim$flow_cfs)^2
  kge_best <- KGE(sim$sim_flow_cfs, sim$flow_cfs)

  message(
    "iter:", format(i, width = nchar(t_iter) + 1, justify = "right"),
    "  Fbest:", format(round(f_best, ifelse(f_best <= -100, 2, 3)), width = 9, nsmall = 3, justify = "right"),
    "  Piter_Len:", format(as.integer(1), width = 6, nsmall = 2, justify = "right"),
    "  Pbias:", format(round(pbias_best, 2), width = 8, nsmall = 2, justify = "right"), "%",
    "  R2:", format(round(R2_best, 2), width = 6, nsmall = 2, justify = "right"),
    "  NSE:", format(round(nse_best, 2), width = 6, nsmall = 2, justify = "right"),
    "  KGE:", format(round(kge_best, 2), width = 6, nsmall = 2, justify = "right")
  )

  # Initialize a list to store seeds for each iteration
  worker_seeds_list <- list()

  # Iteration loop
  while (i < t_iter) {
    # Dynamically change how often the parallel process check in depending on how
    # far along the run is.  Total steps dictate progression
    if (t_iter == 5000) {
      if (i < round(2 / 5 * t_iter, -2)) {
        p_iter <- 1000
      } else if (i < round(7 / 10 * t_iter, -2)) {
        p_iter <- 500
      } else if (i < round(9 / 10 * t_iter, -2)) {
        p_iter <- 100
      } else {
        p_iter <- 10
      }
    } else {
      if (i < round(2 / 5 * t_iter, -2)) {
        p_iter <- 1000
      } else if (i < round(1 / 2 * t_iter, -2)) {
        p_iter <- 500
      } else if (i < round(4 / 5 * t_iter, -2)) {
        p_iter <- 100
      } else {
        p_iter <- 10
      }
    }

    ep_dds_results <- clusterCall(
      cl = my_cluster,
      dds,
      fn = model_wrapper,
      p_bounds = p_bounds[, c("name", "max", "min")],
      f_best = f_best,
      p_best = p_best,
      f_trace = f_trace,
      p_trace = p_trace,
      c_iter = i,
      t_iter = t_iter,
      p_iter = p_iter,
      r = 0.2,
      p_names = names(lower),
      dt_hours = dt_hours,
      default_pars = default_pars,
      obs_daily = obs_daily,
      obs_inst = obs_inst,
      forcing_raw = forcing_raw,
      upflow = upflow,
      obj_fun = obj_fun,
      n_zones = n_zones,
      cu_zones = cu_zones
    )


    f_best_values <- sapply(ep_dds_results, function(x) x$f_best)
    par_best <- which.min(f_best_values)
    f_best <- ep_dds_results[[par_best]]$f_best
    p_best <- ep_dds_results[[par_best]]$p_best
    p_trace <- ep_dds_results[[par_best]]$p_trace
    f_trace <- ep_dds_results[[par_best]]$f_trace

    worker_seeds <- lapply(ep_dds_results, function(x) x$seed)
    worker_seeds_list[[as.character(i)]] <- unlist(worker_seeds)

    # Update counter
    i <- i + p_iter
    # if (!exists("stream") || length(stream) < 2 || !is.numeric(stream)) {
    #   stop("Invalid RNG state detected in `stream` before calling nextRNGStream()")
    # }
    nextRNGStream(stream)

    # Add a message regarding iteration and best run if i is divisable by 100
    if (i %% 100 == 0) {
      # print statement in case there are concerns about the random number generation
      # print(unlist(ep_dds['f_best',]))

      sim <- model_wrapper(p = p_best, ..., return_flow = TRUE)[-(1:365 * 4)]

      nse_best <- NSE(sim$sim_flow_cfs, sim$flow_cfs)
      pbias_best <- pbias(sim$sim_flow_cfs, sim$flow_cfs)
      R2_best <- rPearson(sim$sim_flow_cfs, sim$flow_cfs)^2
      kge_best <- KGE(sim$sim_flow_cfs, sim$flow_cfs)

      message(
        "iter:", format(i, width = nchar(t_iter) + 1, justify = "right"),
        "  Fbest:", format(round(f_best, ifelse(f_best <= -100, 2, 3)), width = 9, nsmall = 3, justify = "right"),
        "  Piter_Len:", format(as.integer(p_iter), width = 6, nsmall = 2, justify = "right"),
        "  Pbias:", format(round(pbias_best, 2), width = 8, nsmall = 2, justify = "right"), "%",
        "  R2:", format(round(R2_best, 2), width = 6, nsmall = 2, justify = "right"),
        "  NSE:", format(round(nse_best, 2), width = 6, nsmall = 2, justify = "right"),
        "  KGE:", format(round(kge_best, 2), width = 6, nsmall = 2, justify = "right")
      )
    }
  }
  # Stop the registered parallel forks
  stopCluster(cl = my_cluster)
  return(list(
    f_best = f_best, p_best = p_best, f_trace = f_trace, p_trace = p_trace, seeds = worker_seeds_list
  ))
}

#' DDS code built upon script by the following author:
#' @author David Kneis \email{david.kneis@@tu-dresden.de}
#'   https://github.com/dkneis/mcu/blob/master/R/dds.r
#'
#' @references The original algorithm is described in
#'   Tolson, B. A. and Shoemaker, C. A. (2007): Dynamically
#'   dimensioned search algorithm for computationally efficient watershed model
#'   calibration, Water Resources Research, 43, W01413,
#'   doi:10.1029/2005WR004723.
#'
#' An inplementation of the dynamically dimensioned search (DDS) algorithm for
#' stochastic optimization of computationally expensive objective functions. The
#' algorithm gradually shifts from global to local search.

dds <- function(fn, p_bounds, f_best, p_best, f_trace, p_trace, c_iter,
                t_iter = 10000, p_iter = 100, r = 0.2, ...) {
  i <- 1
  # Capture seed for verification
  worker_seed <- .Random.seed

  # Iteration loop
  while (i < (p_iter + 1)) {
    # Step 1: Select parameters for perturbation
    inds_perturb <- which(runif(nrow(p_bounds)) <= (1 - log(i + c_iter) / log(t_iter)))
    if (length(inds_perturb) == 0) {
      inds_perturb <- sample.int(n = nrow(p_bounds), size = 1)
    }

    # Step 2: Perturbation with reflection
    p_new <- p_best
    p_new[inds_perturb] <- p_best[inds_perturb] +
      r * (p_bounds$max[inds_perturb] - p_bounds$min[inds_perturb]) *
        rnorm(n = length(inds_perturb), mean = 0, sd = 1)

    # Reflection at min
    inds_out <- inds_perturb[p_new[inds_perturb] < p_bounds$min[inds_perturb]]
    if (length(inds_out) > 0) {
      p_new[inds_out] <- p_bounds$min[inds_out] + (p_bounds$min[inds_out] - p_new[inds_out])
      inds_out <- inds_out[p_new[inds_out] > p_bounds$max[inds_out]]
      if (length(inds_out) > 0) {
        p_new[inds_out] <- p_bounds$min[inds_out]
      }
    }
    # Reflection at max
    inds_out <- inds_perturb[p_new[inds_perturb] > p_bounds$max[inds_perturb]]
    if (length(inds_out) > 0) {
      p_new[inds_out] <- p_bounds$max[inds_out] - (p_new[inds_out] - p_bounds$max[inds_out])
      inds_out <- inds_out[p_new[inds_out] < p_bounds$min[inds_out]]
      if (length(inds_out) > 0) {
        p_new[inds_out] <- p_bounds$max[inds_out]
      }
    }

    # Step 3: Evaluate objective function and possibly update best solution
    f_new <- fn(p_new, ...)
    if (!is.na(f_new[1])) {
      if (f_new[1] <= f_best[1]) {
        f_best <- f_new
        p_best <- p_new
      }
    }

    # Step 4: Update counter
    i <- i + 1
    # Step 5: Save current state
    p_trace <- rbind(p_trace, p_best)
    f_trace <- rbind(f_trace, f_best)
  }
  # list(f_best,p_best,f_trace,p_trace)
  return(list(
    f_best = f_best, p_best = p_best, f_trace = f_trace, p_trace = p_trace, seed = worker_seed
  ))
}
