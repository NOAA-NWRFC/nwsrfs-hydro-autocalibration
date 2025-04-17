#!/usr/bin/env Rscript

# Written by Cameron Bracken and Geoffrey Walters (2025)
# Please see the LICENSE file for license information

# install.packages(c('xfun','import','devtools'))
xfun::pkg_load2(
  "magrittr", "dplyr", "data.table", "dtplyr", "hydroGOF",
  "digest", "lubridate", "readr", "tibble", "ggthemes",
  "gtable", "crayon", "rmarkdown", "stringr", "tidyr", "argparser",
  "plotly"
)
xfun::pkg_attach2("ggplot2")
import::from(magrittr, "%<>%", "%T>%", "%$%")
import::from(
  dplyr, filter, select, summarise, group_by, ungroup, mutate, arrange, rename, bind_rows,
  inner_join, anti_join, semi_join, distinct, right_join, lag, lead, pull, left_join
)
select <- dplyr::select
import::from(
  data.table, as.data.table, data.table, fread, merge.data.table, copy,
  melt.data.table, rbindlist, dcast.data.table, fwrite, setkey, setnames, nafill, setDT
)
import::from(dtplyr, lazy_dt)
import::from(hydroGOF, gof, NSE, pbias)
import::from(digest, digest)
import::from(lubridate, ymd_hms, yday, ymd_hm, year, "year<-", years, hours)
import::from(readr, read_delim, write_csv, read_csv, cols)
import::from(tibble, tibble, as_tibble, rownames_to_column)
import::from(ggthemes, colorblind_pal)
import::from(crayon, "%+%", green, bold, bgGreen, blue, inverse, red, bgCyan)
import::from(gridExtra, grid.arrange)
import::from(grid, unit.pmax, grid.newpage, grid.draw)
import::from(rmarkdown, render)
import::from(tools, file_path_sans_ext)
import::from(rlang, set_names)
import::from(stringr, str_subset, str_detect, str_locate, str_replace, str_replace_all)
import::from(tidyr, pivot_wider, pivot_longer, fill)
import::from(argparser, arg_parser, add_argument, parse_args)
import::from(plotly, ggplotly)
import::from(dplyr, n)
import::from(vctrs, vec_fill_missing)
##############################################################################
# !!To ensure correct parameter/forcings/upstream flow get mapped correctly to
# intended zone/process. The default or optimal parameter file, forcing list, and
# upstream flow list MUST be alphabetized by zone or basin names prior to call wrapper.r
# or rfchydromodels functions
import::from(
  rfchydromodels, sac_snow_uh, sac_snow, sac_snow_states, lagk, sac_snow_uh_lagk, lagk,
  forcing_adjust_map_pet_ptps, forcing_adjust_mat, uh, pet_hs, uh2p_get_scale, uh2p_cfs_in,
  consuse, fa_nwrfc, fa_adj_nwrfc, chanloss, rsnwelev
)
source("wrappers.R")
source("obj_fun.R")

parser <- arg_parser("Auto-calibration postprocessor", hide.opts = TRUE)

# by default ArgumentParser will add an help option
parser <- add_argument(parser, "--dir", help = "Input directory path containing basin directories")
parser <- add_argument(parser, "--basins", default = NA_character_, help = "Basins to run", nargs = Inf)
parser <- add_argument(parser, "--run", default = NA_character_, help = "Output directories to postprocess", nargs = Inf)
args <- parse_args(parser)

# args = list(dir='runs/opt_run/1zone',reportdir=NA,basins='TLMO3',run='results_por_01')
# args = list(dir='runs/cluster20-testing/2zone',reportdir=NA,basins='ELWW1',run='results_por_03')
# args = list(dir='runs/opt_run/2zone',reportdir=NA,basins='AUBW1',run='results_cv_1_01')

# all basins are sub dirs under here
output_dir <- args$dir
basins <- args$basins
run_dirs <- args$run

# inside of the basin dir
results_dir_prefix <- "results"
# at the same level as the basin dir
reports_dir_prefix <- ""

# will reside inside of the results dir
plot_dir <- "plots"
plot_dir_valid <- "plots_valid"
report_template <- "report-scripts/report-single-basin.Rmd"

dt_hours <- 6
tz <- "UTC"

# Full simulation for a calibration of a partial period of full record
full_sim_for_partial_calb <- FALSE
partial_calb_por <- NULL # c(2002,2030)

# turn off colored output because it doesn't output to files well
options("cli.num_colors" = 1)

# check if the user specified the basins, if not run them all
basins <- if (any(is.na(basins))) {
  list.dirs(output_dir, full.names = FALSE, recursive = FALSE)
} else {
  basins
}

for (basin in basins) {
  basin_dir <- file.path(output_dir, basin)
  results_dirs <- list.files(basin_dir, paste0(results_dir_prefix, "*"))

  # if the user specified an output directory, then use that
  if (!is.na(run_dirs)) results_dirs <- run_dirs

  for (results_dir in results_dirs) {
    # are we doing cross validation? if so, find which fold
    rd_split <- strsplit(results_dir, "_")[[1]]
    cv <- ifelse(rd_split[2] == "cv", TRUE, FALSE)
    por <- ifelse(rd_split[2] == "por", TRUE, FALSE)
    if (cv) {
      cvfold <- as.integer(rd_split[3])
    }

    output_path <- file.path(basin_dir, results_dir)
    if (!file.exists(output_path)) {
      message("Results dir ", output_path, " does not exist, skipping")
      next
    }
    plot_path <- file.path(output_path, plot_dir)
    plot_path_valid <- file.path(output_path, plot_dir_valid)
    rc_output_path <- file.path(output_path, "run_settings.txt")
    dir.create(plot_path_valid, showWarnings = FALSE)

    first_underscore <- str_locate(results_dir, "_")[1, "start"]
    report_num <- substr(results_dir, first_underscore + 1, nchar(results_dir))

    # html report name
    html_report_name <- gsub("Rmd", "html", gsub("single-basin", paste0(basin, "_", report_num), report_template)) |>
      basename()

    if (!file.exists(rc_output_path)) next

    cat(inverse$bold("Postprocessing results:", output_path), "\n")

    #########################################
    # Configuration and autocalb output files
    #########################################

    # parameter: limits, default, and optimized values
    limits <- fread(file.path(basin_dir, "pars_limits.csv"))
    optimal_pars <- fread(file.path(output_path, "pars_optimal.csv"))
    optimal_pars <- optimal_pars[order(zone, name)]

    # Get zone configuration
    zones <- optimal_pars$zone |>
      unique() |>
      as.character() |>
      sort() |>
      str_subset("-")
    n_zones <- length(zones)
    if (n_zones > 0) {
      zone_lookup <- setNames(paste0("_", zones), paste0("_", as.character(1:n_zones)))
      zone_ids <- strsplit(zones, "-") |> sapply("[", 2)
    } else {
      zone_lookup <- zone_ids <- NULL
      zones <- paste0(basin, "-1")
    }

    # consuse zones
    cu_zones <- zones |>
      str_subset("-CU") |>
      sort()
    n_cu_zones <- length(cu_zones)
    cu <- ifelse(length(cu_zones) > 0, TRUE, FALSE)

    # forcing data
    forcing_raw <- list()
    for (zone in zones) {
      forcing_raw[[zone]] <- fread(file.path(basin_dir, paste0("forcing_por", "_", zone, ".csv")))
    }
    # Special process to get a full simulation for a autocalibration run on subset of the AORC POR
    if (full_sim_for_partial_calb & n_zones > 0) {
      forcing_raw_subset <- list()
      for (zone in zones) {
        forcing_raw_subset[[zone]] <- fread(file.path(basin_dir, paste0("forcing_por", "_", zone, ".csv"))) |>
          mutate(wyear = ifelse(month >= 10, year + 1, year)) |>
          filter(wyear >= partial_calb_por[1] & wyear <= partial_calb_por[2]) |>
          select(-wyear) |>
          as.data.table()
      }
      forcing_subset_climo <- fa_adj_nwrfc(dt_hours, forcing_raw_subset, optimal_pars, return_climo = TRUE) |>
        data.matrix()
      if (n_zones > 0) {
        forcing <- fa_nwrfc(dt_hours, forcing_raw, optimal_pars, climo = forcing_subset_climo)
      } else {
        forcing <- forcing_raw
      }
    } else {
      forcing_raw_subset <- forcing_subset_climo <- NULL
      if (n_zones > 0) {
        forcing <- fa_nwrfc(dt_hours, forcing_raw, optimal_pars)
      } else {
        forcing <- forcing_raw
      }
    }

    # If no zones, then the run is a routing reach only.  Revert dummy zone and zones back to NULL
    if (n_zones == 0) {
      zones <- NULL
    }

    # Daily streamflow data
    obs_daily <- fread(file.path(basin_dir, paste0("flow_daily_", basin, ".csv"))) |>
      as_tibble() |>
      mutate(imputed = ifelse(str_detect(Source, "Data_Ranger"), T, F)) |>
      select(-Source)

    # Instantaneous streamflow data (if exists)
    if (file.exists(file.path(basin_dir, paste0("flow_instantaneous_", basin, ".csv")))) {
      obs_inst <- fread(file.path(basin_dir, paste0("flow_instantaneous_", basin, ".csv")))[, -"Source"]

      # obs_daily_q95 = quantile(obs_daily$flow_cfs,.95,na.rm = TRUE,names = FALSE)

      # Add all timesteps to the instantaneous record
      obs_inst <- obs_inst |>
        # filter(flow_cfs>=obs_daily_q95) |>
        right_join(
          forcing_raw[[1]] |>
            as.data.table() |>
            dplyr::select(-c("map_mm", "mat_degc", "ptps")),
          by = c("year", "month", "day", "hour")
        ) |>
        as_tibble()
    } else {
      obs_inst <- obs_daily_q95 <- NULL
    }

    # flow data for upstream points to route down (if it exists)
    upflow_files <- list.files(basin_dir, "upflow_*", full.names = TRUE) |> sort()
    n_upstream <- length(upflow_files)
    upflow <- NULL
    if (n_upstream > 0) {
      upstream_lids <- gsub("upflow_", "", gsub(".csv", "", basename(upflow_files)))
      upflow <- list()
      for (u in 1:n_upstream) {
        upflow[[upstream_lids[u]]] <- fread(upflow_files[u]) |>
          right_join(
            forcing_raw[[1]] |>
              as.data.table() |>
              dplyr::select(-c("map_mm", "mat_degc", "ptps")),
            by = c("year", "month", "day", "hour")
          ) |>
          mutate(flow_cfs = vctrs::vec_fill_missing(flow_cfs, max_fill = 4)) |>
          as_tibble()
      }
      upflow_lookup <- setNames(paste0("_", upstream_lids), paste0("_", as.character(1:n_upstream)))
    }

    #########################################
    # Compute and write simulation timeseries
    ##########################################

    if (n_upstream == 0) {
      # no routing, includes chanloss
      optimal_sim_raw <- sac_snow_uh(dt_hours, forcing, optimal_pars)
    } else if (n_zones == 0) {
      # route only
      optimal_sim_raw <- lagk(dt_hours, upflow, optimal_pars)
      optimal_sim_raw <- chanloss(optimal_sim_raw, forcing, dt_hours, optimal_pars)
    } else {
      # local+routing, includes chanloss
      optimal_sim_raw <- sac_snow_uh_lagk(dt_hours, forcing, upflow, optimal_pars)
    }

    # Get optimal sim runs for different time steps
    optimal_sim_6hr_inst <- forcing[[1]][, c("year", "month", "day", "hour")]
    optimal_sim_6hr_inst$sim_flow_cfs <- optimal_sim_raw

    optimal_sim_6hr_ave <- inst_to_ave(forcing, optimal_sim_raw) #|>
    # merge(obs_daily,by=c('year','month','day'))
    optimal_sim_daily <- inst_to_ave(forcing, optimal_sim_raw, agg_to_daily = TRUE) #|>
    # merge(obs_daily,by=c('year','month','day'))

    # Get sac and snow states
    if (n_zones > 0) {
      # re-run the optimal state sim
      optimal_states_6hr <- sac_snow_states(dt_hours, forcing, optimal_pars) |>
        as.data.table() |>
        set_names(~ str_replace_all(.x, zone_lookup))
      tci <- sac_snow(6, forcing, optimal_pars)
      if (n_zones > 1) {
        zone_flow_cfs_raw <- uh(6, tci, optimal_pars, sum_zones = FALSE)
      } else {
        zone_flow_cfs_raw <- uh(6, tci, optimal_pars) |>
          as.data.frame()
      }
      colnames(zone_flow_cfs_raw) <- paste0("flow_", zones)

      optimal_states_6hr <- cbind(optimal_states_6hr, zone_flow_cfs_raw)
    } else {
      optimal_states_6hr <- tci <- zone_flow_cfs_raw <- NULL
    }

    # Get Routing States
    if (n_upstream > 0) {
      optimal_routed_raw <- lagk(dt_hours, upflow, optimal_pars, sum_routes = FALSE) |> cbind()
      colnames(optimal_routed_raw) <- upstream_lids
      optimal_routed <- cbind(
        data.table(
          year = upflow[[1]]$year,
          month = upflow[[1]]$month,
          day = upflow[[1]]$day,
          hour = upflow[[1]]$hour
        ),
        data.table(optimal_routed_raw)
      )
      optimal_routed_states <- lagk(dt_hours, upflow, optimal_pars, return_states = TRUE) |>
        as.data.table() |>
        set_names(~ str_replace_all(.x, upflow_lookup))
    } else {
      optimal_routed_raw <- optimal_routed <- optimal_routed_states <- NULL
    }

    # stop()
    if (cu) {
      # To match CHPS PET has to be shifted back so it is included in previous day
      optimal_sim_6hr_ave$pet <- forcing[[1]] |>
        mutate(pet_mm = lead(pet_mm)) |>
        as_tibble() |>
        fill(pet_mm) |>
        select(pet_mm)

      # consuse needs a different flow name and a pet timeseries
      optimal_sim_daily <- optimal_sim_6hr_ave[, list(flow = mean(sim_flow_cfs), pet = sum(pet)), by = .(year, month, day)]
      # run consuse, replace flow with adjusted flow from cu model
      optimal_cu_states <- consuse(optimal_sim_daily, optimal_pars)

      # If there are multiple CU zone then sum each zone into a single datatable
      if (length(cu_zones) > 1) {
        optimal_cu_daily <- rbindlist(optimal_cu_states, idcol = "zone") |>
          group_by(month, day, year) |>
          summarise(
            qadj = sum(qadj),
            qdiv = sum(qdiv),
            qrfout = sum(qrfout),
            .groups = "drop"
          ) |>
          as.data.table()
      } else {
        optimal_cu_daily <- optimal_cu_states
      }

      # Dont need flow or pet columns anymore, so remove it
      # assign consuse qadj output
      optimal_sim_daily$sim_flow_cfs <- optimal_cu_daily$qadj
      optimal_sim_daily <- subset(optimal_sim_daily, select = -c(flow, pet))

      # apply daily consuse to 6 hour instantaneous data
      optimal_sim_6hr_inst <- optimal_cu_daily |>
        as_tibble() |>
        mutate(qnetdiv = qdiv - qrfout) |>
        select(year, month, day, qnetdiv) |>
        right_join(optimal_sim_6hr_inst, by = c("year", "month", "day")) |>
        # to match CHPS the 00:00 value needs to be from the previous day
        mutate(qnetdiv = lag(qnetdiv)) |>
        fill(qnetdiv, .direction = "up") |>
        mutate(sim_flow_cfs = sim_flow_cfs - qnetdiv) |>
        select(-qnetdiv) |>
        as.data.table()

      # re compute the 6 hour ave based with the consuse baked in
      optimal_sim_6hr_ave <- inst_to_ave(forcing, optimal_sim_6hr_inst$sim_flow_cfs)

      # Format the consuse states varible
      if (length(cu_zones) == 1) {
        setnames(optimal_cu_states,
          old = names(optimal_cu_states)[-1:-3],
          paste0(names(optimal_cu_states)[-1:-3], "_", cu_zones)
        )
      } else {
        optimal_cu_states <- rbindlist(optimal_cu_states, idcol = "zone") |>
          as_tibble() |>
          pivot_wider(
            names_from = zone,
            values_from = c(qnat, qadj, qdiv, qrfin, qrfout, qol, qcd, ce, rfstor)
          ) |>
          as.data.table()
      }
    } else {
      optimal_cu_daily <- optimal_cu_states <- NULL
    }

    # merge flow observations
    optimal_sim_daily <- optimal_sim_daily |>
      merge(obs_daily, by = c("year", "month", "day"), all = TRUE)

    if (!is.null(obs_inst)) {
      optimal_sim_6hr_inst <- optimal_sim_6hr_inst |>
        merge(obs_inst, by = c("year", "month", "day", "hour"), all = TRUE)
    } else {
      optimal_sim_6hr_inst$flow_cfs <- NA
    }

    if (cv) {
      forcing_valid_fold <- fread(file.path(basin_dir, paste0("forcing_validation_cv_", cvfold, "_", zone, ".csv")))
      forcing_valid_fold_daily <- forcing_valid_fold |>
        group_by(year, month, day) |>
        summarise(hour = mean(hour), .groups = "drop") |>
        select(-hour) |>
        as.data.table()

      # include only the validation period
      valid_sim_daily <- semi_join(as_tibble(optimal_sim_daily), as_tibble(forcing_valid_fold_daily),
        by = c("year", "month", "day")
      ) |> as.data.table()
      # valid_sim_6hr_ave = semi_join(as_tibble(optimal_sim_6hr_ave), as_tibble(forcing_valid_fold),
      #                               by=c("year", "month", "day", "hour")) |> as.data.table()
      valid_sim_6hr_inst <- semi_join(as_tibble(optimal_sim_6hr_inst), as_tibble(forcing_valid_fold),
        by = c("year", "month", "day", "hour")
      ) |> as.data.table()

      # write the validation period simulations
      write_csv(valid_sim_daily, file.path(output_path, "validation_daily.csv"))
      write_csv(valid_sim_6hr_inst, file.path(output_path, "validation_6hr_inst.csv"))

      # cut out the validation period
      optimal_sim_daily <- anti_join(optimal_sim_daily, forcing_valid_fold_daily,
        by = c("year", "month", "day")
      ) |> as.data.table()
      # optimal_sim_6hr_ave = anti_join(optimal_sim_6hr_ave, forcing_valid_fold,
      #                                 by=c("year", "month", "day", "hour")) |> as.data.table()
      optimal_sim_6hr_inst <- anti_join(optimal_sim_6hr_inst, forcing_valid_fold,
        by = c("year", "month", "day", "hour")
      ) |> as.data.table()

      # Create sac and snow calb/valid timeseries and write valid period simulation
      if (n_zones > 0) {
        valid_states_6hr <- semi_join(as_tibble(optimal_states_6hr), as_tibble(forcing_valid_fold),
          by = c("year", "month", "day", "hour")
        ) |> as.data.table()

        write_csv(valid_states_6hr, file.path(output_path, "validation_states_6hr.csv"))

        optimal_states_6hr <- anti_join(optimal_states_6hr, forcing_valid_fold,
          by = c("year", "month", "day", "hour")
        ) |> as.data.table()
      } else {
        valid_states_6hr <- optimal_states_6h <- NULL
      }

      # Create Lagk calb/valid timeseries and write valid period simulation
      if (n_upstream != 0) {
        valid_routed <- semi_join(as_tibble(optimal_routed), as_tibble(forcing_valid_fold),
          by = c("year", "month", "day", "hour")
        ) |> as.data.table()
        valid_lagk_states <- semi_join(as_tibble(optimal_routed_states), as_tibble(forcing_valid_fold),
          by = c("year", "month", "day", "hour")
        ) |> as.data.table()

        write_csv(valid_routed, file.path(output_path, "validation_routed_flow.csv"))
        write_csv(valid_lagk_states, file.path(output_path, "validation_routed_states.csv"))

        optimal_routed <- anti_join(optimal_routed, forcing_valid_fold,
          by = c("year", "month", "day", "hour")
        ) |> as.data.table()
        optimal_routed_states <- anti_join(optimal_routed_states, forcing_valid_fold,
          by = c("year", "month", "day", "hour")
        ) |> as.data.table()
      } else {
        valid_routed <- valid_lagk_states <- NULL
      }

      # Create consuse calb/valid timeseries and write valid period simulation
      if (cu) {
        valid_cu_states <- semi_join(as_tibble(optimal_cu_states), as_tibble(forcing_valid_fold_daily),
          by = c("year", "month", "day")
        ) |> as.data.table()

        write_csv(valid_cu_states, file.path(output_path, "validation_consuse_states.csv"))

        optimal_cu_states <- anti_join(optimal_cu_states, forcing_valid_fold_daily,
          by = c("year", "month", "day")
        ) |> as.data.table()
      } else {
        valid_cu_daily <- NULL
      }
    }

    write_csv(optimal_sim_daily, file.path(output_path, "optimal_daily.csv"))
    write_csv(optimal_sim_6hr_inst, file.path(output_path, "optimal_6hr_inst.csv"))
    if (n_zones > 0) {
      write_csv(optimal_states_6hr, file.path(output_path, "optimal_states_6hr.csv"))
    }
    if (n_upstream != 0) {
      write_csv(optimal_routed, file.path(output_path, "optimal_routed_flow.csv"))
      write_csv(optimal_routed_states, file.path(output_path, "optimal_routed_states.csv"))
    }
    if (cu) {
      write_csv(optimal_cu_states, file.path(output_path, "optimal_consuse_states.csv"))
    }

    ###################
    # ADC table csv
    ###################
    adc_y <- seq(0, 1, by = 0.1)
    adc_step_crawl <- .0001
    adc_sf <- 1 / adc_step_crawl

    for (z in zones) {
      # skip if zone is a glacier
      # if(str_detect(z,'-G')) next

      adc_x <- rep(0, 11)
      adc_pars <- optimal_pars[name %in% paste0("adc_", letters[1:3]) & zone == z]$value
      adc_x_crawl <- 1 + adc_step_crawl
      adc_y_crawl <- 1 + adc_step_crawl
      for (j in length(adc_y):1) {
        while (as.integer(adc_y_crawl * adc_sf) > as.integer(adc_y[j] * adc_sf) &
          as.integer(adc_x_crawl * adc_sf) > as.integer(.05 * adc_sf)) {
          adc_x_crawl <- adc_x_crawl - adc_step_crawl
          adc_y_crawl <- adc_pars[1] * adc_x_crawl^adc_pars[2] + (1 - adc_pars[1]) * adc_x_crawl^adc_pars[3]
        }
        adc_x[j] <- adc_x_crawl
      }
      adc <- data.frame(`We_Ai` = adc_y, As = adc_x)
      fwrite(adc, file.path(output_path, sprintf("adc_%s.csv", z)))
    }

    ################
    # UH table csv
    #################
    for (z in zones) {
      # convert sqkm to sqmi
      area_sqmi <- optimal_pars[zone == z & name == "zone_area"]$value * 0.386102
      shape <- optimal_pars[zone == z & name == "unit_shape"]$value
      toc_gis <- optimal_pars[zone == z & name == "unit_toc"]$value
      toc_adj <- optimal_pars[zone == z & name == "unit_toc_adj"]$value

      toc <- toc_gis * toc_adj
      scale <- uh2p_get_scale(shape, toc, 1)

      # generate the UH at 1 hour timestep, it was calibrated at 6 hour
      uh1hr <- uh2p_cfs_in(shape, scale, 1, area_sqmi)
      uh6hr <- uh2p_cfs_in(shape, scale, 6, area_sqmi)
      fwrite(data.frame(y = uh1hr), file.path(output_path, sprintf("uh_1hr_%s.csv", z)))
      fwrite(data.frame(y = uh6hr), file.path(output_path, sprintf("uh_6hr_%s.csv", z)))
    }

    ####################
    # Cold states csv
    #####################
    cold_state_vars <- c("swe", "uztwc", "uzfwc", "lztwc", "lzfsc", "lzfpc", "adimc")
    for (z in zones) {
      # zone_suffix = paste0('_',as.character(which(zones==z)))
      zone_suffix <- paste0("_", z)
      cold_states_zonei <- paste0(cold_state_vars, zone_suffix)
      cold_states <- optimal_states_6hr[365 * 4 + 1, cold_states_zonei, with = F] |>
        melt.data.table(measure.vars = cold_states_zonei)
      cold_states[, variable := gsub(zone_suffix, "", variable)]
      fwrite(cold_states, file.path(output_path, sprintf("cold_states_%s.csv", z)))
    }

    ##########################
    # Forcing factors csv
    ###########################

    if (n_zones > 0) {
      if (full_sim_for_partial_calb) {
        factors <- fa_adj_nwrfc(dt_hours, forcing_raw_subset, optimal_pars, climo = forcing_subset_climo)
      } else {
        factors <- fa_adj_nwrfc(dt_hours, forcing_raw, optimal_pars)
      }

      factors <- factors |>
        as.data.frame() |>
        rename(map_fac = map_adj, mat_fac = mat_adj, ptps_fac = ptps_adj, pet_fac = pet_adj) |>
        mutate(month = 1:12) |>
        select(month, map_fac, mat_fac, pet_fac, ptps_fac)

      write_csv(round(factors, 4), file.path(output_path, "fa_factors.csv"))
    }

    #####################
    # lag-k table csv
    #####################
    if (n_upstream > 0) {
      # output lag and k tables
      lagtbl <- ktbl <- lagktbl <- list()
      for (u in 1:n_upstream) {
        # upflow[[u]] = fread(upflow_files[u])
        ndq <- 0
        lagvec <- kvec <- qvec <- numeric(11)

        lagtbl_a <- optimal_pars[zone == upstream_lids[u] & name == "lagtbl_a"]$value
        lagtbl_b <- optimal_pars[zone == upstream_lids[u] & name == "lagtbl_b"]$value
        lagtbl_c <- optimal_pars[zone == upstream_lids[u] & name == "lagtbl_c"]$value
        lagtbl_d <- optimal_pars[zone == upstream_lids[u] & name == "lagtbl_d"]$value
        ktbl_a <- optimal_pars[zone == upstream_lids[u] & name == "ktbl_a"]$value
        ktbl_b <- optimal_pars[zone == upstream_lids[u] & name == "ktbl_b"]$value
        ktbl_c <- optimal_pars[zone == upstream_lids[u] & name == "ktbl_c"]$value
        ktbl_d <- optimal_pars[zone == upstream_lids[u] & name == "ktbl_d"]$value
        lagk_qmax <- optimal_pars[zone == upstream_lids[u] & name == "lagk_qmax"]$value
        lagk_lagmax <- optimal_pars[zone == upstream_lids[u] & name == "lagk_lagmax"]$value
        lagk_kmax <- optimal_pars[zone == upstream_lids[u] & name == "lagk_kmax"]$value
        lagk_qmin <- optimal_pars[zone == upstream_lids[u] & name == "lagk_qmin"]$value
        lagk_lagmin <- optimal_pars[zone == upstream_lids[u] & name == "lagk_lagmin"]$value
        lagk_kmin <- optimal_pars[zone == upstream_lids[u] & name == "lagk_kmin"]$value

        for (i in 1:11) {
          qvec[i] <- ndq * (lagk_qmax - lagk_qmin) + lagk_qmin

          lag_entry <- lagtbl_a * (ndq - lagtbl_d)**2 + lagtbl_b * ndq + lagtbl_c
          k_entry <- ktbl_a * (ndq - ktbl_d)**2 + ktbl_b * ndq + ktbl_c

          if (lag_entry > 0 && lag_entry < 1) {
            lagvec[i] <- lag_entry * (lagk_lagmax - lagk_lagmin) + lagk_lagmin
          } else if (lag_entry >= 1) {
            lagvec[i] <- lagk_lagmax
          } else {
            lagvec[i] <- lagk_lagmin
          }

          if (k_entry > 0 && k_entry < 1) {
            kvec[i] <- k_entry * (lagk_kmax - lagk_kmin) + lagk_kmin
          } else if (k_entry >= 1) {
            kvec[i] <- lagk_kmax
          } else {
            kvec[i] <- lagk_kmin
          }

          ndq <- ndq + .1
        }
        lagtbl[[upstream_lids[u]]] <- data.table(Q = qvec, lag = lagvec)
        ktbl[[upstream_lids[u]]] <- data.table(Q = qvec, k = kvec)
        lagktbl[[u]] <- data.table(qvec, lagvec, kvec)
        names(lagktbl[[u]]) <- paste0(upstream_lids[u], "_", c("Q", "lag", "k"))
      }
      write_csv(do.call(cbind, lagktbl), file.path(output_path, "lagk_tables.csv"))
    }

    ####################
    # Performance Plots
    ####################

    output_plots <- function(x, plot_path) {
      # plot obs vs predicted scatter plot
      scatter <- ggplot(x) +
        geom_point(aes(flow_cfs, sim_flow_cfs)) +
        geom_abline(aes(slope = 1, intercept = 0)) +
        theme_minimal()
      ggsave(sprintf("%s/scatter.pdf", plot_path), scatter, width = 6, height = 6)

      # plot the hydrograph for each calendar year
      daily <- data.frame(x) |>
        mutate(
          wyear = ifelse(month >= 10, year + 1, year),
          datetime := ISOdatetime(year, month, day, 0, 0, 0)
        )

      # check for a "hanging chad" water year
      wyears_check <- daily |>
        group_by(wyear) |>
        summarise(count = n(), .groups = "drop") |>
        filter(count < 31) |>
        select(wyear)

      daily <- daily |> filter(!(wyear %in% wyears_check$wyear))

      wyears <- unique(daily$wyear)

      # check for a "hanging chad" water year


      # plot the hydrograph for each calendar year
      pdf(sprintf("%s/annual_hydrographs.pdf", plot_path), width = 7, height = 5)
      for (wy in wyears) {
        timeseries <- ggplot(daily[daily$wyear == wy, ]) +
          geom_line(aes(datetime, sim_flow_cfs), color = "steelblue") +
          geom_line(aes(datetime, flow_cfs)) +
          theme_minimal() +
          ggtitle(sprintf("%s (black: obs, blue: sim)", wy)) +
          xlab("") +
          ylab("Flow [cfs]") +
          scale_x_datetime(
            date_labels = "%b", date_breaks = "1 month",
            minor_breaks = NULL
          )

        print(timeseries)
        # old_fn = sprintf('%s/%s.pdf',plot_path,y)
        # if(file.exists(old_fn))unlink(old_fn)
        # ggsave(sprintf('%s/%s.pdf',plot_path,y),timeseries,width=7,height=5)
      }
      dev.off()

      cyclical_daily <- data.frame(x) |>
        group_by(month, day) |>
        summarise(
          sim_flow_cfs = mean(sim_flow_cfs),
          flow_cfs = mean(flow_cfs),
          .groups = "drop"
        ) |>
        mutate(
          wyear = ifelse(month >= 10, 1999, 2000),
          plot_date = ISOdatetime(wyear, month, day, 0, 0, 0)
        )
      p_cyclical <- cyclical_daily |> ggplot() +
        geom_line(aes(plot_date, flow_cfs)) +
        geom_line(aes(plot_date, sim_flow_cfs), color = "steelblue") +
        scale_x_datetime(
          date_labels = "%b", date_breaks = "1 month",
          minor_breaks = NULL
        ) +
        ggtitle("black: obs, blue: sim") +
        xlab("") +
        ylab("Flow [cfs]") +
        theme_minimal()
      p_cyclical_log <- cyclical_daily |> ggplot() +
        geom_line(aes(plot_date, flow_cfs)) +
        geom_line(aes(plot_date, sim_flow_cfs), color = "steelblue") +
        scale_x_datetime(
          date_labels = "%b", date_breaks = "1 month",
          minor_breaks = NULL
        ) +
        scale_y_continuous(trans = "log10") +
        ggtitle("log scale") +
        xlab("") +
        ylab("Flow [cfs]") +
        theme_minimal()
      pdf(sprintf("%s/daily_cyclical.pdf", plot_path), width = 7, height = 7)
      cyc_g1 <- ggplotGrob(p_cyclical)
      cyc_g2 <- ggplotGrob(p_cyclical_log)
      cyc_g <- rbind(cyc_g1, cyc_g2, size = "first")
      cyc_g$widths <- unit.pmax(cyc_g1$widths, cyc_g2$widths)
      # grid.newpage()
      grid.draw(cyc_g)
      dev.off()

      cyclical_monthly <- x |>
        group_by(month) |>
        summarise(
          nse = NSE(sim_flow_cfs, flow_cfs),
          pbias = pbias(sim_flow_cfs, flow_cfs),
          flow_cfs = mean(flow_cfs),
          sim_flow_cfs = mean(sim_flow_cfs),
          .groups = "drop"
        ) |>
        data.frame() |>
        mutate(
          wyear = ifelse(month >= 10, 1999, 2000),
          plot_date = ISOdatetime(wyear, month, 15, 0, 0, 0)
        )
      pbias_range <- diff(range(cyclical_monthly$pbias, na.rm = T))
      nse_range <- diff(range(cyclical_monthly$nse, na.rm = T))
      scale <- nse_range / pbias_range
      p_cyclical_monthly <- cyclical_monthly |> ggplot() +
        geom_line(aes(plot_date, flow_cfs)) +
        geom_line(aes(plot_date, sim_flow_cfs), color = "steelblue") +
        scale_x_datetime(date_labels = "%b", date_breaks = "1 month") +
        ggtitle("black: obs, blue: sim") +
        xlab("") +
        ylab("Flow [cfs]") +
        theme_minimal()
      p_cyclical_monthly_bias <- ggplot(cyclical_monthly) +
        geom_col(aes(plot_date, pbias)) +
        geom_line(aes(plot_date, nse / scale)) +
        scale_y_continuous(
          sec.axis =
            sec_axis(~ . * scale, name = "NSE (line)")
        ) +
        scale_x_datetime(
          date_labels = "%b", date_breaks = "1 month",
          minor_breaks = NULL
        ) +
        xlab("") +
        ylab("PBIAS (bars) [%]") +
        theme_minimal()
      pdf(sprintf("%s/monthly.pdf", plot_path), width = 7, height = 7)
      mthly_g1 <- ggplotGrob(p_cyclical_monthly)
      mthly_g2 <- ggplotGrob(p_cyclical_monthly_bias)
      mthly_g <- rbind(mthly_g1, mthly_g2, size = "first")
      mthly_g$widths <- unit.pmax(mthly_g1$widths, mthly_g2$widths)
      # grid.newpage()
      grid.draw(mthly_g)
      dev.off()
    }

    output_plots(na.omit(optimal_sim_daily), plot_path)
    if (!por) {
      output_plots(na.omit(valid_sim_daily), plot_path_valid)
    }
  }
} # end multi basin loop

options("cli.num_colors" = NULL)
