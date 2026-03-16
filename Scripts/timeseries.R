############################################################################################
# This script creates ForEdgeClim's model output and observational TOMST TMS-4 sensors for a long timeseries.
#
# Author: Emma Van de Walle - Q-ForestLab
############################################################################################

# ---- Packages ----
suppressPackageStartupMessages({
  library(ForEdgeClimOnPanama)
  library(dplyr)
  library(Metrics)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(future)
  library(rlang)
  library(future.apply)
  library(RANN)
})

# Parallelization on HPC (supercomputer)
cores <- as.numeric(Sys.getenv("PBS_NUM_PPN", unset = 1))
plan(multicore, workers = cores)
print(paste("Using", availableCores(), "cores."))

start_analysis <- Sys.time()

# -----------------
# INPUT
# -----------------

plotnr = "R11" # C13 or R11

params = list(
  # Calibrated parameters (HEAT)
  #infl_macro = 32.01, g_macro = 38.88, infl_soil = 3.31,
  # C13: infl_macro = 16.03, g_macro = 38.87, infl_soil = 3.05,
  # R11: infl_macro = 32.01, g_macro = 38.88, infl_soil = 3.31,
  # uncalibrated, fixed values: infl_macro = 32.5, g_macro = 25, infl_soil = 5,

  # Fixed parameters (SW RTM)
  betad = 0.325, beta0 = 0.325, omega = 0.52, Kd_v = 0.775, Kb_v = 1.25, omega_g_v = 0.13,
  Kd_h = 0.725, Kb_h = 1.15, omega_g_h = 0.15,
  # Fixed parameters (LW RTM)
  e_forest = 0.965, beta_lw = 0.325, omega_lw = 0.035, Kd_lw_v = 0.3, omega_g_lw_v = 0.055, Kd_lw_h = 0.3, omega_g_lw_h = 0.035,
  # Fixed parameters (HEAT)
  h = 10, p_ground = 0.225, k_soil = 1.225,

    # g and i
  g_macro_ref = 25,
  g_forest_ref = 12.5,
  g_soil_ref = 10,
  infl_macro_ref = 32.5,
  infl_soil_ref = 5,
  infl_forest_ref = 5,

    # PM LE parameters
  gs_max = 0.05,
  stomata_closure = 0.6,

    # g(u) and i(u) parameter
  wind_factor = 3.5 # 5.5/10+

)

## timeseries
all_datetimes <- c(
  seq(
    from = as.POSIXct("2023-01-01 00:00:00", tz = "UTC"),
    to   = as.POSIXct("2023-10-31 23:00:00", tz = "UTC"),
    by   = "hour"
  )
)



# plot specifics
if(plotnr == "C13"){
  height_canopy = 33
  scale_factor = 1
  voxel_TLS <- readRDS("Data/TLS_scaled_DTM_and_grid_C13_Panama.rds")
}else{
  height_canopy = 35
  scale_factor = 6.20/6.22
  voxel_TLS <- readRDS("Data/TLS_scaled_DTM_and_grid_R11_Panama.rds")
}



# Parallelize
parts <- future_lapply(
  all_datetimes,
  FUN = function(dt) {

    # make params accessible in the worker-GlobalEnv
    list2env(params, envir = .GlobalEnv)

    # model inputs
    create_input_drivers()
    create_physical_constants()

    datetime <- as.POSIXct(dt, tz = "UTC")

    out <- tryCatch({

      # ---- Imports ----
      import_temp_observations(datetime)
      import_rad_observations(datetime)
      import_RH_u_observations(datetime)

      # ---- Model run ----
      grid_local <- voxel_TLS$grid
      grid_local$density <- grid_local$density * scale_factor
      res <- run_foredgeclim(grid_local, datetime)
      micro_grid <- res$micro_grid
      air_temp   <- res$air_temperature
      if (any(!is.finite(air_temp))) stop("NaN/Inf in air_temp")

      # ---- Observations ----
      key <- format(datetime, "%Y%m%d_%H")
      TOMST <- read.csv(paste0("Data/microtemp_air_obs_", key, "_", plotnr, ".csv"))

      # ---- Model extraction ----
      temp_air_grid <- micro_grid
      temp_air_grid$temperature <- air_temp - 273.15
      temp_air_grid <- temp_air_grid %>%
        filter(z == 1) %>%
        select(x, y, temperature)

      # search nearest neighbour
      # matrix with grid coords
      grid_coords <- as.matrix(temp_air_grid[, c("x","y")])

      # matrix with sensor coords
      sensor_coords <- as.matrix(TOMST[, c("x","y")])

      # nearest neighbour search
      nn <- nn2(
        data = grid_coords,
        query = sensor_coords,
        k = 1
      )

      idx <- nn$nn.idx[, 1]

      # add modelled temperatures + bookkeeping
      TOMST <- TOMST %>%
        dplyr::mutate(
          datetime = datetime,
          temperature_model = temp_air_grid$temperature[idx]
        ) %>%
        dplyr::mutate(
          resid = temperature_model - value_mean
        )


      TOMST

    }, error = function(e) {
      warning(paste("Crash during model run:", format(datetime), "-", e$message))
      NULL
    })
    out
  },
  future.seed = TRUE
)

# combine all hours into one big table
tomst_all <- dplyr::bind_rows(parts)

# save
out_dir  <- file.path("Output", "timeseries")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, paste0("timeseries_", plotnr, ".rds"))
saveRDS(tomst_all, out_file)


plan(sequential)


end_analysis <- Sys.time()
cat(sprintf('\nTotal running time = %.2f min\n', as.numeric(end_analysis - start_analysis, units = 'mins')))
