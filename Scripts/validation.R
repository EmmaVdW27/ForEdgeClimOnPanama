############################################################################################
# This script validates ForEdgeClim's model output with observational TOMST TMS-4 sensors for the whole year.
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

plotnr = "C13" # C13 or R11

season_string = "C13_uncalibrated" # C13 or R11 | _calibrated or _uncalibrated

params = list(
  # (Un)calibrated parameters
  infl_macro_ref = 32.5, infl_soil_ref = 5, g_macro_ref = 25, e_forest = 0.965, u_ref = 1.6, k_soil = 1.225,
  # C13:
  # R11:

  # Fixed parameters (SW RTM)
  betad = 0.325, beta0 = 0.325, omega = 0.52, Kd_v = 0.775, Kb_v = 1.25, omega_g_v = 0.13,
  # Fixed parameters (LW RTM)
  beta_lw = 0.325, omega_lw = 0.035, Kd_lw_v = 0.3, omega_g_lw_v = 0.055,
  # Fixed parameters (HEAT)
  h = 10, p_ground = 0.225,

  # g and i
  g_forest_ref = 12.5,
  g_soil_ref = 10,
  infl_forest_ref = 5,

  # PM LE and wind parameters
  gs_max = 0.004,
  stomata_closure = 0.75,
  k_u = 0.7


)


## year
all_datetimes <- c(
  seq(
    from = as.POSIXct("2023-01-01 00:00:00", tz = "UTC"),
    to   = as.POSIXct("2023-12-31 23:00:00", tz = "UTC"),
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




# -----------------
# VALIDATION: STATISTICAL COMPARISON MODEL VS OBSERVATIONS
# -----------------
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

      # add modelled temperatures to TOMST frame
      TOMST$temperature_model <- temp_air_grid$temperature[nn$nn.idx[,1]]


      diffs <- TOMST$value_mean - TOMST$temperature_model
      valid <- is.finite(diffs)
      if (!any(valid)) stop("No valid differences")

      list(
        sse = sum(diffs[valid]^2),
        n   = sum(valid),
        data = tibble::tibble(
          datetime = datetime,
          obs = TOMST$value_mean[valid],
          mod = TOMST$temperature_model[valid],
          resid = TOMST$temperature_model[valid] - TOMST$value_mean[valid],
          TMS_x = TOMST$x[valid],
          TMS_y = TOMST$y[valid]
        )
      )

    }, error = function(e) {
      warning(paste("Crash during model run:", format(datetime), "-", e$message))
      list(sse = NA_real_, n = 0L)
    })
    out
  },
  future.seed = TRUE
)

# combine data
val_data <- dplyr::bind_rows(lapply(parts, `[[`, "data"))

# save validated data
out_dir  <- file.path("Output", "validation", "data")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(out_dir, paste0("val_data_", season_string, ".rds"))
saveRDS(val_data, out_file)

# calculate statistics
mean_error <- mean(val_data$resid, na.rm = TRUE)
sd_resid   <- sd(val_data$resid, na.rm = TRUE)
rmse       <- sqrt(mean(val_data$resid^2, na.rm = TRUE))
r2         <- cor(val_data$obs, val_data$mod)^2
nse        <- 1 - sum(val_data$resid^2) / sum((val_data$obs - mean(val_data$obs))^2)



plan(sequential)

# Console summaries
cat("\nValidation statistics:\n")
cat(sprintf("Mean Error (Bias): %.2f °C\n", mean_error))
cat(sprintf("SD (residuals): %.2f °C\n", sd_resid))
cat(sprintf("RMSE: %.2f °C\n", rmse))
cat(sprintf("R²: %.2f\n", r2))
cat(sprintf("NSE: %.2f\n", nse))




end_analysis <- Sys.time()
cat(sprintf('\nTotal running time = %.2f min\n', as.numeric(end_analysis - start_analysis, units = 'mins')))
