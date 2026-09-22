################################################################################
# ForEdgeClim annual simulation and air-temperature validation
#
# Runs ForEdgeClim for one forest plot for every hourly timestep of 2023 using
# the calibrated shortwave-radiation and microclimate parameter sets.
#
# For each timestep, the script:
#   - runs ForEdgeClim using the plot-specific TLS forest structure;
#   - extracts mean microclimate and energy-balance variables at 1 m height
#     (Tair, Tsurf, LE, RH, VPD, and direct/diffuse shortwave radiation);
#   - matches simulated air temperature to TOMST sensor locations using
#     nearest-neighbour matching;
#   - stores the resulting annual model timeseries.
#
# Air-temperature performance is evaluated across all available hourly TOMST
# observations using mean error (ME), residual standard deviation (SD), RMSE,
# R2, and Nash-Sutcliffe efficiency (NSE).
#
# Author: Emma Van de Walle – Q-ForestLab
################################################################################

# ---- Packages ----
suppressPackageStartupMessages({
  library(ForEdgeClimOnPanama)
  library(dplyr)
  library(Metrics)
  library(tibble)
  library(tidyr)
  library(ggplot2)
  library(future)
  library(future.apply)
  library(RANN)
  library(progressr)
})

# ---- Parallel ----
handlers(global = TRUE)
handlers("txtprogressbar")
options(progressr.enable = TRUE)

cores <- 50
plan(multicore, workers = cores)
print(paste("Using", cores, "cores."))

start_analysis <- Sys.time()

# -----------------
# INPUT
# -----------------

plotnr = "C13" # C13 or R11
season_string = "C13_calibrated"

params = list(
  # uncalibrated values:
  # infl_macro_ref = 32.5, infl_soil_ref = 5, g_macro_ref = 25, k_soil = 1.225, u_ref = 1.6,
  # omega = 0.52, Kb_v = 1.25,

  # calibrated values:
  infl_macro_ref = 21.0959353, infl_soil_ref = 6.4076838, g_macro_ref = 20.4131636, k_soil = 1.0114076, u_ref = 0.3800997,
  omega = 0.5266289, Kb_v = 0.5058042,

  # fixed values:
  betad = 0.325, beta0 = 0.325,
  Kd_v = 0.775, omega_g_v = 0.13,

  e_forest = 0.965, beta_lw = 0.325, omega_lw = 0.035,
  Kd_lw_v = 0.3, omega_g_lw_v = 0.055,

  h = 10,
  p_ground = 0.225,

  infl_forest_ref = 5,
  g_soil_ref = 10,
  g_forest_ref = 12.5,

  gs_max = 0.004, stomata_closure = 0.75,
  k_u = 0.7
)

# ---- Time ----
all_datetimes <- seq(
  from = as.POSIXct("2023-01-01 00:00:00", tz = "UTC"),
  to   = as.POSIXct("2023-12-31 23:00:00", tz = "UTC"),
  by   = "hour"
)

# plot specifics
if(plotnr == "C13"){
  plot_map = "C13"
  height_canopy = 33
  scale_factor = 1
  voxel_TLS <- readRDS("Data/TLS_scaled_DTM_and_grid_C13_Panama.rds")
} else { # plotnr == "R11"
  plot_map = "R11"
  height_canopy = 35
  scale_factor = 6.20/6.22
  voxel_TLS <- readRDS("Data/TLS_scaled_DTM_and_grid_R11_Panama.rds")
}

# -----------------
# MAIN LOOP
# -----------------
process_datetime <- function(dt, i) {

  list2env(params, envir = .GlobalEnv)
  create_input_drivers()
  create_physical_constants()

  datetime <- as.POSIXct(dt, tz = "UTC")

  out <- tryCatch({

    # ---- Imports ----
    import_temp_observations(datetime)
    import_rad_observations(datetime)
    import_RH_u_observations(datetime)

    # ---- Model ----
    grid_local <- voxel_TLS$grid
    grid_local$density <- grid_local$density * scale_factor

    res <- run_foredgeclim(grid_local, datetime)

    micro_grid <- res$micro_grid
    air_temp   <- res$air_temperature

    if (any(!is.finite(air_temp))) stop("NaN/Inf in air_temp")

    # ---- Extract z = 1 ----
    sel <- micro_grid$z == 1

    surf_temp_mean <- mean(micro_grid$temperature[sel] - 273.15, na.rm = TRUE)
    LE_mean <- mean(res$latent_flux[sel], na.rm = TRUE)
    RH_mean <- mean(res$RH[sel], na.rm = TRUE)
    VPD_mean <- mean(res$VPD[sel], na.rm = TRUE)
    Fd_mean <- mean(res$sw_rad$F_d_down[sel], na.rm = TRUE)
    Fb_mean <- mean(res$sw_rad$F_b_down[sel], na.rm = TRUE)

    # ---- Observations ----
    key <- format(datetime, "%Y%m%d_%H")
    TOMST <- read.csv(paste0("Data/microtemp_air_obs_", key, "_", plotnr, ".csv"))

    # ---- Model grid ----
    temp_air_grid <- micro_grid %>%
      mutate(temperature = air_temp - 273.15) %>%
      filter(z == 1) %>%
      select(x, y, temperature)

    # ---- Mean model temperature at 1 m height ----
    Tair_mean_1m <- mean(temp_air_grid$temperature, na.rm = TRUE)

    # ---- Nearest neighbour ----
    nn <- nn2(
      data = as.matrix(temp_air_grid[, c("x","y")]),
      query = as.matrix(TOMST[, c("x","y")]),
      k = 1
    )

    idx <- nn$nn.idx[,1]

    # ---- Combine ----
    TOMST <- TOMST %>%
      mutate(
        datetime = datetime,
        temperature_model = temp_air_grid$temperature[idx],
        resid = temperature_model - value_mean,
        Tair_mean_1m = Tair_mean_1m,
        LE_mean = LE_mean,
        RH_mean = RH_mean,
        VPD_mean = VPD_mean,
        Fd_down_mean = Fd_mean,
        Fb_down_mean = Fb_mean,
        Tsurf_mean = surf_temp_mean
      )

    # ---- Validation data ----
    diffs <- TOMST$value_mean - TOMST$temperature_model
    valid <- is.finite(diffs)

    list(
      timeseries = TOMST,
      validation = tibble(
        datetime = datetime,
        obs = TOMST$value_mean[valid],
        mod = TOMST$temperature_model[valid],
        resid = TOMST$resid[valid],
        TMS_x = TOMST$x[valid],
        TMS_y = TOMST$y[valid]
      )
    )

  }, error = function(e) {
    warning(paste("Crash:", format(datetime), "-", e$message))
    NULL
  })

  out
}

# run parallel, but keep index
with_progress({
  p <- progressor(along = all_datetimes)

  parts <- future_lapply(
    seq_along(all_datetimes),
    FUN = function(i) {
      p()  # <-- update progress
      process_datetime(all_datetimes[i], i)
    },
    future.seed = TRUE
  )
})

# -----------------
# COMBINE OUTPUT
# -----------------

timeseries_all <- bind_rows(lapply(parts, `[[`, "timeseries"))
val_data       <- bind_rows(lapply(parts, `[[`, "validation"))

# ---- Save ----
dir.create("Output/timeseries", recursive = TRUE, showWarnings = FALSE)
saveRDS(timeseries_all, paste0("Output/timeseries/timeseries_", season_string, ".rds"))

dir.create("Output/validation/data", recursive = TRUE, showWarnings = FALSE)
saveRDS(val_data, paste0("Output/validation/data/val_data_", season_string, ".rds"))

# -----------------
# VALIDATION STATS
# -----------------

mean_error <- mean(val_data$resid, na.rm = TRUE)
sd_resid   <- sd(val_data$resid, na.rm = TRUE)
rmse       <- sqrt(mean(val_data$resid^2, na.rm = TRUE))
r2         <- cor(val_data$obs, val_data$mod, use = "complete.obs")^2
nse        <- 1 - sum(val_data$resid^2, na.rm = TRUE) /
  sum((val_data$obs - mean(val_data$obs, na.rm = TRUE))^2, na.rm = TRUE)

cat("\nValidation statistics:\n")
cat(sprintf("Mean Error: %.2f °C\n", mean_error))
cat(sprintf("SD: %.2f °C\n", sd_resid))
cat(sprintf("RMSE: %.2f °C\n", rmse))
cat(sprintf("R²: %.2f\n", r2))
cat(sprintf("NSE: %.2f\n", nse))

# ---- Cleanup ----
plan(sequential)

end_analysis <- Sys.time()
cat(sprintf('\nTotal running time = %.2f min\n',
            as.numeric(end_analysis - start_analysis, units = 'mins')))
