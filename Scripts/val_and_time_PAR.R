############################################################################################
# Combined script: ForEdgeClim timeseries + validation + txtProgressBar
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
  library(future.apply)
  library(RANN)
  library(progressr)
})

# ---- Parallel ----
handlers(global = TRUE)
handlers("txtprogressbar")
options(progressr.enable = TRUE)
cores <- as.numeric(Sys.getenv("PBS_NUM_PPN", unset = 1))
plan(multisession, workers = cores)

print(paste("Using", availableCores(), "cores."))

start_analysis <- Sys.time()

# -----------------
# INPUT
# -----------------

plotnr = "C8"
season_string = "C8_50_gens_3_params_uncalibrated"

params = list(
  # calibrated values:
  beta0 = 0.325, Kd_v = 0.775, Kb_v = 1.25,

  # fixed values:
  betad = 0.325, omega = 0.52,
  omega_g_v = 0.13,

  e_forest = 0.965, beta_lw = 0.325, omega_lw = 0.035,
  Kd_lw_v = 0.3, omega_g_lw_v = 0.055,

  h = 10,
  p_ground = 0.225,
  k_soil = 1.225,

  infl_macro_ref = 32.5,
  infl_soil_ref = 5,
  infl_forest_ref = 5,
  g_macro_ref = 25,
  g_soil_ref = 10,
  g_forest_ref = 12.5,

  gs_max = 0.004,
  stomata_closure = 0.75,

  k_u = 0.7,
  u_ref = 1.6
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
} else if(plotnr == "R11"){
  plot_map = "R11"
  height_canopy = 35
  scale_factor = 6.20/6.22
  voxel_TLS <- readRDS("Data/TLS_scaled_DTM_and_grid_R11_Panama.rds")
} else if(plotnr == "C8"){
  plot_map = "C8"
  height_canopy = 38
  scale_factor = 5.81/6.22
  voxel_TLS <- readRDS("Data/TLS_scaled_DTM_and_grid_C8_Panama.rds")
} else { # plotnr == "R7"
  plot_map = "R7"
  height_canopy = 37
  scale_factor = 6.19/6.22
  voxel_TLS <- readRDS("Data/TLS_scaled_DTM_and_grid_R7_Panama.rds")
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

    # ---- Extract z = 1 ----
    micro_grid <- res$micro_grid
    sel <- micro_grid$z == 1

    surf_temp_mean <- mean(micro_grid$temperature[sel] - 273.15, na.rm = TRUE)
    LE_mean <- mean(res$latent_flux[sel], na.rm = TRUE)
    RH_mean <- mean(res$RH[sel], na.rm = TRUE)
    VPD_mean <- mean(res$VPD[sel], na.rm = TRUE)
    Fd_mean <- mean(res$sw_rad$F_d_down[sel], na.rm = TRUE)
    Fb_mean <- mean(res$sw_rad$F_b_down[sel], na.rm = TRUE)

    # ---- Observations ----
    key <- format(datetime, "%Y%m%d_%H")
    file_par <- paste0(
      "Data/PAR_obs_",
      key,
      "_",
      plotnr,
      ".csv"
    )

    if(!file.exists(file_par))
      return(NULL)

    PAR_obs <- read.csv(file_par)

    # ---- Model grid ----
    par_grid <- res$sw_rad

    par_grid <- par_grid %>%
      filter(Z == 1) %>%
      mutate(
        PAR_model =
          (F_d_down + F_b_down) * 2.02
      ) %>%
      select(
        X,
        Y,
        PAR_model
      )

    # ---- Nearest neighbour ----
    nn <- nn2(
      data = as.matrix(par_grid[, c("X","Y")]),
      query = as.matrix(PAR_obs[, c("x","y")]),
      k = 1
    )

    idx <- nn$nn.idx[,1]

    # ---- Combine ----
    PAR_obs <- PAR_obs %>%
      mutate(
        datetime = datetime,
        PAR_model = par_grid$PAR_model[idx],
        resid = PAR_model - radiation_mean,
        LE_mean = LE_mean,
        RH_mean = RH_mean,
        VPD_mean = VPD_mean,
        Fd_down_mean = Fd_mean,
        Fb_down_mean = Fb_mean,
        Tsurf_mean = surf_temp_mean
      )

    # ---- Validation data ----
    valid <-
      is.finite(PAR_obs$radiation_mean) &
      is.finite(PAR_obs$PAR_model)
    n_valid <- sum(valid)
    if(n_valid == 0)
      return(NULL)

    list(
      timeseries = PAR_obs,
      validation = tibble(
        datetime = datetime,
        obs = PAR_obs$radiation_mean[valid],
        mod = PAR_obs$PAR_model[valid],
        resid = PAR_obs$resid[valid],
        sensor_x = PAR_obs$x[valid],
        sensor_y = PAR_obs$y[valid]
      )
    )

  }, error = function(e) {
    warning(paste("Crash:", format(datetime), "-", e$message))
    NULL
  })

  out
}

# run parallel, maar index meegeven
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
cat(sprintf("Mean Error: %.2f µmol m-2 s-1\n", mean_error))
cat(sprintf("SD: %.2f µmol m-2 s-1\n", sd_resid))
cat(sprintf("RMSE: %.2f µmol m-2 s-1\n", rmse))
cat(sprintf("R²: %.2f\n", r2))
cat(sprintf("NSE: %.2f\n", nse))

# ---- Cleanup ----
plan(sequential)

end_analysis <- Sys.time()
cat(sprintf('\nTotal running time = %.2f min\n',
            as.numeric(end_analysis - start_analysis, units = 'mins')))
