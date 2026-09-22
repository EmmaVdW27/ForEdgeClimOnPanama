################################################################################
#
# ForEdgeClim - Annual timeseries simulation - treatment plots
#
# Runs ForEdgeClim for all plots within a specified treatment at every hourly
# timestep of 2023 using the calibrated parameter set.
#
# For each plot and timestep, the script:
#   - runs ForEdgeClim using the plot-specific TLS-derived forest structure;
#   - calculates spatial mean air temperature, surface temperature, latent heat
#     flux, relative humidity, vapour pressure deficit, and direct and diffuse
#     downward shortwave radiation at 1 m height;
#   - matches simulated air temperature at 1 m height to TOMST sensor locations
#     using nearest-neighbour matching in the horizontal plane;
#   - calculates model-observation air temperature residuals at TOMST locations;
#   - stores TOMST observations, matched model temperatures, and spatially
#     averaged model outputs in a combined annual timeseries.
#
# For C12 and C16, where the TLS-derived model domain represents only the
# central 40 x 40 m of the original 60 x 60 m plot, TOMST sensors outside this
# domain are excluded and the remaining sensor coordinates are transformed to
# the local model coordinate system.
#
# Output:
#   Output/timeseries/timeseries_<treatment>_calibrated.rds
#
# Author: Emma Van de Walle - Q-ForestLab
#
################################################################################

# ---- Packages ----
suppressPackageStartupMessages({
  library(ForEdgeClimOnPanama)
  library(dplyr)
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

cores <- 20
plan(multicore, workers = cores)
print(paste("Using", cores, "cores."))

start_analysis <- Sys.time()

# ---------------------------------------------------------------
# INPUT
# ---------------------------------------------------------------

season_string <- "control_calibrated"

# plot specifics
plot_info <- tibble(
  plotnr = c("C8","C12","C13","C16"),

  canopy_height = c(38, 29, 33, 30),

  scale_factor = c(5.81/6.70, 6.37/6.70, 6.22/6.70, 5.98/6.70),

  tls_file = c(
    "Data/TLS_scaled_DTM_and_grid_C8_Panama.rds",
    "Data/TLS_scaled_DTM_and_grid_C12_Panama.rds",
    "Data/TLS_scaled_DTM_and_grid_C13_Panama.rds",
    "Data/TLS_scaled_DTM_and_grid_C16_Panama.rds"
  )
)

# ---------------------------------------------------------------

# ---- Time ----
all_datetimes <- seq(
  from = as.POSIXct("2023-01-01 00:00:00", tz = "UTC"),
  to   = as.POSIXct("2023-12-31 23:00:00", tz = "UTC"),
  by   = "hour"
)

params = list(
  # calibrated values:
  infl_macro_ref = 21.0959353, infl_soil_ref = 6.4076838, g_macro_ref = 20.4131636, k_soil = 1.0114076, u_ref = 0.3800997,
  omega = 0.5266289, Kb_v = 0.5058042,
  # fixed values:
  betad = 0.325, beta0 = 0.325, Kd_v = 0.775, omega_g_v = 0.13,
  e_forest = 0.965, beta_lw = 0.325, omega_lw = 0.035, Kd_lw_v = 0.3, omega_g_lw_v = 0.055,
  h = 10, p_ground = 0.225, infl_forest_ref = 5, g_soil_ref = 10, g_forest_ref = 12.5,
  gs_max = 0.004, stomata_closure = 0.75,
  k_u = 0.7
)

# -----------------
# MAIN LOOP
# -----------------
process_datetime <- function(dt, plotnr, voxel_TLS, scale_factor) {

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

    # C12 and C16: keep only sensors in the central 40 x 40 m plot
    # and convert coordinates to the local model grid
    if (plotnr %in% c("C12", "C16")) {

      TOMST <- TOMST %>%
        filter(
          x >= 10, x <= 50,
          y >= 10, y <= 50
        ) %>%
        mutate(
          x = x - 10,
          y = y - 10
        )
    }

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
        plot = plotnr,
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

    list(
      timeseries = TOMST
    )

  }, error = function(e) {
    warning(paste("Crash:", format(datetime), "-", e$message))
    NULL
  })

  out
}

# run parallel
all_parts <- list()

for (p in seq_len(nrow(plot_info))) {

  info <- plot_info[p, ]

  plotnr <- info$plotnr

  height_canopy <<- info$canopy_height

  scale_factor <- info$scale_factor

  voxel_TLS <- readRDS(info$tls_file)

  with_progress({

    prog <- progressor(along = all_datetimes)

    parts <- future_lapply(

      seq_along(all_datetimes),

      function(i){

        prog()

        process_datetime(
          dt = all_datetimes[i],
          plotnr = plotnr,
          voxel_TLS = voxel_TLS,
          scale_factor = scale_factor
        )

      },

      future.seed = TRUE

    )

  })

  all_parts[[plotnr]] <- parts

}

# -----------------
# COMBINE OUTPUT
# -----------------

timeseries_all <- bind_rows(
  lapply(
    unlist(all_parts, recursive = FALSE),
    `[[`,
    "timeseries"
  )
)

# ---- Save ----
dir.create("Output/timeseries", recursive = TRUE, showWarnings = FALSE)
saveRDS(timeseries_all, paste0("Output/timeseries/timeseries_", season_string, ".rds"))



# ---- Cleanup ----
plan(sequential)

end_analysis <- Sys.time()
cat(sprintf('\nTotal running time = %.2f min\n',
            as.numeric(end_analysis - start_analysis, units = 'mins')))
