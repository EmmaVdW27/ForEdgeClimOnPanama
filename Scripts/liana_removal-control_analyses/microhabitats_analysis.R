################################################################################
# ForEdgeClim 3D thermal microhabitat (= niche) analysis
#
# Quantifies the availability, vertical distribution, and spatial connectivity
# of thermal microhabitats within the simulated three-dimensional forest
# microclimate.
#
# For each plot and timestep, the script:
#   - runs ForEdgeClim using the plot-specific TLS forest structure;
#   - calculates the air-temperature anomaly of every 1 m3 voxel relative to
#     the mean temperature of the complete 3D forest volume;
#   - classifies each voxel into one of five thermal microhabitat classes:
#     Very cool, Cool, Near-average, Warm, or Very warm;
#   - quantifies the availability and vertical distribution of each thermal
#     microhabitat;
#   - quantifies spatial connectivity using the six face-sharing neighbours
#     of each voxel;
#   - calculates thermal microhabitat composition within predefined vertical
#     forest zones;
#   - calculates the combined availability of the Very cool and Very warm
#     microhabitats.
#
# Analysis periods:
#   - End dry season:  1 April - 30 April 2023
#   - End wet season: 15 November - 15 December 2023
#
# Daytime and nighttime conditions are both included.
# Every model voxel represents a 1 x 1 x 1 m volume.
#
# Author: Emma Van de Walle – Q-ForestLab
################################################################################


# ---- Packages ----
suppressPackageStartupMessages({
  library(ForEdgeClimOnPanama)
  library(dplyr)
  library(tibble)
  library(tidyr)
  library(future)
  library(future.apply)
  library(progressr)
})


# ==========================================================================================
# PARALLEL
# ==========================================================================================

handlers(global = TRUE)
handlers("txtprogressbar")
options(progressr.enable = TRUE)

cores <- 20

plan(
  multicore,
  workers = cores
)

cat("Using", cores, "cores.\n")

start_analysis <- Sys.time()


# ==========================================================================================
# INPUT
# ==========================================================================================

treatment_string <- "control_calibrated"


# ------------------------------------------------------------------------------------------
# Plot specifics
# ------------------------------------------------------------------------------------------

plot_info <- tibble(

  plotnr = c(
    "C8",
    "C12",
    "C13",
    "C16"
  ),

  canopy_height = c(
    38,
    29,
    33,
    30
  ),

  scale_factor = c(
    5.81 / 6.70,
    6.37 / 6.70,
    6.22 / 6.70,
    5.98 / 6.70
  ),

  tls_file = c(
    "Data/TLS_scaled_DTM_and_grid_C8_Panama.rds",
    "Data/TLS_scaled_DTM_and_grid_C12_Panama.rds",
    "Data/TLS_scaled_DTM_and_grid_C13_Panama.rds",
    "Data/TLS_scaled_DTM_and_grid_C16_Panama.rds"
  )
)


# ==========================================================================================
# TIME
# ==========================================================================================

# ---- End dry season: April 2023 ----

datetimes_end_dry <- seq(

  from = as.POSIXct(
    "2023-04-01 00:00:00",
    tz = "UTC"
  ),

  to = as.POSIXct(
    "2023-04-30 23:00:00",
    tz = "UTC"
  ),

  by = "hour"
)


# ---- End wet season: mid-November to mid-December 2023 ----

datetimes_end_wet <- seq(

  from = as.POSIXct(
    "2023-11-15 00:00:00",
    tz = "UTC"
  ),

  to = as.POSIXct(
    "2023-12-15 23:00:00",
    tz = "UTC"
  ),

  by = "hour"
)


# ---- Combine periods ----

all_datetimes <- c(
  datetimes_end_dry,
  datetimes_end_wet
)

cat(
  "Number of timesteps per plot:",
  length(all_datetimes),
  "\n"
)


# ==========================================================================================
# MODEL PARAMETERS
# ==========================================================================================

params <- list(

  # calibrated values

  infl_macro_ref = 21.0959353,
  infl_soil_ref  = 6.4076838,
  g_macro_ref    = 20.4131636,
  k_soil         = 1.0114076,
  u_ref          = 0.3800997,

  omega = 0.5266289,
  Kb_v  = 0.5058042,


  # fixed values

  betad        = 0.325,
  beta0        = 0.325,
  Kd_v         = 0.775,
  omega_g_v    = 0.13,

  e_forest     = 0.965,
  beta_lw      = 0.325,
  omega_lw     = 0.035,
  Kd_lw_v      = 0.3,
  omega_g_lw_v = 0.055,

  h               = 10,
  p_ground        = 0.225,
  infl_forest_ref = 5,
  g_soil_ref      = 10,
  g_forest_ref    = 12.5,

  gs_max          = 0.004,
  stomata_closure = 0.75,

  k_u = 0.7
)


# ==========================================================================================
# THERMAL NICHE DEFINITIONS
#
# For every timestep:
#
#     T_anomaly = Tair_voxel - mean(Tair_3D_volume)
#
# Thermal niches:
#
# Very cool       T' < -1.0 °C
# Cool            -1.0 <= T' < -0.5 °C
# Near-average    -0.5 <= T' < +0.5 °C
# Warm            +0.5 <= T' < +1.0 °C
# Very warm       T' >= +1.0 °C
# ==========================================================================================

niche_levels <- c(
  "Very cool",
  "Cool",
  "Near-average",
  "Warm",
  "Very warm"
)


classify_thermal_niche <- function(T_anomaly) {

  cut(

    T_anomaly,

    breaks = c(
      -Inf,
      -1,
      -0.5,
      0.5,
      1,
      Inf
    ),

    labels = niche_levels,

    right = FALSE
  )
}


# ==========================================================================================
# SPATIAL NICHE CONNECTIVITY
#
# Connectivity is quantified using the six face-sharing neighbours of every voxel:
#
#       x - 1
#       x + 1
#       y - 1
#       y + 1
#       z - 1
#       z + 1
#
# Diagonal neighbours are not included.
#
# For every voxel:
#
#   n_existing_neighbours
#       number of face-sharing neighbours that exist in the 3D grid
#
#   n_same_niche_neighbours
#       number of those neighbours belonging to the same thermal niche
#
# For every niche:
#
#   mean_same_niche_neighbours
#       mean number of same-niche neighbours per niche voxel
#
#   same_niche_neighbour_fraction
#       total same-niche neighbour connections /
#       total existing neighbour connections for niche voxels
#
# Higher values indicate greater local spatial aggregation / connectivity.
# ==========================================================================================

calculate_niche_connectivity <- function(temp_grid) {


  # ----------------------------------------------------------------------------------------
  # Integer voxel coordinates
  # ----------------------------------------------------------------------------------------

  x0 <- min(
    temp_grid$x,
    na.rm = TRUE
  )

  y0 <- min(
    temp_grid$y,
    na.rm = TRUE
  )

  z0 <- min(
    temp_grid$z,
    na.rm = TRUE
  )


  xi <- as.integer(
    round(
      temp_grid$x - x0
    )
  )

  yi <- as.integer(
    round(
      temp_grid$y - y0
    )
  )

  zi <- as.integer(
    round(
      temp_grid$z - z0
    )
  )


  # ----------------------------------------------------------------------------------------
  # Dimensions of 3D voxel array
  # ----------------------------------------------------------------------------------------

  nx <- max(xi) + 1L
  ny <- max(yi) + 1L
  nz <- max(zi) + 1L


  # ----------------------------------------------------------------------------------------
  # Convert niche factor to integer
  #
  # 0 = no voxel / NA
  # 1 = Very cool
  # 2 = Cool
  # 3 = Near-average
  # 4 = Warm
  # 5 = Very warm
  # ----------------------------------------------------------------------------------------

  niche_id <- as.integer(
    temp_grid$thermal_niche
  )

  niche_array <- array(
    0L,
    dim = c(
      nx,
      ny,
      nz
    )
  )


  niche_array[
    cbind(
      xi + 1L,
      yi + 1L,
      zi + 1L
    )
  ] <- niche_id


  # ----------------------------------------------------------------------------------------
  # Storage
  # ----------------------------------------------------------------------------------------

  n_voxels <- nrow(
    temp_grid
  )

  n_existing_neighbours <- integer(
    n_voxels
  )

  n_same_niche_neighbours <- integer(
    n_voxels
  )


  # ========================================================================================
  # X DIRECTION
  # ========================================================================================

  if (nx > 1) {

    # ---- neighbour at x - 1 ----

    has_neighbour <- xi > 0L

    idx <- which(
      has_neighbour
    )

    neighbour_niche <- niche_array[
      cbind(
        xi[idx],
        yi[idx] + 1L,
        zi[idx] + 1L
      )
    ]

    existing <- neighbour_niche > 0L

    n_existing_neighbours[idx] <-
      n_existing_neighbours[idx] +
      existing

    n_same_niche_neighbours[idx] <-
      n_same_niche_neighbours[idx] +
      (
        existing &
          neighbour_niche == niche_id[idx]
      )


    # ---- neighbour at x + 1 ----

    has_neighbour <- xi < (nx - 1L)

    idx <- which(
      has_neighbour
    )

    neighbour_niche <- niche_array[
      cbind(
        xi[idx] + 2L,
        yi[idx] + 1L,
        zi[idx] + 1L
      )
    ]

    existing <- neighbour_niche > 0L

    n_existing_neighbours[idx] <-
      n_existing_neighbours[idx] +
      existing

    n_same_niche_neighbours[idx] <-
      n_same_niche_neighbours[idx] +
      (
        existing &
          neighbour_niche == niche_id[idx]
      )
  }


  # ========================================================================================
  # Y DIRECTION
  # ========================================================================================

  if (ny > 1) {

    # ---- neighbour at y - 1 ----

    has_neighbour <- yi > 0L

    idx <- which(
      has_neighbour
    )

    neighbour_niche <- niche_array[
      cbind(
        xi[idx] + 1L,
        yi[idx],
        zi[idx] + 1L
      )
    ]

    existing <- neighbour_niche > 0L

    n_existing_neighbours[idx] <-
      n_existing_neighbours[idx] +
      existing

    n_same_niche_neighbours[idx] <-
      n_same_niche_neighbours[idx] +
      (
        existing &
          neighbour_niche == niche_id[idx]
      )


    # ---- neighbour at y + 1 ----

    has_neighbour <- yi < (ny - 1L)

    idx <- which(
      has_neighbour
    )

    neighbour_niche <- niche_array[
      cbind(
        xi[idx] + 1L,
        yi[idx] + 2L,
        zi[idx] + 1L
      )
    ]

    existing <- neighbour_niche > 0L

    n_existing_neighbours[idx] <-
      n_existing_neighbours[idx] +
      existing

    n_same_niche_neighbours[idx] <-
      n_same_niche_neighbours[idx] +
      (
        existing &
          neighbour_niche == niche_id[idx]
      )
  }


  # ========================================================================================
  # Z DIRECTION
  # ========================================================================================

  if (nz > 1) {

    # ---- neighbour at z - 1 ----

    has_neighbour <- zi > 0L

    idx <- which(
      has_neighbour
    )

    neighbour_niche <- niche_array[
      cbind(
        xi[idx] + 1L,
        yi[idx] + 1L,
        zi[idx]
      )
    ]

    existing <- neighbour_niche > 0L

    n_existing_neighbours[idx] <-
      n_existing_neighbours[idx] +
      existing

    n_same_niche_neighbours[idx] <-
      n_same_niche_neighbours[idx] +
      (
        existing &
          neighbour_niche == niche_id[idx]
      )


    # ---- neighbour at z + 1 ----

    has_neighbour <- zi < (nz - 1L)

    idx <- which(
      has_neighbour
    )

    neighbour_niche <- niche_array[
      cbind(
        xi[idx] + 1L,
        yi[idx] + 1L,
        zi[idx] + 2L
      )
    ]

    existing <- neighbour_niche > 0L

    n_existing_neighbours[idx] <-
      n_existing_neighbours[idx] +
      existing

    n_same_niche_neighbours[idx] <-
      n_same_niche_neighbours[idx] +
      (
        existing &
          neighbour_niche == niche_id[idx]
      )
  }


  # ========================================================================================
  # RETURN VOXEL-LEVEL CONNECTIVITY
  # ========================================================================================

  tibble(

    n_existing_neighbours =
      n_existing_neighbours,

    n_same_niche_neighbours =
      n_same_niche_neighbours
  )
}


# ==========================================================================================
# CALCULATE NICHE METRICS
#
# Output:
#
# one row per thermal niche
#
# Metrics:
#
#   niche_availability
#       fraction of complete 3D forest volume belonging to niche
#
#   mean_same_niche_neighbours
#       average number of same-niche face-sharing neighbours
#
#   same_niche_neighbour_fraction
#       fraction of existing neighbour connections that connect to the same niche
#
#   median_height
#       median vertical position of niche
#
#   p10_height / p90_height
#       central 80% vertical distribution
#
#   vertical_range80
#       p90_height - p10_height
# ==========================================================================================

calculate_niche_metrics <- function(temp_grid) {


  total_voxels <- nrow(
    temp_grid
  )


  # ----------------------------------------------------------------------------------------
  # Calculate connectivity once for the complete grid
  # ----------------------------------------------------------------------------------------

  connectivity <-
    calculate_niche_connectivity(
      temp_grid
    )


  temp_grid <- bind_cols(
    temp_grid,
    connectivity
  )


  # ----------------------------------------------------------------------------------------
  # Summarise each thermal niche
  # ----------------------------------------------------------------------------------------

  niche_results <- lapply(

    niche_levels,

    function(current_niche) {


      niche_dat <- temp_grid[

        !is.na(temp_grid$thermal_niche) &

          as.character(
            temp_grid$thermal_niche
          ) == current_niche,

      ]


      n_niche <- nrow(
        niche_dat
      )


      # ====================================================================================
      # NICHE ABSENT
      # ====================================================================================

      if (n_niche == 0) {

        return(

          tibble(

            thermal_niche =
              current_niche,

            n_voxels =
              0L,

            niche_availability =
              0,

            mean_same_niche_neighbours =
              NA_real_,

            same_niche_neighbour_fraction =
              NA_real_,

            median_height =
              NA_real_,

            p10_height =
              NA_real_,

            p90_height =
              NA_real_,

            vertical_range80 =
              NA_real_
          )
        )
      }


      # ====================================================================================
      # SPATIAL CONNECTIVITY
      # ====================================================================================

      mean_same_niche_neighbours <- mean(

        niche_dat$n_same_niche_neighbours,

        na.rm = TRUE
      )


      total_existing_connections <- sum(

        niche_dat$n_existing_neighbours,

        na.rm = TRUE
      )


      total_same_niche_connections <- sum(

        niche_dat$n_same_niche_neighbours,

        na.rm = TRUE
      )


      same_niche_neighbour_fraction <-

        if (
          total_existing_connections > 0
        ) {

          total_same_niche_connections /
            total_existing_connections

        } else {

          NA_real_
        }


      # ====================================================================================
      # VERTICAL POSITION
      # ====================================================================================

      median_height <- median(
        niche_dat$z,
        na.rm = TRUE
      )


      p10_height <- quantile(

        niche_dat$z,

        probs = 0.10,

        na.rm = TRUE,

        names = FALSE
      )


      p90_height <- quantile(

        niche_dat$z,

        probs = 0.90,

        na.rm = TRUE,

        names = FALSE
      )


      # ====================================================================================
      # RETURN
      # ====================================================================================

      tibble(

        thermal_niche =
          current_niche,

        # ---- Availability ----

        n_voxels =
          n_niche,

        niche_availability =
          n_niche /
          total_voxels,

        # ---- Spatial connectivity ----

        mean_same_niche_neighbours =
          mean_same_niche_neighbours,

        same_niche_neighbour_fraction =
          same_niche_neighbour_fraction,

        # ---- Vertical distribution ----

        median_height =
          median_height,

        p10_height =
          p10_height,

        p90_height =
          p90_height,

        vertical_range80 =
          p90_height -
          p10_height
      )
    }
  )


  bind_rows(
    niche_results
  )
}


# ==========================================================================================
# CALCULATE NICHE COMPOSITION WITHIN VERTICAL FOREST ZONES
#
# For each vertical zone:
#
#     fraction = number of voxels belonging to niche
#                ------------------------------------
#                total number of voxels in that zone
#
# Therefore, within each height class:
#
#     sum(niche_fraction_within_zone) = 1
# ==========================================================================================

calculate_vertical_niche_metrics <- function(temp_grid) {


  # ========================================================================================
  # DEFINE HEIGHT CLASSES
  # ========================================================================================

  vertical_dat <- temp_grid %>%

    mutate(

      height_class = cut(

        z,

        breaks = c(
          -Inf,
          2,
          10,
          20,
          30,
          Inf
        ),

        labels = c(
          "0-2 m",
          "2-10 m",
          "10-20 m",
          "20-30 m",
          ">30 m"
        ),

        right = FALSE
      )
    )


  # ========================================================================================
  # COUNT NICHE VOXELS WITHIN EACH HEIGHT CLASS
  # ========================================================================================

  vertical_dat <- vertical_dat %>%

    count(

      height_class,

      thermal_niche,

      name = "n_voxels",

      .drop = FALSE
    )


  # ========================================================================================
  # CALCULATE FRACTION WITHIN EACH VERTICAL ZONE
  # ========================================================================================

  vertical_dat <- vertical_dat %>%

    group_by(
      height_class
    ) %>%

    mutate(

      zone_total_voxels =
        sum(
          n_voxels
        ),

      niche_fraction_within_zone =
        ifelse(

          zone_total_voxels > 0,

          n_voxels /
            zone_total_voxels,

          NA_real_
        )
    ) %>%

    ungroup()


  vertical_dat
}


# ==========================================================================================
# PROCESS ONE DATETIME
# ==========================================================================================

process_datetime <- function(
    dt,
    plotnr,
    voxel_TLS,
    scale_factor
) {


  # ========================================================================================
  # MODEL SETUP
  # ========================================================================================

  list2env(
    params,
    envir = .GlobalEnv
  )

  create_input_drivers()

  create_physical_constants()


  datetime <- as.POSIXct(
    dt,
    tz = "UTC"
  )


  # ========================================================================================
  # RUN
  # ========================================================================================

  out <- tryCatch({


    # --------------------------------------------------------------------------------------
    # IMPORT DRIVERS
    # --------------------------------------------------------------------------------------

    import_temp_observations(
      datetime
    )

    import_rad_observations(
      datetime
    )

    import_RH_u_observations(
      datetime
    )


    # --------------------------------------------------------------------------------------
    # PREPARE TLS GRID
    # --------------------------------------------------------------------------------------

    grid_local <-
      voxel_TLS$grid


    grid_local$density <-
      grid_local$density *
      scale_factor


    # --------------------------------------------------------------------------------------
    # RUN FOREDGECLIM
    # --------------------------------------------------------------------------------------

    res <- run_foredgeclim(

      grid_local,

      datetime
    )


    micro_grid <-
      res$micro_grid


    air_temp <-
      res$air_temperature


    # --------------------------------------------------------------------------------------
    # CHECK OUTPUT
    # --------------------------------------------------------------------------------------

    if (
      any(
        !is.finite(
          air_temp
        )
      )
    ) {

      stop(
        "NaN/Inf in air_temp"
      )
    }


    # ======================================================================================
    # TEMPORARY 3D TEMPERATURE GRID
    # ======================================================================================

    temp_grid <- micro_grid %>%

      transmute(

        x = x,

        y = y,

        z = z,

        Tair =
          air_temp -
          273.15
      )


    # ======================================================================================
    # MEAN TEMPERATURE OF COMPLETE 3D FOREST VOLUME
    # ======================================================================================

    Tair_volume_mean <- mean(

      temp_grid$Tair,

      na.rm = TRUE
    )


    # ======================================================================================
    # TEMPERATURE ANOMALY + NICHE CLASSIFICATION
    # ======================================================================================

    temp_grid <- temp_grid %>%

      mutate(

        T_anomaly =
          Tair -
          Tair_volume_mean,

        thermal_niche =
          classify_thermal_niche(
            T_anomaly
          )
      )


    # ======================================================================================
    # 1. COMPLETE 3D NICHE METRICS
    # ======================================================================================

    niche_metrics <-
      calculate_niche_metrics(
        temp_grid
      )


    niche_metrics <- niche_metrics %>%

      mutate(

        plot =
          plotnr,

        datetime =
          datetime,

        Tair_volume_mean =
          Tair_volume_mean
      ) %>%

      relocate(

        plot,

        datetime,

        Tair_volume_mean,

        thermal_niche
      )


    # ======================================================================================
    # 2. VERTICAL NICHE COMPOSITION
    # ======================================================================================

    vertical_niche_metrics <-
      calculate_vertical_niche_metrics(
        temp_grid
      )


    vertical_niche_metrics <-

      vertical_niche_metrics %>%

      mutate(

        plot =
          plotnr,

        datetime =
          datetime,

        Tair_volume_mean =
          Tair_volume_mean
      ) %>%

      relocate(

        plot,

        datetime,

        Tair_volume_mean,

        height_class,

        thermal_niche
      )


    # ======================================================================================
    # RETURN BOTH OUTPUTS
    # ======================================================================================

    list(

      niche_metrics =
        niche_metrics,

      vertical_niche_metrics =
        vertical_niche_metrics
    )


  }, error = function(e) {


    warning(

      paste(

        "Crash:",

        plotnr,

        format(
          datetime
        ),

        "-",

        e$message
      )
    )


    NULL
  })


  out
}


# ==========================================================================================
# RUN ALL PLOTS
# ==========================================================================================

all_parts <- list()


for (p in seq_len(nrow(plot_info))) {


  info <- plot_info[p, ]

  plotnr <- info$plotnr

  height_canopy <<- info$canopy_height

  scale_factor <- info$scale_factor


  voxel_TLS <- readRDS(
    info$tls_file
  )


  # ========================================================================================
  # PARALLEL RUN OVER ALL DATETIMES
  # ========================================================================================

  cat(
    "\nProcessing plot:",
    plotnr,
    "\n"
  )


  plot_start <- Sys.time()


  with_progress({

    prog <- progressor(
      along = all_datetimes
    )


    parts <- future_lapply(

      seq_along(all_datetimes),

      function(i) {

        prog()


        process_datetime(

          dt =
            all_datetimes[i],

          plotnr =
            plotnr,

          voxel_TLS =
            voxel_TLS,

          scale_factor =
            scale_factor
        )
      },

      future.seed = TRUE
    )
  })


  # ========================================================================================
  # STORE RESULTS
  # ========================================================================================

  all_parts[[plotnr]] <- parts


  plot_end <- Sys.time()


  cat(

    sprintf(

      "Finished %s in %.2f min\n",

      plotnr,

      as.numeric(
        plot_end -
          plot_start,
        units = "mins"
      )
    )
  )
}


# ==========================================================================================
# FLATTEN RESULTS
# ==========================================================================================

results_flat <- unlist(

  all_parts,

  recursive = FALSE
)


# ---- Remove failed simulations ----

results_flat <- results_flat[

  !vapply(

    results_flat,

    is.null,

    logical(1)
  )
]


# ==========================================================================================
# COMBINE NICHE METRICS
# ==========================================================================================

niche_metrics_all <- bind_rows(

  lapply(

    results_flat,

    function(x) {

      x$niche_metrics
    }
  )
)


# ==========================================================================================
# COMBINE VERTICAL NICHE METRICS
# ==========================================================================================

vertical_niche_all <- bind_rows(

  lapply(

    results_flat,

    function(x) {

      x$vertical_niche_metrics
    }
  )
)


# ==========================================================================================
# ADD TIME VARIABLES
# ==========================================================================================

add_time_variables <- function(x) {


  x %>%

    mutate(

      # ---- Date ----

      date = as.Date(
        datetime,
        tz = "UTC"
      ),


      # ---- Hour ----

      hour = as.integer(
        format(
          datetime,
          "%H"
        )
      ),


      # ---- Day of year ----

      doy = as.integer(
        format(
          datetime,
          "%j"
        )
      ),


      # ---- Day / night ----

      period = case_when(

        hour >= 11 &
          hour <= 22 ~ "Day",

        TRUE ~ "Night"
      ),


      # ---- Seasonal analysis period ----

      seasonal_period = case_when(

        date >= as.Date("2023-04-01") &
          date <= as.Date("2023-04-30") ~
          "End dry season",

        date >= as.Date("2023-11-15") &
          date <= as.Date("2023-12-15") ~
          "End wet season",

        TRUE ~ NA_character_
      )
    )
}


niche_metrics_all <-
  add_time_variables(
    niche_metrics_all
  )


vertical_niche_all <-
  add_time_variables(
    vertical_niche_all
  )


# ==========================================================================================
# ADD EXTREME NICHE AVAILABILITY
#
# Fraction of complete 3D forest volume that is either:
#
#       Very cool (< -1 C)
#       OR
#       Very warm (> +1 C)
#
# One row per plot x datetime.
# ==========================================================================================

extreme_niche_all <- niche_metrics_all %>%

  filter(
    thermal_niche %in%
      c(
        "Very cool",
        "Very warm"
      )
  ) %>%

  group_by(
    plot,
    datetime,
    Tair_volume_mean,
    date,
    hour,
    doy,
    period,
    seasonal_period
  ) %>%

  summarise(

    extreme_niche_availability =
      sum(
        niche_availability,
        na.rm = TRUE
      ),

    .groups = "drop"
  )


# ==========================================================================================
# FACTOR ORDER
# ==========================================================================================

niche_metrics_all <- niche_metrics_all %>%

  mutate(

    thermal_niche = factor(

      thermal_niche,

      levels = niche_levels
    ),

    period = factor(

      period,

      levels = c(
        "Day",
        "Night"
      )
    ),

    seasonal_period = factor(

      seasonal_period,

      levels = c(
        "End dry season",
        "End wet season"
      )
    )
  )


vertical_niche_all <- vertical_niche_all %>%

  mutate(

    thermal_niche = factor(

      thermal_niche,

      levels = niche_levels
    ),

    height_class = factor(

      height_class,

      levels = c(
        "0-2 m",
        "2-10 m",
        "10-20 m",
        "20-30 m",
        ">30 m"
      )
    ),

    period = factor(

      period,

      levels = c(
        "Day",
        "Night"
      )
    ),

    seasonal_period = factor(

      seasonal_period,

      levels = c(
        "End dry season",
        "End wet season"
      )
    )
  )


extreme_niche_all <- extreme_niche_all %>%

  mutate(

    period = factor(

      period,

      levels = c(
        "Day",
        "Night"
      )
    ),

    seasonal_period = factor(

      seasonal_period,

      levels = c(
        "End dry season",
        "End wet season"
      )
    )
  )


# ==========================================================================================
# OUTPUT DIRECTORY
# ==========================================================================================

dir.create(

  "Output/thermal_niches",

  recursive = TRUE,

  showWarnings = FALSE
)


# ==========================================================================================
# SAVE 1
#
# COMPLETE 3D NICHE METRICS
#
# Contains:
#
# plot
# datetime
# Tair_volume_mean
# thermal_niche
#
# n_voxels
# niche_availability
#
# mean_same_niche_neighbours
# same_niche_neighbour_fraction
#
# median_height
# p10_height
# p90_height
# vertical_range80
#
# date
# hour
# doy
# period
# seasonal_period
# ==========================================================================================

saveRDS(

  niche_metrics_all,

  paste0(

    "Output/thermal_niches/thermal_niche_metrics_",

    treatment_string,

    ".rds"
  )
)


# ==========================================================================================
# SAVE 2
#
# NICHE COMPOSITION WITHIN VERTICAL FOREST ZONES
# ==========================================================================================

saveRDS(

  vertical_niche_all,

  paste0(

    "Output/thermal_niches/thermal_niche_vertical_",

    treatment_string,

    ".rds"
  )
)


# ==========================================================================================
# SAVE 3
#
# EXTREME NICHE AVAILABILITY
# ==========================================================================================

saveRDS(

  extreme_niche_all,

  paste0(

    "Output/thermal_niches/thermal_niche_extreme_",

    treatment_string,

    ".rds"
  )
)


# ==========================================================================================
# CLEANUP
# ==========================================================================================

plan(
  sequential
)


end_analysis <- Sys.time()


cat(

  sprintf(

    "\nTotal running time = %.2f min\n",

    as.numeric(

      end_analysis -
        start_analysis,

      units = "mins"
    )
  )
)


# ==========================================================================================
# CONSOLE SUMMARY
# ==========================================================================================

cat(
  "\nSaved:\n"
)


cat(

  paste0(

    "  Output/thermal_niches/thermal_niche_metrics_",

    treatment_string,

    ".rds\n"
  )
)


cat(

  paste0(

    "  Output/thermal_niches/thermal_niche_vertical_",

    treatment_string,

    ".rds\n"
  )
)


cat(

  paste0(

    "  Output/thermal_niches/thermal_niche_extreme_",

    treatment_string,

    ".rds\n"
  )
)
