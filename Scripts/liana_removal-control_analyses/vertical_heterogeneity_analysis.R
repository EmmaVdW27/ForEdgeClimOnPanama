################################################################################
# ForEdgeClim annual 3D thermal heterogeneity analysis
#
# Quantifies the spatial heterogeneity of simulated air temperature throughout
# the three-dimensional forest microclimate for every hourly timestep of 2023.
#
# For each plot and timestep, the script:
#   - runs ForEdgeClim using the plot-specific TLS forest structure;
#   - extracts the simulated 3D air temperature field;
#   - calculates spatial temperature heterogeneity separately within each
#     horizontal 1 m forest layer;
#   - quantifies overall heterogeneity using the standard deviation, interquartile
#     range, 90% spatial range (P95 - P05), and mean absolute anomaly;
#   - quantifies fine-scale thermal heterogeneity from the temperature contrast
#     between individual voxels and their horizontal neighbours.
#
# Only layer-level summary statistics are retained; individual voxel
# temperatures are not stored. The analysis uses model output only and does
# not include observational validation.
#
# Output:
#   Output/heterogeneity/heterogeneity_<treatment>_calibrated.rds
#
# Author: Emma Van de Walle – Q-ForestLab
################################################################################


# ==========================================================================================
# PACKAGES
# ==========================================================================================

suppressPackageStartupMessages({
  library(ForEdgeClimOnPanama)
  library(dplyr)
  library(tibble)
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
plan(multicore, workers = cores)

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

all_datetimes <- seq(

  from = as.POSIXct(
    "2023-01-01 00:00:00",
    tz = "UTC"
  ),

  to = as.POSIXct(
    "2023-12-31 23:00:00",
    tz = "UTC"
  ),

  by = "hour"
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
# FUNCTION: LOCAL THERMAL CONTRAST
# ==========================================================================================

# For each voxel:
#
# local contrast = Tair voxel - mean(Tair of its 8 horizontal neighbours)
#
# This quantifies fine-scale spatial temperature differentiation within
# each horizontal 1 m layer.


calculate_local_contrast <- function(x, y, temperature) {


  # ----------------------------------------------------------------------------------------
  # Coordinates
  # ----------------------------------------------------------------------------------------

  x_unique <- sort(unique(x))
  y_unique <- sort(unique(y))

  nx <- length(x_unique)
  ny <- length(y_unique)


  ix <- match(x, x_unique)
  iy <- match(y, y_unique)


  # ----------------------------------------------------------------------------------------
  # Temperature matrix
  # ----------------------------------------------------------------------------------------

  temp_matrix <- matrix(
    NA_real_,
    nrow = ny,
    ncol = nx
  )


  temp_matrix[
    cbind(iy, ix)
  ] <- temperature


  # ----------------------------------------------------------------------------------------
  # Neighbour sum and number of neighbours
  # ----------------------------------------------------------------------------------------

  neighbour_sum <- matrix(
    0,
    nrow = ny,
    ncol = nx
  )


  neighbour_n <- matrix(
    0,
    nrow = ny,
    ncol = nx
  )


  # ----------------------------------------------------------------------------------------
  # Eight surrounding cells
  # ----------------------------------------------------------------------------------------

  offsets <- expand.grid(
    dy = -1:1,
    dx = -1:1
  )


  offsets <- offsets[
    !(offsets$dx == 0 & offsets$dy == 0),
  ]


  # ----------------------------------------------------------------------------------------
  # Loop over neighbour directions
  # ----------------------------------------------------------------------------------------

  for (k in seq_len(nrow(offsets))) {


    dx <- offsets$dx[k]
    dy <- offsets$dy[k]


    # Define source and target indices

    source_rows <- max(1, 1 - dy):min(ny, ny - dy)
    source_cols <- max(1, 1 - dx):min(nx, nx - dx)

    target_rows <- source_rows + dy
    target_cols <- source_cols + dx


    neighbour <-
      temp_matrix[
        target_rows,
        target_cols,
        drop = FALSE
      ]


    valid <- is.finite(neighbour)


    neighbour_sum[
      source_rows,
      source_cols
    ] <-
      neighbour_sum[
        source_rows,
        source_cols
      ] +
      ifelse(
        valid,
        neighbour,
        0
      )


    neighbour_n[
      source_rows,
      source_cols
    ] <-
      neighbour_n[
        source_rows,
        source_cols
      ] +
      valid
  }


  # ----------------------------------------------------------------------------------------
  # Mean neighbouring temperature
  # ----------------------------------------------------------------------------------------

  neighbour_mean <-
    neighbour_sum /
    neighbour_n


  neighbour_mean[
    neighbour_n == 0
  ] <- NA_real_


  # ----------------------------------------------------------------------------------------
  # Local contrast
  # ----------------------------------------------------------------------------------------

  local_contrast <-
    temp_matrix -
    neighbour_mean


  # Return values in original voxel order

  local_contrast[
    cbind(
      iy,
      ix
    )
  ]
}


# ==========================================================================================
# FUNCTION: HETEROGENEITY FOR ONE HORIZONTAL LAYER
# ==========================================================================================

calculate_layer_heterogeneity <- function(dat) {


  # ----------------------------------------------------------------------------------------
  # Mean temperature
  # ----------------------------------------------------------------------------------------

  T_mean <-
    mean(
      dat$Tair,
      na.rm = TRUE
    )


  # ----------------------------------------------------------------------------------------
  # Temperature anomaly relative to horizontal layer mean
  # ----------------------------------------------------------------------------------------

  T_anomaly <-
    dat$Tair -
    T_mean


  # ----------------------------------------------------------------------------------------
  # Local thermal contrast
  # ----------------------------------------------------------------------------------------

  local_contrast <-
    calculate_local_contrast(

      x =
        dat$x,

      y =
        dat$y,

      temperature =
        dat$Tair
    )


  # ----------------------------------------------------------------------------------------
  # Summary statistics
  # ----------------------------------------------------------------------------------------

  tibble(

    n_voxels =
      sum(
        is.finite(dat$Tair)
      ),


    # --------------------------------------------------------------------------------------
    # Mean state
    # --------------------------------------------------------------------------------------

    Tair_mean =
      T_mean,


    # --------------------------------------------------------------------------------------
    # Overall spatial temperature heterogeneity
    # --------------------------------------------------------------------------------------

    Tair_sd =
      sd(
        dat$Tair,
        na.rm = TRUE
      ),


    Tair_p05 =
      quantile(
        dat$Tair,
        0.05,
        na.rm = TRUE,
        names = FALSE
      ),


    Tair_p25 =
      quantile(
        dat$Tair,
        0.25,
        na.rm = TRUE,
        names = FALSE
      ),


    Tair_p75 =
      quantile(
        dat$Tair,
        0.75,
        na.rm = TRUE,
        names = FALSE
      ),


    Tair_p95 =
      quantile(
        dat$Tair,
        0.95,
        na.rm = TRUE,
        names = FALSE
      ),


    Tair_IQR =
      quantile(
        dat$Tair,
        0.75,
        na.rm = TRUE,
        names = FALSE
      ) -
      quantile(
        dat$Tair,
        0.25,
        na.rm = TRUE,
        names = FALSE
      ),


    Tair_range90 =
      quantile(
        dat$Tair,
        0.95,
        na.rm = TRUE,
        names = FALSE
      ) -
      quantile(
        dat$Tair,
        0.05,
        na.rm = TRUE,
        names = FALSE
      ),


    mean_absolute_anomaly =
      mean(
        abs(T_anomaly),
        na.rm = TRUE
      ),


    # --------------------------------------------------------------------------------------
    # Fine-scale/local temperature heterogeneity
    # --------------------------------------------------------------------------------------

    mean_local_contrast =
      mean(
        abs(local_contrast),
        na.rm = TRUE
      ),


    median_local_contrast =
      median(
        abs(local_contrast),
        na.rm = TRUE
      ),


    p95_local_contrast =
      quantile(
        abs(local_contrast),
        0.95,
        na.rm = TRUE,
        names = FALSE
      )
  )
}


# ==========================================================================================
# FUNCTION: PROCESS ONE DATETIME
# ==========================================================================================

process_datetime <- function(
    dt,
    plotnr,
    voxel_TLS,
    scale_factor
) {


  # ----------------------------------------------------------------------------------------
  # Model setup
  # ----------------------------------------------------------------------------------------

  list2env(
    params,
    envir = .GlobalEnv
  )


  create_input_drivers()
  create_physical_constants()


  datetime <-
    as.POSIXct(
      dt,
      tz = "UTC"
    )


  # ----------------------------------------------------------------------------------------
  # Run model safely
  # ----------------------------------------------------------------------------------------

  out <- tryCatch({


    # ======================================================================================
    # IMPORT METEOROLOGICAL DRIVERS
    # ======================================================================================

    import_temp_observations(
      datetime
    )


    import_rad_observations(
      datetime
    )


    import_RH_u_observations(
      datetime
    )


    # ======================================================================================
    # PREPARE TLS GRID
    # ======================================================================================

    grid_local <-
      voxel_TLS$grid


    grid_local$density <-
      grid_local$density *
      scale_factor


    # ======================================================================================
    # RUN FOREDGECLIM
    # ======================================================================================

    res <-
      run_foredgeclim(
        grid_local,
        datetime
      )


    micro_grid <-
      res$micro_grid


    air_temp <-
      res$air_temperature


    if (
      any(
        !is.finite(
          air_temp
        )
      )
    ) {

      stop(
        "NaN/Inf in air temperature"
      )
    }


    # ======================================================================================
    # CREATE TEMPORARY 3D AIR TEMPERATURE GRID
    #
    # This exists only in memory for this simulation.
    # Individual voxel values are NOT saved.
    # ======================================================================================

    temp_grid <-
      micro_grid %>%

      transmute(

        x = x,

        y = y,

        z = z,

        Tair =
          air_temp -
          273.15
      )


    # ======================================================================================
    # CALCULATE HETEROGENEITY FOR EACH HORIZONTAL 1 M LAYER
    # ======================================================================================

    heterogeneity <-
      temp_grid %>%

      group_by(
        z
      ) %>%

      group_modify(
        ~ calculate_layer_heterogeneity(
          .x
        )
      ) %>%

      ungroup() %>%

      mutate(

        plot =
          plotnr,

        datetime =
          datetime

      ) %>%

      relocate(
        plot,
        datetime,
        z
      )


    # ======================================================================================
    # RETURN ONLY SUMMARIZED HETEROGENEITY
    # ======================================================================================

    heterogeneity


  }, error = function(e) {


    warning(
      paste(
        "Crash:",
        plotnr,
        format(datetime),
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


for (
  p in seq_len(
    nrow(plot_info)
  )
) {


  # ----------------------------------------------------------------------------------------
  # Plot information
  # ----------------------------------------------------------------------------------------

  info <-
    plot_info[p, ]


  plotnr <-
    info$plotnr


  height_canopy <<-
    info$canopy_height


  scale_factor <-
    info$scale_factor


  voxel_TLS <-
    readRDS(
      info$tls_file
    )


  cat(
    "\nProcessing plot:",
    plotnr,
    "\n"
  )


  # ----------------------------------------------------------------------------------------
  # Parallel run
  # ----------------------------------------------------------------------------------------

  with_progress({


    prog <-
      progressor(
        along =
          all_datetimes
      )


    parts <-
      future_lapply(

        seq_along(
          all_datetimes
        ),

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


  all_parts[[plotnr]] <-
    parts
}


# ==========================================================================================
# COMBINE OUTPUT
# ==========================================================================================

heterogeneity_all <-
  bind_rows(
    unlist(
      all_parts,
      recursive = FALSE
    )
  )


# ==========================================================================================
# ADD SOME USEFUL TIME VARIABLES
# ==========================================================================================

heterogeneity_all <-
  heterogeneity_all %>%

  mutate(

    date =
      as.Date(
        datetime,
        tz = "UTC"
      ),

    hour =
      as.integer(
        format(
          datetime,
          "%H"
        )
      ),

    doy =
      as.integer(
        format(
          datetime,
          "%j"
        )
      )
  )


# ==========================================================================================
# SAVE
# ==========================================================================================

dir.create(
  "Output/heterogeneity",
  recursive = TRUE,
  showWarnings = FALSE
)


saveRDS(

  heterogeneity_all,

  paste0(
    "Output/heterogeneity/heterogeneity_",
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


end_analysis <-
  Sys.time()


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
