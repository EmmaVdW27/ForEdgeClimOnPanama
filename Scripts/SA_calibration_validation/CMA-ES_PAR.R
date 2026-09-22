################################################################################
# ForEdgeClim shortwave-radiation calibration – CMA-ES
#
# Calibrates the shortwave-radiation parameters omega and Kb_v for one forest
# plot using PAR observations from four selected calibration days.
#
# For each candidate parameter set, the script:
#   - runs ForEdgeClim for daytime hours (11–22 UTC);
#   - extracts simulated shortwave radiation at 1 m height;
#   - matches simulated values to PAR sensor locations;
#   - aggregates observations and simulations to daily medians per sensor;
#   - evaluates temporal and spatial Pearson and Spearman correlations.
#
# CMA-ES minimises the objective:
#   4 - (temporal Pearson + temporal Spearman +
#        spatial Pearson + spatial Spearman)
#
# The optimisation therefore gives equal weight to temporal and spatial
# agreement in both absolute correlation and rank correlation.
#
# Author: Emma Van de Walle – Q-ForestLab
################################################################################

# ---- Packages ----
suppressPackageStartupMessages({
  library(ForEdgeClimOnPanama)
  library(cmaesr)
  library(smoof)
  library(ParamHelpers)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(purrr)
  library(tibble)
  library(future)
  library(future.apply)
  library(lubridate)
  library(RANN)
})

# Parallelization on HPC (supercomputer)
cores <- 32
plan(multicore, workers = cores)
print(paste("Using", cores, "cores."))

start_analysis <- Sys.time()

#########
# INPUT #
#########

plotnr = "C8" # R7 or C8

date_string <- 'C8_50_generations_2_params'
max_it        <- 49     # generations (there will be max_it + 1 generations)

# 4 selected calibration days, resp. most sunny and sun-variable day, most humid and diffuse light day, most windy day, and
# day with most day-night temperature difference
selected_days <- as.Date(c(
  "2023-03-28",
  "2023-11-02",
  "2023-12-20",
  "2023-06-07"
))

# For each day: hours 11–22 UTC
datetimes <- unlist(
  lapply(selected_days, function(d) {
    seq(
      from = as.POSIXct(paste0(d, " 11:00:00"), tz = "UTC"),
      to   = as.POSIXct(paste0(d, " 22:00:00"), tz = "UTC"),
      by   = "hour"
    )
  })
)

stop_fitness  <- 2      # corr target
lambda        <- NULL   # if NULL we set a heuristic later

#########

# Output paths
log_dir     <- file.path("Output", "calibration", "data")
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)
logfile_offspring   <- file.path(log_dir, paste0("logging_offspring_", date_string, ".csv"))
logfile_gener       <- file.path(log_dir, paste0("logging_generation_", date_string, ".csv"))
results_output_file <- file.path(log_dir, paste0("calibration_results_", date_string, ".rds"))
metrics_output_file <- file.path(log_dir, paste0("calibration_metrics_", date_string, ".rds"))

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
# PARAMETER NAMES & DEFAULTS
# -----------------
# Key parameters
param_names <- c(
  "omega",
  "Kb_v"
)
# Fixed parameters (SW RTM)
betad <<- 0.325; beta0 <<- 0.325; omega_g_v <<- 0.13; Kd_v <<- 0.775;
# Fixed parameters (LW RTM)
e_forest <<- 0.965; beta_lw <<- 0.325; omega_lw <<- 0.035; Kd_lw_v <<- 0.3; omega_g_lw_v <<- 0.055;
# Fixed parameters (OTHER)
h <<- 10; g_macro_ref <<- 25; infl_macro_ref <<- 32.5; infl_soil_ref <<- 5; infl_forest_ref <<- 5; g_forest_ref <<- 12.5;
p_ground <<- 0.225; g_soil_ref <<- 10; k_soil <<- 1.225; gs_max <<- 0.004; stomata_closure <<- 0.75; k_u <<- 0.7; u_ref <<- 1.6

# Initial values (mean value from uniform parameter distributions)
start_vals <- c(0.52, 1.25)


# Bounds
lower_bounds <- c(
  0.43,
  0.50
)

upper_bounds <- c(
  0.61,
  2.00
)


# -----------------
# MODEL INPUTS
# -----------------
create_input_drivers()

create_physical_constants()

# State for logging without monitor callbacks
.eval_counter <- 0L # variable accessible in multiple functions '.' is a convention: 'do not touch'

# Helper to set parameters
set_params_from_vec <- function(par) {

  omega <<- par[1]
  Kb_v  <<- par[2]

}


# -----------------
# PRECOMPUTE & CACHE TOMST OBSERVATIONS
# -----------------
cat("\nPrecomputing PAR observation files...\n")

obs_data <- list()

grid_base <- voxel_TLS$grid
grid_base$density <- grid_base$density * scale_factor

# use start_vals as temporary parameters
set_params_from_vec(start_vals)

for (dt in datetimes) {

  # convert only if needed
  if (inherits(dt, "POSIXct")) {
    datetime <- dt
  } else {
    datetime <- as.POSIXct(dt, origin = "1970-01-01", tz = "UTC")
  }

  key <- format(datetime, "%Y%m%d_%H")

  file <- paste0("Data/PAR_obs_", key, "_", plotnr, ".csv")

  # Only generate if not yet existing
  if (!file.exists(file)) {
    cat("Generating:", key, "\n")

    # imports
    import_temp_observations(datetime)
    import_rad_observations(datetime)
    import_RH_u_observations(datetime)

    # run model once (with temporary variables)
    run_foredgeclim(grid_base, datetime)
  }

  # read once into memory
  obs_data[[key]] <- data.table::fread(file)
}

cat("PAR data preloaded in memory.\n")

# evaluator that does everything
evaluate_par <- function(par) {

  parts <- future_lapply(
    datetimes,
    FUN = function(dt) {

      set_params_from_vec(par)
      # convert only if needed
      if (inherits(dt, "POSIXct")) {
        datetime <- dt
      } else {
        datetime <- as.POSIXct(dt, origin = "1970-01-01", tz = "UTC")
      }

      tryCatch({

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

        # ---- Observations en model outputs ----
        key <- format(datetime, "%Y%m%d_%H")
        PAR_obs <- obs_data[[key]]

        if(is.null(PAR_obs))
          stop("PAR file missing")

        sw_grid <- res$sw_rad

        sw_grid <- sw_grid[
          sw_grid$Z == 1,
        ]

        sw_grid$PAR_model <-
          (sw_grid$F_d_down +
             sw_grid$F_b_down) * 2.02 # no futher conversion needed since we check correlations and not absolute values

        grid_coords <-
          as.matrix(sw_grid[, c("X","Y")])

        sensor_coords <-
          as.matrix(PAR_obs[, c("x","y")])

        # nearest neighbour
        nn <- nn2(
          data = grid_coords, query = sensor_coords, k = 1 # amount of nearest neighbours
        )
        PAR_obs$PAR_model <- sw_grid$PAR_model[nn$nn.idx[,1]]

        valid <-
          is.finite(PAR_obs$radiation_mean) &
          is.finite(PAR_obs$PAR_model)

        if(!any(valid))
          stop("No valid PAR values")

        obs <- PAR_obs$radiation_mean[valid]

        sim <- PAR_obs$PAR_model[valid]

        sensor_ids <- PAR_obs$location_ID_number[valid]

        data.frame(
          date = as.Date(datetime),
          sensor_id = sensor_ids,
          obs = obs,
          sim = sim
        )
      }, error = function(e) {
        warning(paste("Crash during evaluation:", format(datetime), "-", e$message))
        NULL
      })
    },
    future.seed = TRUE
  )

  all_data <- bind_rows(parts)

  n_failed <- sum(vapply(parts, is.null, logical(1)))

  if(nrow(all_data) == 0){
    return(list(
      objective = 1e6,
      metrics = list(
        temporal_pearson = NA,
        temporal_spearman = NA,
        spatial_pearson = NA,
        spatial_spearman = NA
      ),
      n_total = 0,
      n_failed = length(parts)
    ))
  }

  daily_data <- all_data %>%
    group_by(date, sensor_id) %>%
    summarise(
      obs = median(obs, na.rm = TRUE),
      sim = median(sim, na.rm = TRUE),
      n_hours = n(),
      .groups = "drop"
    )

  # only keep days with min 8 hourly observations
  daily_data <- daily_data %>%
    filter(n_hours >= 8)

  if(nrow(daily_data) == 0){
    return(list(
      objective = 1e6,
      metrics = list(
        temporal_pearson = NA,
        temporal_spearman = NA,
        spatial_pearson = NA,
        spatial_spearman = NA
      ),
      n_total = 0,
      n_failed = length(parts)
    ))
  }

  temporal_metrics <- daily_data %>%
    group_by(sensor_id) %>%
    summarise(

      temporal_pearson =
        ifelse(
          n() > 2,
          cor(
            obs,
            sim,
            method = "pearson",
            use = "complete.obs"
          ),
          NA_real_
        ),

      temporal_spearman =
        ifelse(
          n() > 2,
          cor(
            obs,
            sim,
            method = "spearman",
            use = "complete.obs"
          ),
          NA_real_
        ),

      .groups = "drop"
    )

  temporal_pearson <-
    median(
      temporal_metrics$temporal_pearson,
      na.rm = TRUE
    )

  temporal_spearman <-
    median(
      temporal_metrics$temporal_spearman,
      na.rm = TRUE
    )

  spatial_metrics <- daily_data %>%
    group_by(date) %>%
    summarise(

      spatial_pearson =
        ifelse(
          n() > 2,
          cor(
            obs,
            sim,
            method = "pearson",
            use = "complete.obs"
          ),
          NA_real_
        ),

      spatial_spearman =
        ifelse(
          n() > 2,
          cor(
            obs,
            sim,
            method = "spearman",
            use = "complete.obs"
          ),
          NA_real_
        ),

      .groups = "drop"
    )

  spatial_pearson <-
    median(
      spatial_metrics$spatial_pearson,
      na.rm = TRUE
    )

  spatial_spearman <-
    median(
      spatial_metrics$spatial_spearman,
      na.rm = TRUE
    )
  objective <-
    4 -
    (
      temporal_pearson +
        temporal_spearman +
        spatial_pearson +
        spatial_spearman
    )

  list(
    objective = objective,

    metrics = list(

      temporal_pearson =
        temporal_pearson,

      temporal_spearman =
        temporal_spearman,

      spatial_pearson =
        spatial_pearson,

      spatial_spearman =
        spatial_spearman

    ),

    n_total = nrow(daily_data),

    n_failed = n_failed
  )
}

# -----------------
# COST = OBJECTIVE FUNCTION
# -----------------
compute_objective <- function(par) {
  evaluate_par(par)$objective
}

# -----------------
# EXTRA METRICS
# -----------------
compute_metrics <- function(par) {
  evaluate_par(par)$metrics
}

# -------------------
# Objective wrapper to log each offspring evaluation
# -------------------
# Wrapper around objective function so logging is performed
wrapped_obj_fn <- function(x) {
  objective <- compute_objective(x)
  .eval_counter <<- .eval_counter + 1L
  setTxtProgressBar(pb, .eval_counter)
  # infer generation and offspring id from eval count and lambda
  gen <- 1L + floor((.eval_counter - 1L) / lambda)
  off <- 1L + ((.eval_counter - 1L) %% lambda)
  # build row
  row <- tibble::tibble(
    eval_id = .eval_counter,
    iter = gen,
    offspring_id = off,
    objective = objective,
    time = Sys.time()
  ) %>% dplyr::bind_cols(tibble::as_tibble_row(stats::setNames(as.list(as.numeric(x)), param_names)))
  offspring_log <<- dplyr::bind_rows(offspring_log, row)
  objective
}

# -----------------
# SMOOF WRAPPER with bounds
# -----------------
# Specific Smoof Wrapper around wrapped ojective function
# Necessary to set the objective function (compute_rmse) into a readable format for cmaesr, ie, cmaesr needs a smoof object
dim_par <- length(start_vals)
obj_logged <- smoof::makeSingleObjectiveFunction(
  name = "ForEdgeClim_correlation_logged",
  fn   = function(x) wrapped_obj_fn(x), # function to be optimized
  has.simple.signature = TRUE, # fn accepts a simple vector as input, not a list
  par.set = ParamHelpers::makeNumericParamSet(len = dim_par, lower = lower_bounds, upper = upper_bounds), # number of params and search space (lower/upper bounds)
  minimize = TRUE # minimise the correlation-based objective
)


# -----------------
# DEFINING LOGGING STRUCTURES
# -----------------
# offspring log: one row per objective evaluation (i.e. per offspring)
offspring_log <- tibble::tibble(
  eval_id = integer(),
  iter    = integer(),           # generation index (1-based)
  offspring_id = integer(),      # lambda within generation
  objective    = numeric(),
  time    = as.POSIXct(NA, tz = "UTC")
) %>%
  dplyr::bind_cols(tibble::as_tibble_row(stats::setNames(as.list(rep(NA_real_, dim_par)), param_names))) %>%
  dplyr::slice(0) # no data yet

# generation log (to be derived post-run)
generation_log <- tibble::tibble(
  iter = integer(),
  best_value = numeric(),
  sigma_proxy = numeric(),
  spread_max = numeric(),
  spread_mean = numeric(),
  mean_step = numeric()
)

# -----------------
# RUN CMA-ES
# -----------------
# Choose lambda if not provided
if (is.null(lambda)) lambda <- 4 + floor(3 * log(dim_par))

# progress bar settings
total_evals <- lambda * (max_it + 1)
pb <- txtProgressBar(min = 0, max = total_evals, style = 3)

# initial sigma heuristic (initial global step size)
init_sigma <- 0.3 * mean(upper_bounds - lower_bounds)

# Build a robust stop list
# Stop after max_it iterations
stop_list <- list(cmaesr::stopOnMaxIters(max_it))
# Add default stopping conditions if available in your cmaesr version
# These include e.g. too small step sizes, no variance in population any more...
if ("getDefaultStoppingConditions" %in% getNamespaceExports("cmaesr")) {
  stop_list <- c(stop_list, cmaesr::getDefaultStoppingConditions())
}
# Target fitness as additional criterion if supported
if ("stopOnOptFitness" %in% getNamespaceExports("cmaesr")) {
  stop_list <- c(stop_list, list(cmaesr::stopOnOptFitness(stop_fitness)))
}

# make results reproducible
set.seed(1)
res <- cmaesr::cmaes(
  obj_logged, # Smoof object with logging
  monitor = cmaesr::makeSimpleMonitor(), # simple progressing output in console
  control = list(
    lambda   = lambda,
    sigma    = init_sigma,
    stop.ons = stop_list
  )
)

saveRDS(res, file = results_output_file)

close(pb)

# -----------------
# DERIVE GENERATION-LEVEL LOG FROM OFFSPRING LOG
# -----------------

# Function to get spread values to investigate convergence
# M is a matrix of dimension lambda x d (number of params) with the offspring params for several iterations.
# Convergence would lead to small sigma_proxy, spread_max and spread_mean values.
get_spreads <- function(M) {
  sds <- apply(M, 2, stats::sd) # get sd for each param (= each column)
  # Covariance matrix representing the multivariate spread (including correlations between params)
  # The eigenvalues of this matrix are the variances along the main axes (principal components, PC).
  # The diagonal of a covariance matrix contains the variance of each parameter; the off-diagonal elements contain covariances between parameters.
  # CMA-ES learns correlations between parameters: the algorithm eg detects improvements when changing combinations of parameters,
  # then the matrix is adapted, meaning the search directions (= PCs) are updated. Main search is in the direction of the main PC,
  # since it has the highest eigenvalue.
  # There are as many searching directions as parameters.
  covM <- stats::cov(M)
  ev <- tryCatch(eigen(covM, symmetric = TRUE)$values, error = function(e) rep(NA_real_, ncol(M)))
  list(
    sigma_proxy = mean(sds, na.rm = TRUE), # mean of param sds is a proxy for the global step size (= distribution, spread)
    spread_max  = sqrt(max(ev, na.rm = TRUE)), # largest deviation (sd) along main PC = spread in main searching direction
    spread_mean = sqrt(mean(ev, na.rm = TRUE)) # mean deviation in the PC-space = mean spread over all PCs/searching directions = general value for spread
  )
}

param_mat_cols <- param_names

generation_log <- offspring_log %>%
  dplyr::group_by(iter) %>%
  dplyr::summarise(
    best_value = min(objective, na.rm = TRUE), # # get lowest objective value within each generation
    params = list(as.matrix(dplyr::pick(dplyr::all_of(param_mat_cols)))), # make a list with for every generation a lambda x d matrix with paramvalues
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    spreads = purrr::map(params, get_spreads), # get spread values per generation
    sigma_proxy = purrr::map_dbl(spreads, "sigma_proxy"), # get sigma_proxy from spreads
    spread_max  = purrr::map_dbl(spreads, "spread_max"),
    spread_mean = purrr::map_dbl(spreads, "spread_mean"),
    mean_param  = purrr::map(params, ~ colMeans(.x, na.rm = TRUE)) # mean param value per generation
  ) %>%
  # Get change in population mean between generations
  # This measures how strong the population center shifts from one generation to the next.
  # This change is defined as the Euclidean distance between consecutive generation means.
  # A low value means convergence.
  dplyr::mutate(mean_step = c(NA_real_,
                              purrr::map2_dbl(mean_param[-length(mean_param)], mean_param[-1],
                                              ~ sqrt(sum((.y - .x)^2))))) %>%
  # Keep columns of interest
  dplyr::select(iter, best_value, sigma_proxy, spread_max, spread_mean, mean_step)

# -----------------
# SAVE LOGS
# -----------------
readr::write_csv(offspring_log, logfile_offspring)
readr::write_csv(generation_log, logfile_gener)

# ---- Search best offspring ----
best_offspring <- offspring_log %>%
  dplyr::filter(is.finite(objective)) %>%
  dplyr::slice_min(objective, n = 1)

# Get corresponding parameter-vector
best_offspring_par <- unlist(
  best_offspring %>%
    dplyr::select(dplyr::all_of(param_names)),
  use.names = FALSE
)

# Console summaries
cat("\nOptimised parameters (final):\n")
print(res$best.param)
cat("\nBest objective value:\n")
print(res$best.fitness)
cat("\nBest offspring based on objective in offspring_log (should be same as above):\n")
print(best_offspring$objective)

# -----------------
# EXTRA METRICS FOR BEST PARAMETERS SET
# -----------------
best_metrics <- compute_metrics(best_offspring_par)
plan(sequential)


cat("\nPerformance metrics for best parameter set:\n")
print(best_metrics)

saveRDS(best_metrics, file = metrics_output_file)


end_analysis <- Sys.time()
cat(sprintf('\nTotal running time = %.2f h\n', as.numeric(end_analysis - start_analysis, units = 'hours')))
