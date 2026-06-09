############################################################################################
# This script optimizes ForEdgeClim parameters for the whole year with the CMA-ES algorithm (via cmaesr).
# In the objective function, we minimize RMSE (Root Mean Square Error).
#
# Author: Emma Van de Walle - Q-ForestLab
############################################################################################

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
cores <- as.numeric(Sys.getenv("PBS_NUM_PPN", unset = 1))
plan(multicore, workers = cores)
print(paste("Using", availableCores(), "cores."))

start_analysis <- Sys.time()

#########
# INPUT #
#########

plotnr = "C13" # C13 or R11

# We calibrate on 12 days - selected over the 12 months, capturing both seasons

date_string <- 'C13_year_30_generations'
max_it        <- 29     # generations (there will be max_it + 1 generations)

# 4 selected calibration days, resp. most sunny and sun-variable day, most humid and diffuse light day, most windy day, and
# day with most day-night temperature difference
selected_days <- as.Date(c(
  "2023-03-28",
  "2023-11-12",
  "2023-12-21",
  "2023-06-05"
))

# For each day: 24 hours (0–23h)
datetimes <- unlist(
  lapply(selected_days, function(d) {
    seq(
      from = as.POSIXct(paste0(d, " 00:00:00"), tz = "UTC"),
      to   = as.POSIXct(paste0(d, " 23:00:00"), tz = "UTC"),
      by   = "hour"
    )
  })
)

stop_fitness  <- 1      # RMSE target (°C)
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
  height_canopy = 33
  scale_factor = 1
  voxel_TLS <- readRDS("Data/TLS_scaled_DTM_and_grid_C13_Panama.rds")
}else{
  height_canopy = 35
  scale_factor = 6.20/6.22
  voxel_TLS <- readRDS("Data/TLS_scaled_DTM_and_grid_R11_Panama.rds")
}

# -----------------
# PARAMETER NAMES & DEFAULTS
# -----------------
# Key parameters
param_names <- c("infl_macro_ref", "infl_soil_ref", "g_macro_ref", "e_forest", "u_ref", "k_soil")
# Fixed parameters (SW RTM)
betad <<- 0.325; beta0 <<- 0.325; omega <<- 0.52; Kd_v <<- 0.775; Kb_v <<- 1.25; omega_g_v <<- 0.13;
# Fixed parameters (LW RTM)
beta_lw <<- 0.325; omega_lw <<- 0.035; Kd_lw_v <<- 0.3; omega_g_lw_v <<- 0.055;
# Fixed parameters (HEAT)
h <<- 10; infl_forest_ref <<- 5; g_forest_ref <<- 12.5; p_ground <<- 0.225; g_soil_ref <<- 10;
gs_max <<- 0.004; stomata_closure <<- 0.75; k_u <<- 0.7

# Initial values (mean value from uniform parameter distributions)
start_vals <- c(32.5, 5, 25, 0.965, 1.6, 1.225)

# Bounds
lower_bounds <- c(5, 0, 10, 0.94, 0.2, 0.25)
upper_bounds <- c(60, 10, 40, 0.99, 3, 2.2)


# -----------------
# MODEL INPUTS
# -----------------
create_input_drivers()

create_physical_constants()

# State for logging without monitor callbacks
.eval_counter <- 0L # variable accessible in multiple functions '.' is a convention: 'do not touch'

# Helper to set parameters
set_params_from_vec <- function(par) {
  infl_macro_ref <<- par[1]; infl_soil_ref <<- par[2]; g_macro_ref <<- par[3]; e_forest <<- par[4]; u_ref <<- par[5]; k_soil <<- par[6]
}


# -----------------
# PRECOMPUTE & CACHE TOMST OBSERVATIONS
# -----------------
cat("\nPrecomputing TOMST observation files...\n")

obs_data <- list()

grid_base <- voxel_TLS$grid
grid_base$density <- grid_base$density * scale_factor

# start_vals as temporary vals
set_params_from_vec(start_vals)

for (dt in datetimes) {

  datetime <- as.POSIXct(dt, tz = "UTC")
  key <- format(datetime, "%Y%m%d_%H")

  file <- paste0("Data/microtemp_air_obs_", key, "_", plotnr, ".csv")

  # Only generate if not yet existing
  if (!file.exists(file)) {
    cat("Generating:", key, "\n")

    # imports
    import_temp_observations(datetime)
    import_rad_observations(datetime)
    import_RH_u_observations(datetime)

    # run model once (met tijdelijke parameters)
    run_foredgeclim(grid_base, datetime)
  }

  # read once into memory
  obs_data[[key]] <- data.table::fread(file)
}

cat("TOMST data preloaded in memory.\n")

# evaluator that does everything (RMSE, R2, ME, NSE and SD)
evaluate_par <- function(par) {

  parts <- future_lapply(
    datetimes,
    FUN = function(dt) {

      set_params_from_vec(par)
      datetime <- as.POSIXct(dt, tz = "UTC")

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

        # ---- Observations ----
        key <- format(datetime, "%Y%m%d_%H")
        TOMST <- obs_data[[key]]
        #TOMST <- read.csv(paste0("Data/microtemp_air_obs_", key, "_", plotnr, ".csv"))

        # ---- Model extraction ----
        temp_air_grid <- micro_grid
        temp_air_grid$temperature <- air_temp - 273.15
        temp_air_grid <- temp_air_grid[temp_air_grid$z == 1, c("x","y","temperature")]

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

        sim <- TOMST$temperature_model[valid]
        obs <- TOMST$value_mean[valid]

        list(
          ok          = TRUE,
          sse         = sum((sim - obs)^2),
          n           = length(sim),
          sum_sim     = sum(sim),
          sum_obs     = sum(obs),
          sum_sim2    = sum(sim^2),
          sum_obs2    = sum(obs^2),
          sum_sim_obs = sum(sim * obs),
          sum_diff    = sum(diffs[valid]),
          sum_diff2   = sum(diffs[valid]^2)
        )

      }, error = function(e) {
        warning(paste("Crash during evaluation:", format(datetime), "-", e$message))
        list(
          ok = FALSE,
          sse = NA_real_, n = 0L,
          sum_sim = 0, sum_obs = 0,
          sum_sim2 = 0, sum_obs2 = 0,
          sum_sim_obs = 0,
          sum_diff = 0,
          sum_diff2 = 0
        )
      })
    },
    future.seed = TRUE
  )

  # ---- Aggregate ----
  sse_total         <- sum(vapply(parts, `[[`, numeric(1), "sse"),         na.rm = TRUE)
  n_total           <- sum(vapply(parts, `[[`, integer(1), "n"))
  sum_sim_total     <- sum(vapply(parts, `[[`, numeric(1), "sum_sim"),     na.rm = TRUE)
  sum_obs_total     <- sum(vapply(parts, `[[`, numeric(1), "sum_obs"),     na.rm = TRUE)
  sum_sim2_total    <- sum(vapply(parts, `[[`, numeric(1), "sum_sim2"),    na.rm = TRUE)
  sum_obs2_total    <- sum(vapply(parts, `[[`, numeric(1), "sum_obs2"),    na.rm = TRUE)
  sum_sim_obs_total <- sum(vapply(parts, `[[`, numeric(1), "sum_sim_obs"), na.rm = TRUE)
  sum_diff_total  <- sum(vapply(parts, `[[`, numeric(1), "sum_diff"),  na.rm = TRUE)
  sum_diff2_total <- sum(vapply(parts, `[[`, numeric(1), "sum_diff2"), na.rm = TRUE)
  n_failed          <- sum(!vapply(parts, `[[`, logical(1), "ok"))

  if (n_total == 0) {
    return(list(
      RMSE = 1e6,
      metrics = list(RMSE = NA_real_, R2 = NA_real_, NSE = NA_real_, ME = NA_real_),
      n_total = 0L,
      n_failed = n_failed
    ))
  }

  rmse <- sqrt(sse_total / n_total)

  mean_sim <- sum_sim_total / n_total
  mean_obs <- sum_obs_total / n_total

  denom_nse <- sum_obs2_total - n_total * mean_obs^2
  nse <- if (denom_nse > 0) 1 - sse_total / denom_nse else NA_real_

  cov_num     <- sum_sim_obs_total - n_total * mean_sim * mean_obs
  var_sim_num <- sum_sim2_total    - n_total * mean_sim^2
  var_obs_num <- sum_obs2_total    - n_total * mean_obs^2

  r2 <- if (var_sim_num > 0 && var_obs_num > 0) (cov_num^2) / (var_sim_num * var_obs_num) else NA_real_

  me <- (sum_sim_total - sum_obs_total) / n_total

  # SD of residuals (sample SD)
  sd_res <- if (n_total > 1) {
    sqrt((sum_diff2_total - (sum_diff_total^2) / n_total) / (n_total - 1))
  } else {
    NA_real_
  }


  list(
    RMSE = rmse,
    metrics = list(RMSE = rmse, R2 = r2, NSE = nse, ME = me, SD = sd_res),
    n_total = n_total,
    n_failed = n_failed
  )
}

# -----------------
# COST = OBJECTIVE FUNCTION (RMSE)
# -----------------
compute_rmse <- function(par) {
  evaluate_par(par)$RMSE
}


# -----------------
# EXTRA METRICS: RMSE, R2, NSE, ME, SD
# -----------------
compute_metrics <- function(par) {
  evaluate_par(par)$metrics
}

# -------------------
# Objective wrapper to log each offspring evaluation
# -------------------
# Wrapper around objective function so logging is performed
wrapped_obj_fn <- function(x) {
  rmse <- compute_rmse(x)
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
    rmse = rmse,
    time = Sys.time()
  ) %>% dplyr::bind_cols(tibble::as_tibble_row(stats::setNames(as.list(as.numeric(x)), param_names)))
  offspring_log <<- dplyr::bind_rows(offspring_log, row)
  rmse
}

# -----------------
# SMOOF WRAPPER with bounds
# -----------------
# Specific Smoof Wrapper around wrapped ojective function
# Necessary to set the objective function (compute_rmse) into a readable format for cmaesr, ie, cmaesr needs a smoof object
dim_par <- length(start_vals)
obj_logged <- smoof::makeSingleObjectiveFunction(
  name = "ForEdgeClim_RMSE_logged",
  fn   = function(x) wrapped_obj_fn(x), # function to be optimized
  has.simple.signature = TRUE, # fn accepts a simple vector as input, not a list
  par.set = ParamHelpers::makeNumericParamSet(len = dim_par, lower = lower_bounds, upper = upper_bounds), # number of params and search space (lower/upper bounds)
  minimize = TRUE # minimize RMSE, direction of optimalization
)


# -----------------
# DEFINING LOGGING STRUCTURES
# -----------------
# offspring log: one row per objective evaluation (i.e. per offspring)
offspring_log <- tibble::tibble(
  eval_id = integer(),
  iter    = integer(),           # generation index (1-based)
  offspring_id = integer(),      # lambda within generation
  rmse    = numeric(),
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
# Stop if best RMSE =< stop_fitness
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
  # The diagonal of a covariance matrix contains the variance of 1 param, the other elements the correlations between 2 params.
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
    best_value = min(rmse, na.rm = TRUE), # get lowest RMSE within each generation
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

# ---- Search best offspring based on RMSE ----
best_offspring <- offspring_log %>%
  dplyr::filter(is.finite(rmse)) %>%
  dplyr::slice_min(rmse, n = 1)   # row with lowest RMSE

# Get corresponding parameter-vector
best_offspring_par <- best_offspring %>%
  dplyr::select(dplyr::all_of(param_names)) %>%
  as.numeric()

# Console summaries
cat("\nOptimised parameters (final):\n")
print(res$best.param)
cat("\nMinimal RMSE observed:\n")
print(res$best.fitness)
cat("\nBest offspring based on RMSE in offspring_log (should be same as above):\n")
print(best_offspring$rmse)

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
