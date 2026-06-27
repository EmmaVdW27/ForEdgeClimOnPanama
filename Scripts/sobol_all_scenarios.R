
###############################################################################
# Sobol sensitivity analysis for ForEdgeClim
# Treatments:
#   control = C8 + C13
#   removal = R7 + R11
#
# Seasons:
#   dry, wet
#
# Periods:
#   noon, night
#
# QoI:
#   1. mean_understorey_T
#   2. sd_understorey_T
#   3. vertical_T_gradient
#   4. mean_understorey_SW
#   5. sd_understorey_SW
#   6. vertical_SW_extinction
###############################################################################

library(ForEdgeClimOnPanama)
library(dplyr)
library(lhs)
library(future.apply)
library(future)
library(progressr)
library(purrr)
library(sensitivity)

handlers("progress")

cores <- as.numeric(Sys.getenv("PBS_NUM_PPN", unset = 1))
plan(multicore, workers = cores)

n <- 400
n_boot <- 1000

output_root <- "Output/sensitivity_analysis_treatments"

###############################################################################
# PARAMETER RANGES
###############################################################################

param_ranges <- data.frame(
  parameter = c(
    "betad","beta0","omega","Kd_v","Kb_v","omega_g_v",
    "e_forest","beta_lw","omega_lw","Kd_lw_v","omega_g_lw_v",
    "h","g_macro_ref","infl_macro_ref","infl_soil_ref",
    "infl_forest_ref","g_forest_ref","p_ground",
    "g_soil_ref","k_soil","gs_max","stomata_closure",
    "k_u","u_ref"
  ),
  min = c(
    0.3,0.2,0.43,0.6,0.5,0.08,
    0.94,0.3,0.01,0.2,0.04,
    0,10,5,0,0,5,0.1,5,0.25,
    0.002,0.3,0.2,0.2
  ),
  max = c(
    0.35,0.45,0.61,0.95,2,0.18,
    0.99,0.35,0.06,0.4,0.07,
    20,40,60,10,10,20,0.35,15,2.2,
    0.006,1.2,1.2,3
  )
)

param_names <- param_ranges$parameter
k <- nrow(param_ranges)

###############################################################################
# TREATMENTS
###############################################################################

get_treatment_settings <- function(treatment){

  if(treatment == "control"){

    return(list(

      structures = c(
        "Data/TLS_scaled_DTM_and_grid_C8_Panama.rds",
        "Data/TLS_scaled_DTM_and_grid_C13_Panama.rds"
      ),

      plots = c(
        "C8",
        "C13"
      ),

      heights = c(
        38,
        33
      ),

      scaling = c(
        5.81/6.22,
        1
      )

    ))
  }

  if(treatment == "removal"){

    return(list(

      structures = c(
        "Data/TLS_scaled_DTM_and_grid_R7_Panama.rds",
        "Data/TLS_scaled_DTM_and_grid_R11_Panama.rds"
      ),

      plots = c(
        "R7",
        "R11"
      ),

      heights = c(
        37,
        35
      ),

      scaling = c(
        6.19/6.22,
        6.20/6.22
      )

    ))
  }
}

###############################################################################
# DATETIMES
###############################################################################

get_datetime <- function(season, period){

  if(season == "dry" & period == "noon")
    return(as.POSIXct("2023-04-26 19:00:00", tz="UTC"))

  if(season == "dry" & period == "night")
    return(as.POSIXct("2023-04-26 07:00:00", tz="UTC"))

  if(season == "wet" & period == "noon")
    return(as.POSIXct("2023-11-26 19:00:00", tz="UTC"))

  if(season == "wet" & period == "night")
    return(as.POSIXct("2023-11-26 07:00:00", tz="UTC"))
}

###############################################################################
# METRICS
###############################################################################

compute_metrics <- function(res){

  understorey <- res$micro_grid$z == 1

  Ta_under <- res$air_temperature[understorey] - 273.15
  SW_under <- res$sw_rad$net_sw[understorey]

  mean_understorey_T <- mean(Ta_under, na.rm = TRUE)
  sd_understorey_T   <- sd(Ta_under, na.rm = TRUE)

  temp_profile <- data.frame(
    z = res$micro_grid$z,
    T = res$air_temperature - 273.15
  ) |>
    group_by(z) |>
    summarise(T = mean(T, na.rm = TRUE),
              .groups = "drop")

  vertical_T_gradient <-
    coef(lm(T ~ z, data = temp_profile))[2]

  mean_understorey_SW <- mean(SW_under, na.rm = TRUE)
  sd_understorey_SW   <- sd(SW_under, na.rm = TRUE)

  sw_profile <- data.frame(
    z = res$micro_grid$z,
    SW = pmax(res$sw_rad$net_sw, 0.001)
  ) |>
    group_by(z) |>
    summarise(SW = mean(SW, na.rm = TRUE),
              .groups = "drop")

  vertical_SW_extinction <-
    coef(lm(log(SW) ~ z, data = sw_profile))[2]

  c(
    mean_understorey_T = mean_understorey_T,
    sd_understorey_T = sd_understorey_T,
    vertical_T_gradient = vertical_T_gradient,
    mean_understorey_SW = mean_understorey_SW,
    sd_understorey_SW = sd_understorey_SW,
    vertical_SW_extinction = vertical_SW_extinction
  )
}

###############################################################################
# MODEL FUNCTION
###############################################################################

model_function <- function(
    param_values,
    datetime,
    treatment_settings
){

  create_input_drivers()
  create_physical_constants()
  create_model_parameters()

  for(i in seq_along(param_names)){
    assign(
      param_names[i],
      param_values[i],
      envir = .GlobalEnv
    )
  }

  metrics_all <- list()

  for(p in seq_along(treatment_settings$structures)){

    assign(
      "plotnr",
      treatment_settings$plots[p],
      envir = .GlobalEnv
    )

    assign(
      "height_canopy",
      treatment_settings$heights[p],
      envir = .GlobalEnv
    )

    import_temp_observations(datetime)
    import_rad_observations(datetime)
    import_RH_u_observations(datetime)

    voxel_TLS <- readRDS(
      treatment_settings$structures[p]
    )

    voxel_TLS$grid$density <-
      voxel_TLS$grid$density *
      treatment_settings$scaling[p]

    res <- run_foredgeclim(
      voxel_TLS$grid,
      datetime
    )

    metrics_all[[p]] <-
      compute_metrics(res)
  }

  metrics_matrix <- do.call(
    rbind,
    metrics_all
  )

  colMeans(metrics_matrix)
}

###############################################################################
# SOBOL DESIGN
###############################################################################

set.seed(42)

X1 <- randomLHS(n, k)
X2 <- randomLHS(n, k)

colnames(X1) <- param_names
colnames(X2) <- param_names

sobol_design <- soboljansen(
  model = NULL,
  X1 = X1,
  X2 = X2,
  nboot = n_boot,
  conf = 0.95
)

scale_sobol <- function(X){

  for(i in seq_len(k)){

    X[,i] <-
      param_ranges$min[i] +
      X[,i] *
      (param_ranges$max[i] - param_ranges$min[i])
  }

  colnames(X) <- param_names

  as.data.frame(X)
}

scaled_input <- scale_sobol(sobol_design$X)

###############################################################################
# SCENARIOS
###############################################################################

scenarios <- expand.grid(
  treatment = c("control","removal"),
  season = c("dry","wet"),
  period = c("noon","night"),
  stringsAsFactors = FALSE
)

summary_table <- list()

###############################################################################
# RUN SOBOL
###############################################################################

for(s in seq_len(nrow(scenarios))){

  treatment <- scenarios$treatment[s]
  season <- scenarios$season[s]
  period <- scenarios$period[s]

  scenario_name <- paste(
    treatment,
    season,
    period,
    sep = "_"
  )

  cat("\nRunning:", scenario_name, "\n")

  datetime <- get_datetime(
    season,
    period
  )

  settings <- get_treatment_settings(
    treatment
  )

  outdir <- file.path(
    output_root,
    scenario_name
  )

  dir.create(
    outdir,
    recursive = TRUE,
    showWarnings = FALSE
  )

  sobol_outputs_all <- with_progress({

    p <- progressor(
      steps = nrow(scaled_input)
    )

    future_lapply(

      seq_len(nrow(scaled_input)),

      function(i){

        metrics <- model_function(

          param_values =
            unlist(scaled_input[i,]),

          datetime = datetime,

          treatment_settings =
            settings
        )

        p()

        metrics
      },

      future.seed = TRUE
    )
  })

  metric_names <- names(
    sobol_outputs_all[[1]]
  )

  temperature_metrics <- c(
    "mean_understorey_T",
    "sd_understorey_T",
    "vertical_T_gradient"
  )

  sw_metrics <- c(
    "mean_understorey_SW",
    "sd_understorey_SW",
    "vertical_SW_extinction"
  )

  if(period == "night"){

    metric_names <- intersect(
      metric_names,
      temperature_metrics
    )

  }else{

    metric_names <- c(
      temperature_metrics,
      sw_metrics
    )

  }

  for(metric in metric_names){

    y <- map_dbl(
      sobol_outputs_all,
      metric
    )

    sobol_res <- tell(
      sobol_design,
      y
    )

    saveRDS(
      sobol_res,
      file.path(
        outdir,
        paste0(metric, ".rds")
      )
    )

    tmp <- data.frame(
      scenario = scenario_name,
      metric = metric,
      parameter = rownames(sobol_res$S),
      S1 = sobol_res$S$original,
      ST = sobol_res$T$original
    )

    summary_table[[length(summary_table)+1]] <- tmp
  }
}

###############################################################################
# SAVE SUMMARY
###############################################################################

summary_table <- bind_rows(
  summary_table
)

saveRDS(
  summary_table,
  file.path(
    output_root,
    "sobol_summary_table.rds"
  )
)

write.csv(
  summary_table,
  file.path(
    output_root,
    "sobol_summary_table.csv"
  ),
  row.names = FALSE
)

plan(sequential)

cat("\nFinished.\n")
