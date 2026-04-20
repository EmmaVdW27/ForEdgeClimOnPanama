###############################################################################
# This is a script to perform a Sobol sensitivity analysis on the ForEdgeClim
# model - version without edge effects.
# Sobol indices per quantity of interest are calculated and saved as well.
#
# Author: Emma Van de Walle - Q-ForestLab
###############################################################################

library(ForEdgeClimOnPanama)
library(dplyr)
library(ggplot2)
library(scales)
library(lhs)
library(future.apply)
library(future)
library(tidyr)
library(forcats)
library(purrr)
library(ggridges)
library(sensitivity)
library(ggtext)
library(progressr)

#########
# INPUT #
#########

plotnr = "C13" # C13 or R11

# dry season noon: 2023-04-26 19:00:00
# dry season night: 2023-04-26 07:00:00
# wet season noon: 2023-11-26 19:00:00
# wet season night: 2023-11-26 07:00:00
datehour = "2023-04-26 19:00:00" # in UTC not LC!

#########

# Parallelization on HPC (supercomputer)
handlers("progress")
cores <- as.numeric(Sys.getenv("PBS_NUM_PPN", unset = 1))
plan(multicore, workers = cores)
print(paste("Using", availableCores(), "cores."))

if(plotnr == "C13"){
  height_canopy = 33
  structure_PAI_scaling_factor = 1
}else{
  height_canopy = 35
  structure_PAI_scaling_factor = 6.20/6.22
}


datetime <- as.POSIXct(datehour, tz = "UTC")
structure <- paste0("Data/TLS_scaled_DTM_and_grid_", plotnr, "_Panama.rds")
n <- 400
n_boot <- 1000
output_path <- "Output/sensitivity_analysis/"

# Parameter ranges
param_ranges <- data.frame(
  parameter = c(
    "betad", "beta0", "omega", "Kd_v", "Kb_v", "omega_g_v",
    "e_forest", "beta_lw", "omega_lw", "Kd_lw_v", "omega_g_lw_v",
    "h", "g_macro_ref", "infl_macro_ref", "infl_soil_ref", "infl_forest_ref", "g_forest_ref", "p_ground", "g_soil_ref", "k_soil",
    "gs_max", "stomata_closure", "wind_factor"
  ),
  min = c(
    0.3, 0.2, 0.43, 0.6, 0.5, 0.08,
    0.94, 0.3, 0.01, 0.2, 0.04,
    0, 10, 5, 0, 0, 5, 0.1, 5, 0.25,
    0.2, 0.3, 0.5
  ),
  max = c(
    0.35, 0.45, 0.61, 0.95, 2, 0.18,
    0.99, 0.35, 0.06, 0.4, 0.07,
    20, 40, 60, 10, 10, 20, 0.35, 15, 2.2,
    0.8, 1.2, 10
  ),
  stringsAsFactors = FALSE
)

param_names <- param_ranges$parameter
k <- nrow(param_ranges)

formatted_datetime <- format(datetime, "%Hh_%d%m%Y")

output_name <- paste0(n, "samples_", k, "parameters_", formatted_datetime)

##################
# SOBOL SAMPLING #
##################

# Sobol design
set.seed(42)
X1 <- randomLHS(n, k)
X2 <- randomLHS(n, k)
colnames(X1) <- param_names
colnames(X2) <- param_names
sobol_design <- soboljansen(NULL, X1 = X1, X2 = X2, nboot = n_boot, conf = 0.95)

# Scaling
scale_sobol <- function(X) {
  for (i in seq_len(k)) {
    X[, i] <- param_ranges$min[i] + X[, i] * (param_ranges$max[i] - param_ranges$min[i])
  }
  colnames(X) <- param_ranges$parameter
  as.data.frame(X)
}
scaled_input <- scale_sobol(sobol_design$X)

#############################################
# MODEL DEFINITION WITH METRIC CALCULATIONS #
#############################################

model_function <- function(param_values) {
  create_input_drivers()
  create_physical_constants()
  create_model_parameters()
  for (i in seq_along(param_names)) {
    assign(param_names[i], param_values[i], envir = .GlobalEnv)
  }

  import_temp_observations(datetime)
  import_rad_observations(datetime)
  import_RH_u_observations(datetime)

  voxel_TLS <- readRDS(structure)
  voxel_TLS$grid$density = voxel_TLS$grid$density * structure_PAI_scaling_factor
  res <- run_foredgeclim(voxel_TLS$grid, datetime)

  # mean and sd on whole grid temperatures
  ta <- res$air_temperature - 273.15 # air temperature
  tf <- res$micro_grid$temperature - 273.15 # forest surface temperature

  ## air temperature metrics
  avTa <- mean(ta, na.rm = TRUE)
  SDTa <- sd(ta, na.rm = TRUE)

  ## forest surface temperature metrics
  avTf <- mean(tf, na.rm = TRUE)
  SDTf <- sd(tf, na.rm = TRUE)

  # vertical gradients
  sel_v <- res$micro_grid$x == 30 & res$micro_grid$y == 30
  zs_v <- res$micro_grid$z[sel_v]
  ta_v <- res$air_temperature[sel_v] - 273.15 # air temperature
  tf_v <- res$micro_grid$temperature[sel_v] - 273.15 # forest surface temperature

  # air temperature metrics
  avTa_v <- mean(ta_v, na.rm = TRUE)
  gradTa_v <- (ta_v[zs_v == max(zs_v)][1] - ta_v[zs_v == min(zs_v)][1]) / (max(zs_v) - min(zs_v))
  SDTa_v <- sd(ta_v, na.rm = TRUE)

  # forest surface temperature metrics
  avTf_v <- mean(tf_v, na.rm = TRUE)
  gradTf_v <- (tf_v[zs_v == max(zs_v)][1] - tf_v[zs_v == min(zs_v)][1]) / (max(zs_v) - min(zs_v))
  SDTf_v <- sd(tf_v, na.rm = TRUE)

  return(c(avTa = avTa, SDTa = SDTa, avTa_v = avTa_v, gradTa_v = gradTa_v, SDTa_v = SDTa_v, avTf = avTf, SDTf = SDTf, avTf_v = avTf_v, gradTf_v = gradTf_v, SDTf_v = SDTf_v))
}


##############
# SOBOL RUNS #
##############

sobol_outputs_all <- with_progress({
  p <- progressor(steps = nrow(scaled_input))
  future_lapply(seq_len(nrow(scaled_input)), function(i) {
    metrics <- model_function(unlist(scaled_input[i, ]))
    p()
    return(metrics)
  }, future.seed = TRUE)
})

plan(sequential)

# Metric extraction
# whole grid temperatures
# air temperature
sobol_outputs_metric1 <- map_dbl(sobol_outputs_all, "avTa")
sobol_outputs_metric2 <- map_dbl(sobol_outputs_all, "SDTa")

sobol_metric1 <- tell(sobol_design, sobol_outputs_metric1)
sobol_metric2 <- tell(sobol_design, sobol_outputs_metric2)

# forest surface temperature
sobol_outputs_metric3 <- map_dbl(sobol_outputs_all, "avTf")
sobol_outputs_metric4 <- map_dbl(sobol_outputs_all, "SDTf")

sobol_metric3 <- tell(sobol_design, sobol_outputs_metric3)
sobol_metric4 <- tell(sobol_design, sobol_outputs_metric4)

# vertical gradients
# air temperature
sobol_outputs_metric5 <- map_dbl(sobol_outputs_all, "avTa_v")
sobol_outputs_metric6 <- map_dbl(sobol_outputs_all, "gradTa_v")
sobol_outputs_metric7 <- map_dbl(sobol_outputs_all, "SDTa_v")

sobol_metric5 <- tell(sobol_design, sobol_outputs_metric5)
sobol_metric6 <- tell(sobol_design, sobol_outputs_metric6)
sobol_metric7 <- tell(sobol_design, sobol_outputs_metric7)

# forest surface temperature
sobol_outputs_metric8 <- map_dbl(sobol_outputs_all, "avTf_v")
sobol_outputs_metric9 <- map_dbl(sobol_outputs_all, "gradTf_v")
sobol_outputs_metric10 <- map_dbl(sobol_outputs_all, "SDTf_v")

sobol_metric8 <- tell(sobol_design, sobol_outputs_metric8)
sobol_metric9 <- tell(sobol_design, sobol_outputs_metric9)
sobol_metric10 <- tell(sobol_design, sobol_outputs_metric10)


##########
# SAVING #
##########

# whole grid temperatures
# air temperature
saveRDS(sobol_metric1, file = file.path(output_path, paste0(output_name,"_", plotnr, "_avTa.rds")))
saveRDS(sobol_metric2, file = file.path(output_path, paste0(output_name,"_", plotnr, "_SDTa.rds")))

# forest surface temperature
saveRDS(sobol_metric3, file = file.path(output_path, paste0(output_name,"_", plotnr, "_avTf.rds")))
saveRDS(sobol_metric4, file = file.path(output_path, paste0(output_name,"_", plotnr, "_SDTf.rds")))

# vertical gradient
# air temperature
saveRDS(sobol_metric5, file = file.path(output_path, paste0(output_name,"_", plotnr, "_avTa_v.rds")))
saveRDS(sobol_metric6, file = file.path(output_path, paste0(output_name,"_", plotnr, "_gradTa_v.rds")))
saveRDS(sobol_metric7, file = file.path(output_path, paste0(output_name,"_", plotnr, "_SDTa_v.rds")))

# forest surface temperature
saveRDS(sobol_metric8, file = file.path(output_path, paste0(output_name,"_", plotnr, "_avTf_v.rds")))
saveRDS(sobol_metric9, file = file.path(output_path, paste0(output_name,"_", plotnr, "_gradTf_v.rds")))
saveRDS(sobol_metric10, file = file.path(output_path, paste0(output_name,"_", plotnr, "_SDTf_v.rds")))
