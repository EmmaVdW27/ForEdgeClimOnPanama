###############################################################################
#
# This is an example script to show you how the microclimate model ForEdgeClim
# can be run. To run the model, you need to define several variables, input and
# output files.
#
# Author: Emma Van de Walle - Q-ForestLab
#
###############################################################################


library(ForEdgeClimOnPanama)
start_script = Sys.time()

#########
# INPUT #
#########
plotnr <- "C13" # C13 or R11

start_time <- as.POSIXct("2023-08-12 19:00:00", tz = "UTC") # first hour you want to run the model for
end_time <- as.POSIXct("2023-08-12 19:00:00", tz = "UTC")   # last hour you want to run the model for
datetime_series <- seq(start_time, end_time, by = "hour")   # ForEdgeClim is run for all the hours between the first and last hour

#########

if(plotnr == "C13"){
  plot_map = "C13"
  height_canopy = 33
  # original TLS las file:
  TLS_input_file <- 'Data/Data_Kasper/plot13_trimmed_5cm.las'
  # Once the structural voxel grid has been made using the function generate_DTM_grid_TLS() (see further below),
  # you can also immediately use the TLS_filtered_file:
  TLS_filtered_file <- 'Data/TLS_scaled_DTM_and_grid_C13_Panama.rds'
} else { # plotnr == "R11"
  plot_map = "R11"
  height_canopy = 35
  TLS_input_file <- 'Data/Data_Kasper/plot11_trimmed_5cm.las'
  TLS_filtered_file <- 'Data/TLS_scaled_DTM_and_grid_R11_Panama.rds'
}

for (current_datetime in datetime_series) {

  # Make sure the current datetime is a POSIXct object
  current_datetime <- as.POSIXct(current_datetime, origin = "1970-01-01", tz = "UTC")
  cat(paste0('Running model for plot ', plotnr, ' at ', current_datetime, ' UTC...\n'))

  #################################################
  # INPUT DRIVERS, CONSTANTS AND MODEL PARAMETERS #
  #################################################

  create_input_drivers()
  create_physical_constants()
  create_model_parameters()

  # Create output path
  output_path <- file.path("Output", format(current_datetime, "%Y-%m-%d_%H"), plot_map)
  # Make sure the output path exists
  if (!dir.exists(output_path)) {
    dir.create(output_path, recursive = TRUE, showWarnings = FALSE)
  }

  ########################################################
  # EXTRACT THE OBSERVATIONS FROM THE INPUT DRIVER FILES #
  ########################################################

  # Import observations as input variables and as variables to compare the model with
  import_temp_observations(current_datetime)
  import_rad_observations(current_datetime)

  #######################
  # CREATION OF 3D GRID #
  #######################

  # Structure from TLS las file
  # voxel_TLS = generate_DTM_grid_TLS(las_file = TLS_input_file, voxel_size = voxel_length)
  # saveRDS(voxel_TLS, TLS_filtered_file)
  voxel_TLS = readRDS(TLS_filtered_file)

  #############
  # RUN MODEL #
  #############

  res = run_foredgeclim(voxel_TLS$grid, current_datetime)
  saveRDS(res, paste0(output_path, '/model_results.rds'))

  #--------------------
  ########################################################
  # COMPARE MODELLED AND OBSERVED TEMPERATURES (TOMST)  #
  ########################################################

  library(RANN)
  library(dplyr)

  # ---- Import TOMST observations ----
  key <- format(current_datetime, "%Y%m%d_%H")
  TOMST <- read.csv(paste0("Data/microtemp_air_obs_", key, "_", plotnr, ".csv"))

  # ---- Extract modelled temperatures at z = 1 ----
  temp_air_grid <- res$micro_grid
  temp_air_grid$temperature <- res$air_temperature - 273.15

  temp_air_grid <- temp_air_grid %>%
    filter(z == 1) %>%
    select(x, y, temperature)

  # ---- Nearest neighbour search ----

  # matrix with grid coordinates
  grid_coords <- as.matrix(temp_air_grid[, c("x", "y")])

  # matrix with sensor coordinates
  sensor_coords <- as.matrix(TOMST[, c("x", "y")])

  # nearest neighbour
  nn <- nn2(
    data = grid_coords,
    query = sensor_coords,
    k = 1
  )

  idx <- nn$nn.idx[, 1]

  # ---- Add modelled temperatures to observations ----
  TOMST <- TOMST %>%
    mutate(
      datetime = current_datetime,
      temperature_model = temp_air_grid$temperature[idx],
      resid = temperature_model - value_mean
    )

  # ---- Calculate mean temperatures ----
  mean_model <- mean(TOMST$temperature_model, na.rm = TRUE)
  mean_obs   <- mean(TOMST$value_mean, na.rm = TRUE)

  cat("Mean model temperature:", round(mean_model, 2), "°C\n")
  cat("Mean observed temperature:", round(mean_obs, 2), "°C\n")

  #--------------------------


  # ############
  # # PLOTTING #
  # ############
  #
  # # Digital Terrain Model & structural grid plots
  # plots_dtm_struct(dtm = voxel_TLS$dtm, grid = voxel_TLS$grid, output_path)
  #
  # # Shortwave radiation plots
  # plots_sw(sw_rad = res$sw_rad, output_path)
  #
  # # Longwave radiation plots
  # plots_lw(lw_rad = res$lw_rad, output_path)
  #
  # # Flux plots
  # plots_flux(res$micro_grid, res$net_radiation, res$sensible_flux, res$latent_flux, res$ground_flux, output_path, current_datetime)
  #
  # # Temperature plots
  # plots_temp(res$micro_grid, res$air_temperature, output_path, current_datetime)

}
end_script = Sys.time()
cat(paste0('Total run time = ', round(as.numeric(end_script - start_script, units = "secs"), 2), ' s'))
