library(ForEdgeClimOnPanama)

# replace all "R7" by "C8" or vice versa

folders <- paste0(
  "Output/",
  sprintf("2023-04-26_%02d", 11:22), "/R7"
)

results <- lapply(
  file.path(folders, "model_results.Rds"),
  readRDS
)

n_timesteps <- length(results)

# Use first result as template
mean_results <- results[[1]]


mean_results$micro_grid$temperature <- Reduce(
  "+",
  lapply(results, \(x) x$micro_grid$temperature)
) / n_timesteps

mean_results$micro_grid$T_soil <- Reduce(
  "+",
  lapply(results, \(x) x$micro_grid$T_soil)
) / n_timesteps

mean_results$air_temperature <- Reduce(
  "+",
  lapply(results, \(x) x$air_temperature)
) / n_timesteps

mean_results$sw_rad <- Reduce(
  "+",
  lapply(results, \(x) x$sw_rad)
) / n_timesteps

den <<- mean_results$sw_rad$density

dir.create(
  "Output/2023-04-26_mean/R7",
  recursive = TRUE,
  showWarnings = FALSE
)

saveRDS(
  mean_results,
  "Output/2023-04-26_mean/R7/model_results.Rds"
)

res <- readRDS("Output/2023-04-26_mean/R7/model_results.Rds")

plots_sw(
  sw_rad = res$sw_rad,
  output_path = "Output/2023-04-26_mean/R7"
)

plots_temp(
  micro_grid = res$micro_grid,
  T_air_vec = res$air_temperature,
  output_path = "Output/2023-04-26_mean/R7",
  datetime = "2023-04-26"
)
