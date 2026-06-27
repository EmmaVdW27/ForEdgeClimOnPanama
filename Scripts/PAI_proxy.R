files <- c(
  "Output/2023-04-26_19/C8/model_results.rds",
  "Output/2023-04-26_19/C13/model_results.rds",
  "Output/2023-04-26_19/R7/model_results.rds",
  "Output/2023-04-26_19/R11/model_results.rds"
)
library(dplyr)
library(stringr)
library(purrr)

get_pai_proxy <- function(file){

  model_results <- readRDS(file)

  lai_proxy <- model_results$sw_rad %>%
    group_by(X, Y) %>%
    summarise(
      lai_proxy = sum(density, na.rm = TRUE),
      .groups = "drop"
    )

  plot_name <- str_extract(file, "(C8|C13|R7|R11)")

  tibble(
    plot = plot_name,
    mean_PAI_proxy   = mean(lai_proxy$lai_proxy, na.rm = TRUE),
    median_PAI_proxy = median(lai_proxy$lai_proxy, na.rm = TRUE),
    sd_PAI_proxy     = sd(lai_proxy$lai_proxy, na.rm = TRUE),
    min_PAI_proxy    = min(lai_proxy$lai_proxy, na.rm = TRUE),
    max_PAI_proxy    = max(lai_proxy$lai_proxy, na.rm = TRUE)
  )
}

results <- map_dfr(files, get_pai_proxy)

print(results)


# # global statistics
# df <- readRDS("Output/timeseries/data/timeseries_R7_uncalibrated.rds")
#
# valid <- complete.cases(df$radiation_mean, df$PAR_model)
#
# obs <- df$radiation_mean[valid]
# sim <- df$PAR_model[valid]
#
# RMSE <- sqrt(mean((sim - obs)^2))
#
# R2 <- cor(obs, sim)^2
#
# NSE <- 1 - sum((sim - obs)^2) /
#   sum((obs - mean(obs))^2)
#
# ME <- mean(sim - obs)
#
# SD <- sd(sim - obs)
#
# print(c(RMSE = RMSE,
#   R2 = R2,
#   NSE = NSE,
#   ME = ME,
#   SD = SD))
