###############################################################################
# Script to investigate whether the observed liana effect (cooler temperatures
# in liana plots) can be modelled.
# Liana removal plot R11 and control plot C13 were used for this investigation.
#
# Author: Emma Van de Walle - Q-ForestLab
###############################################################################

library(dplyr)
library(lubridate)
library(ggplot2)

#########
# INPUT #
#########

input_model_C13 <- readRDS("Output/timeseries/data/timeseries_C13_v1.0_calibrated.rds")
input_model_R11 <- readRDS("Output/timeseries/data/timeseries_R11_v1.0_calibrated.rds")

input_model_C13_v2 <- readRDS("Output/timeseries/data/timeseries_C13_calibrated_PM_and_u.rds")
input_model_R11_v2 <- readRDS("Output/timeseries/data/timeseries_R11_calibrated_PM_and_u.rds")

output_path = "Output/liana-effect/"

#########
# CLEAN DATETIME
#########

clean_datetime <- function(df) {
  df %>%
    mutate(
      datetime = trimws(datetime),
      datetime = na_if(datetime, ""),
      datetime = na_if(datetime, "NA"),
      datetime = suppressWarnings(ymd_hms(datetime))
    ) %>%
    filter(!is.na(datetime))
}

input_model_C13    <- clean_datetime(input_model_C13)
input_model_R11    <- clean_datetime(input_model_R11)
input_model_C13_v2 <- clean_datetime(input_model_C13_v2)
input_model_R11_v2 <- clean_datetime(input_model_R11_v2)

#########
# AGGREGATION
#########

C13_plot <- input_model_C13 %>%
  group_by(datetime) %>%
  summarise(
    obs_C13 = mean(value_mean, na.rm = TRUE),
    mod_C13 = mean(temperature_model, na.rm = TRUE)
  )

R11_plot <- input_model_R11 %>%
  group_by(datetime) %>%
  summarise(
    obs_R11 = mean(value_mean, na.rm = TRUE),
    mod_R11 = mean(temperature_model, na.rm = TRUE)
  )

C13_plot_v2 <- input_model_C13_v2 %>%
  group_by(datetime) %>%
  summarise(
    mod_C13_v2 = mean(temperature_model, na.rm = TRUE)
  )

R11_plot_v2 <- input_model_R11_v2 %>%
  group_by(datetime) %>%
  summarise(
    mod_R11_v2 = mean(temperature_model, na.rm = TRUE)
  )

#########
# JOIN
#########

comparison <- C13_plot %>%
  inner_join(R11_plot, by = "datetime") %>%
  inner_join(C13_plot_v2, by = "datetime") %>%
  inner_join(R11_plot_v2, by = "datetime")

#########
# REMOVE PERIOD (NO GAPS)
#########

exclude_start <- ymd_hms("2023-06-21 00:00:00")
exclude_end   <- ymd_hms("2023-06-27 23:59:59")

comparison <- comparison %>%
  mutate(
    across(where(is.numeric), ~ ifelse(datetime >= exclude_start & datetime <= exclude_end, NA, .))
  )

#########
# DELTAS
#########

comparison <- comparison %>%
  mutate(
    delta_obs    = obs_R11 - obs_C13,
    delta_mod    = mod_R11 - mod_C13,
    delta_mod_v2 = mod_R11_v2 - mod_C13_v2
  )

###################
# TIMESERIES PLOT #
###################

mean_obs    <- mean(comparison$delta_obs, na.rm = TRUE)
mean_mod    <- mean(comparison$delta_mod, na.rm = TRUE)
mean_mod_v2 <- mean(comparison$delta_mod_v2, na.rm = TRUE)

xmin <- min(comparison$datetime, na.rm = TRUE)
xmax <- max(comparison$datetime, na.rm = TRUE)

timeseries <- ggplot(comparison, aes(x = datetime)) +

  # lijnen
  geom_line(aes(y = delta_obs, colour = "Observed"),
            linewidth = 0.8, alpha = 0.9) +
  geom_line(aes(y = delta_mod_v2, colour = "Modelled v1.1"),
            linewidth = 0.8, alpha = 0.9) +
  geom_line(aes(y = delta_mod, colour = "Modelled v1.0"),
            linewidth = 0.8, alpha = 0.9) +


  # mean lines → dashed + dikker + bovenop
  # WITTE HALO (achtergrond)
  geom_hline(yintercept = mean_obs,
             colour = "white", linewidth = 3) +
  geom_hline(yintercept = mean_mod,
             colour = "white", linewidth = 3) +
  geom_hline(yintercept = mean_mod_v2,
             colour = "white", linewidth = 3) +

  # ECHTE MEAN LINES (bovenop)
  geom_hline(yintercept = mean_obs,
             colour = "blue", linetype = "dashed", linewidth = 1.5) +
  geom_hline(yintercept = mean_mod,
             colour = "orange", linetype = "dashed", linewidth = 1.5) +
  geom_hline(yintercept = mean_mod_v2,
             colour = "red", linetype = "dashed", linewidth = 1.5) +

  # labels met witte achtergrond + outline
  geom_label(
    data = data.frame(x = xmin, y = mean_obs),
    aes(x = x, y = y,
        label = paste0("Mean obs = ", round(mean_obs, 2), " °C")),
    fill = "white",
    colour = "blue",
    hjust = 0,
    vjust = -0.5,
    size = 6,
    label.size = 0.3,
    fontface = "bold",
    inherit.aes = FALSE
  ) +

  geom_label(
    data = data.frame(x = xmax, y = mean_mod),
    aes(x = x, y = y,
        label = paste0("Mean v1.0 = ", round(mean_mod, 2), " °C")),
    fill = "white",
    colour = "orange",
    hjust = 1,
    vjust = 1.5,
    size = 6,
    label.size = 0.3,
    fontface = "bold",
    inherit.aes = FALSE
  ) +

  geom_label(
    data = data.frame(x = xmax, y = mean_mod_v2),
    aes(x = x, y = y,
        label = paste0("Mean v1.1 = ", round(mean_mod_v2, 2), " °C")),
    fill = "white",
    colour = "red",
    hjust = 1,
    vjust = -0.5,
    size = 6,
    label.size = 0.3,
    fontface = "bold",
    inherit.aes = FALSE
  ) +

  # betere kleuren (colorblind safe)
  scale_colour_manual(values = c(
    "Observed" = "blue",   # blauw
    "Modelled v1.0" = "orange",# oranje
    "Modelled v1.1" = "red" # red
  )) +

  labs(
    title = "Liana effect: observed vs modelled",
    x = "Datetime",
    y = "Temperature difference (R11 - C13, °C)",
    colour = NULL
  ) +

  theme_bw() +
  theme(
    plot.title = element_text(size = 22),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.position = "top",
    legend.key.width = unit(2, "cm")
  )

print(timeseries)

ggsave(paste0(output_path, "liana-effect_timeseries_v1_v2.png"),
       plot = timeseries,
       width = 12, height = 6, dpi = 500)

##########
# METRICS #
##########

rmse_v1 <- sqrt(mean((comparison$delta_mod - comparison$delta_obs)^2, na.rm = TRUE))
rmse_v2 <- sqrt(mean((comparison$delta_mod_v2 - comparison$delta_obs)^2, na.rm = TRUE))

print(paste("RMSE v1.0:", rmse_v1))
print(paste("RMSE v1.1:", rmse_v2))

cor_v1 <- cor(comparison$delta_obs, comparison$delta_mod, use = "complete.obs")
cor_v2 <- cor(comparison$delta_obs, comparison$delta_mod_v2, use = "complete.obs")

print(paste("Correlation v1.0:", cor_v1))
print(paste("Correlation v1.1:", cor_v2))
