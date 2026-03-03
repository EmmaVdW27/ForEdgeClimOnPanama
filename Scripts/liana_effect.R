###############################################################################
# Script to investigate whether the observed liana effect (cooler temperatures in liana plots) can be modelled.
# Liana removal plot R11 and control plot C13 were used for this investigation.
#
# Author: Emma Van de Walle - Q-ForestLab
#
###############################################################################

library(dplyr)
library(lubridate)
library(ggplot2)

#########
# INPUT #
#########

# micro temperature modelled and observed
input_model_C13 <- readRDS("Output/timeseries/data/timeseries_C13_calibrated.rds")
input_model_R11 <- readRDS("Output/timeseries/data/timeseries_R11_calibrated.rds")

# macro temperature
input_temps <- readRDS("Data/TMS_macro_2023_C13_and_R11.rds")

# output path
output_path = "Output/liana-effect/"

#########

# set datetime
input_model_C13 <- input_model_C13 %>%
  mutate(datetime = ymd_hms(datetime))

input_model_R11 <- input_model_R11 %>%
  mutate(datetime = ymd_hms(datetime))

# mean over sensors
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

# join dataframes
comparison <- C13_plot %>%
  inner_join(R11_plot, by = "datetime") #%>%
  # filter on day or night hours
  #filter(hour(datetime) >= 21 | hour(datetime) <= 3) # night
  #filter(hour(datetime) >= 9 & hour(datetime) <= 15) # day

# liana effect comparison
comparison <- comparison %>%
  mutate(
    delta_obs  = obs_R11 - obs_C13,
    delta_mod  = mod_R11 - mod_C13,
    delta_diff = delta_mod - delta_obs   # model bias in the liana-effect
  )

# PLOTS
#######

# Scatter plot
scatter = ggplot(comparison, aes(x = delta_obs, y = delta_mod)) +
  geom_point(alpha = 0.5, colour = "blue") +
  geom_abline(slope = 1, intercept = 0, colour = "red", linewidth = 1.2) +
  labs(
    x = "Observed liana effect (R11 - C13, °C)",
    y = "Modelled liana effect (R11 - C13, °C)"
  ) +
  theme_bw()

print(scatter)
ggsave(paste0(output_path, "liana-effect_scatter.png"), plot = scatter, width = 8, height = 4, dpi = 500)

# Timeseries plot
# Calculate means
mean_obs <- mean(comparison$delta_obs, na.rm = TRUE)
mean_mod <- mean(comparison$delta_mod, na.rm = TRUE)

# Left boundary for placing text
xmin <- min(comparison$datetime, na.rm = TRUE)

timeseries <- ggplot(comparison, aes(x = datetime)) +

  # graphs
  geom_line(aes(y = delta_obs, colour = "Observed"),
            linewidth = 0.6, alpha = 0.4) +
  geom_line(aes(y = delta_mod, colour = "Modelled"),
            linewidth = 0.6, alpha = 0.4) +

  # mean lines
  geom_hline(aes(yintercept = mean_obs, colour = "Observed"),
             linewidth = 1.2) +
  geom_hline(aes(yintercept = mean_mod, colour = "Modelled"),
             linewidth = 1.2) +

  # mean value labels (boven de lijn, links)
  annotate("text",
           x = xmin,
           y = mean_obs,
           label = paste0("Mean = ", round(mean_obs, 2), " °C"),
           colour = "#0047AB",
           hjust = 0,
           vjust = -0.5,
           size = 5) +

  annotate("text",
           x = xmin,
           y = mean_mod,
           label = paste0("Mean = ", round(mean_mod, 2), " °C"),
           colour = "#B22222",
           hjust = 0,
           vjust = +1.5,
           size = 5) +

  # colours
  scale_colour_manual(values = c("Observed" = "#0047AB",
                                 "Modelled" = "#B22222")) +

  labs(y = "Temperature difference (R11 - C13, °C)",
       colour = "") +
  theme_bw()

print(timeseries)

ggsave(paste0(output_path, "liana-effect_timeseries.png"),
       plot = timeseries,
       width = 10, height = 5, dpi = 500)

# RMSE of liana-effect
print(sqrt(mean((comparison$delta_mod - comparison$delta_obs)^2, na.rm = TRUE)))

# correlation of liana-effect
print(cor(comparison$delta_obs, comparison$delta_mod, use = "complete.obs"))
