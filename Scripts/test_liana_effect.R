###############################################################################
# Liana effect analysis (temperature + process variables)
###############################################################################

library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)

#########
# INPUT #
#########

input_model_C13 <- readRDS("Output/timeseries/data/v1.0/timeseries_C13_v1.0_calibrated.rds")
input_model_R11 <- readRDS("Output/timeseries/data/v1.0/timeseries_R11_v1.0_calibrated.rds")

input_model_C13_v2 <- readRDS("Output/timeseries/data/timeseries_C13_semi_calibrated_subset1.rds")
input_model_R11_v2 <- readRDS("Output/timeseries/data/timeseries_R11_semi_calibrated_subset1.rds")

output_path <- "Output/liana-effect/"

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
# AGGREGATION (TEMPERATURE)
#########

C13_plot <- input_model_C13 %>%
  group_by(datetime) %>%
  summarise(
    obs_C13 = mean(value_mean, na.rm = TRUE),
    mod_C13 = mean(temperature_model, na.rm = TRUE),
    .groups = "drop"
  )

R11_plot <- input_model_R11 %>%
  group_by(datetime) %>%
  summarise(
    obs_R11 = mean(value_mean, na.rm = TRUE),
    mod_R11 = mean(temperature_model, na.rm = TRUE),
    .groups = "drop"
  )

C13_plot_v2 <- input_model_C13_v2 %>%
  group_by(datetime) %>%
  summarise(
    mod_C13_v2 = mean(temperature_model, na.rm = TRUE),
    .groups = "drop"
  )

R11_plot_v2 <- input_model_R11_v2 %>%
  group_by(datetime) %>%
  summarise(
    mod_R11_v2 = mean(temperature_model, na.rm = TRUE),
    .groups = "drop"
  )

#########
# EXTRA VARIABLES (v1.1 ONLY)
#########

C13_extra_v2 <- input_model_C13_v2 %>%
  group_by(datetime) %>%
  summarise(
    LE_C13    = first(LE_mean),
    RH_C13    = first(RH_mean),
    VPD_C13   = first(VPD_mean),
    Fd_C13    = first(Fd_down_mean),
    Fb_C13    = first(Fb_down_mean),
    Tsurf_C13 = first(Tsurf_mean),
    .groups = "drop"
  )

R11_extra_v2 <- input_model_R11_v2 %>%
  group_by(datetime) %>%
  summarise(
    LE_R11    = first(LE_mean),
    RH_R11    = first(RH_mean),
    VPD_R11   = first(VPD_mean),
    Fd_R11    = first(Fd_down_mean),
    Fb_R11    = first(Fb_down_mean),
    Tsurf_R11 = first(Tsurf_mean),
    .groups = "drop"
  )

#########
# JOIN
#########

comparison <- C13_plot %>%
  inner_join(R11_plot, by = "datetime") %>%
  inner_join(C13_plot_v2, by = "datetime") %>%
  inner_join(R11_plot_v2, by = "datetime") %>%
  inner_join(C13_extra_v2, by = "datetime") %>%
  inner_join(R11_extra_v2, by = "datetime")

#########
# REMOVE PERIOD
#########

exclude_start <- ymd_hms("2023-06-21 00:00:00")
exclude_end   <- ymd_hms("2023-06-27 23:59:59")

comparison <- comparison %>%
  mutate(
    across(
      where(is.numeric),
      ~ ifelse(
        datetime >= exclude_start &
          datetime <= exclude_end,
        NA,
        .
      )
    )
  )

#########
# DELTAS
#########

comparison <- comparison %>%
  mutate(

    delta_obs    = obs_R11 - obs_C13,
    delta_mod    = mod_R11 - mod_C13,
    delta_mod_v2 = mod_R11_v2 - mod_C13_v2,

    delta_LE    = LE_R11 - LE_C13,
    delta_RH    = RH_R11 - RH_C13,
    delta_VPD   = VPD_R11 - VPD_C13,
    delta_Fd    = Fd_R11 - Fd_C13,
    delta_Fb    = Fb_R11 - Fb_C13,
    delta_Tsurf = Tsurf_R11 - Tsurf_C13

  )

#########
# DAY / NIGHT
#########

comparison <- comparison %>%
  mutate(

    hour_utc = hour(datetime),

    period = case_when(
      hour_utc >= 23 | hour_utc < 11 ~ "Night",
      TRUE ~ "Day"
    ),

    period = factor(
      period,
      levels = c("Day", "Night")
    )

  )

###############################################################################
# TIMESERIES PLOT WITH DAY/NIGHT FACETS
###############################################################################

timeseries_long <- comparison %>%
  select(
    datetime,
    period,
    delta_obs,
    delta_mod,
    delta_mod_v2
  ) %>%
  pivot_longer(
    cols = c(delta_obs, delta_mod, delta_mod_v2),
    names_to = "type",
    values_to = "value"
  ) %>%
  mutate(
    type = recode(
      type,
      delta_obs    = "Observed",
      delta_mod    = "Modelled v1.0",
      delta_mod_v2 = "Modelled v1.1"
    )
  )

means_ts <- timeseries_long %>%
  group_by(type, period) %>%
  summarise(
    mean_val = mean(value, na.rm = TRUE),
    .groups = "drop"
  )

xmin_ts <- min(timeseries_long$datetime, na.rm = TRUE)

timeseries <- ggplot(
  timeseries_long,
  aes(x = datetime, y = value)
) +

  geom_point(
    aes(colour = type),
    size = 1.2,
    alpha = 0.8,
    na.rm = TRUE
  ) +

  # witte achtergrondlijnen
  geom_hline(
    data = means_ts,
    aes(yintercept = mean_val),
    colour = "white",
    linewidth = 3
  ) +

  # dashed mean lijnen
  geom_hline(
    data = means_ts,
    aes(
      yintercept = mean_val,
      colour = type
    ),
    linetype = "dashed",
    linewidth = 1.2
  ) +

  geom_label(
    data = means_ts,
    aes(
      x = xmin_ts,
      y = mean_val,
      label = paste0(
        type,
        ": ",
        round(mean_val, 2),
        " °C"
      ),
      colour = type
    ),
    fill = "white",
    hjust = 0,
    vjust = -0.5,
    size = 4,
    fontface = "bold",
    inherit.aes = FALSE
  ) +

  scale_colour_manual(
    values = c(
      "Observed" = "blue",
      "Modelled v1.0" = "orange",
      "Modelled v1.1" = "red"
    )
  ) +

  facet_wrap(
    ~period,
    ncol = 2
  ) +

  labs(
    x = "Datetime",
    y = "Air temperature difference (R11 - C13, °C)",
    colour = NULL
  ) +

  theme_bw() +

  theme(
    strip.text = element_text(
      size = 15,
      face = "bold"
    ),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12),
    legend.position = "top"
  )

print(timeseries)

ggsave(
  paste0(output_path,
         "liana-effect_timeseries_v1_v2_daynight.png"),
  plot = timeseries,
  width = 14,
  height = 6,
  dpi = 500
)

###############################################################################
# FACET PLOTS PER PROCESS GROUP
###############################################################################

#########
# LONG FORMAT
#########

comparison_long <- comparison %>%
  select(
    datetime,
    period,
    delta_LE,
    delta_RH,
    delta_VPD,
    delta_Fd,
    delta_Fb,
    delta_Tsurf
  ) %>%
  pivot_longer(
    cols = c(
      delta_LE,
      delta_RH,
      delta_VPD,
      delta_Fd,
      delta_Fb,
      delta_Tsurf
    ),
    names_to = "variable",
    values_to = "value"
  )

#########
# LABELS
#########

comparison_long <- comparison_long %>%
  mutate(
    variable = recode(
      variable,
      delta_LE    = "a) Δ LE (W m⁻²)",
      delta_RH    = "a) Δ RH (%)",
      delta_VPD   = "b) Δ VPD (kPa)",
      delta_Fd    = "a) Δ Diffuse SW (W m⁻²)",
      delta_Fb    = "b) Δ Direct SW (W m⁻²)",
      delta_Tsurf = "b) Δ Surface T (°C)"
    )
  )

plot_facet <- function(df, vars, title, filename) {

  df_sub <- df %>%
    filter(variable %in% vars)

  #########
  # MEANS
  #########

  means <- df_sub %>%
    group_by(variable, period) %>%
    summarise(
      mean_val = mean(value, na.rm = TRUE),
      .groups = "drop"
    )

  xmin <- min(df_sub$datetime, na.rm = TRUE)

  #########
  # PLOT
  #########

  p <- ggplot(
    df_sub,
    aes(x = datetime, y = value)
  ) +

    geom_point(
      colour = "blue",
      size = 1.2,
      alpha = 0.8,
      na.rm = TRUE
    ) +

    # witte achtergrondlijn
    geom_hline(
      data = means,
      aes(yintercept = mean_val),
      colour = "white",
      linewidth = 3
    ) +

    # rode dashed lijn
    geom_hline(
      data = means,
      aes(yintercept = mean_val),
      colour = "red",
      linetype = "dashed",
      linewidth = 1
    ) +

    geom_label(
      data = means,
      aes(
        x = xmin,
        y = mean_val,
        label = paste0(
          "Mean = ",
          round(mean_val, 2)
        )
      ),
      fill = "white",
      colour = "red",
      hjust = 0,
      vjust = -0.5,
      size = 4,
      fontface = "bold",
      inherit.aes = FALSE
    ) +

    # variabelen onder elkaar
    # day/night als kolommen
    facet_grid(
      variable ~ period,
      scales = "free_y"
    ) +

    scale_x_datetime(
      limits = c(
        min(df_sub$datetime),
        max(df_sub$datetime) + days(10)
      )
    ) +

    theme_bw() +

    labs(
      title = title,
      x = "Datetime",
      y = NULL
    ) +

    theme(
      strip.text.y = element_text(
        size = 14,
        face = "bold",
        angle = 0
      ),
      strip.text.x = element_text(
        size = 14,
        face = "bold"
      ),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 12),
      plot.title = element_text(size = 18)
    )

  print(p)

  ggsave(
    paste0(output_path, filename),
    plot = p,
    width = 13,
    height = 8,
    dpi = 300
  )
}

###############################################################################
# 1. RADIATION
###############################################################################

plot_facet(
  comparison_long,
  vars = c(
    "a) Δ Diffuse SW (W m⁻²)",
    "b) Δ Direct SW (W m⁻²)"
  ),
  title = NULL,
  filename = "facet_radiation_daynight.png"
)

###############################################################################
# 2. HUMIDITY
###############################################################################

plot_facet(
  comparison_long,
  vars = c(
    "a) Δ RH (%)",
    "b) Δ VPD (kPa)"
  ),
  title = NULL,
  filename = "facet_humidity_daynight.png"
)

###############################################################################
# 3. LATENT HEAT & SURFACE TEMPERATURE
###############################################################################

plot_facet(
  comparison_long,
  vars = c(
    "a) Δ LE (W m⁻²)",
    "b) Δ Surface T (°C)"
  ),
  title = NULL,
  filename = "facet_LE_Tf_daynight.png"
)

###############################################################################
# METRICS
###############################################################################

rmse_v2 <- sqrt(
  mean(
    (comparison$delta_mod_v2 - comparison$delta_obs)^2,
    na.rm = TRUE
  )
)

cor_v2 <- cor(
  comparison$delta_obs,
  comparison$delta_mod_v2,
  use = "complete.obs"
)

print(paste("RMSE v1.1:", rmse_v2))
print(paste("Correlation v1.1:", cor_v2))
