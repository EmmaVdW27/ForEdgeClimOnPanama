###############################################################################
# Liana effect analysis (temperature + process variables)
# Separate day/night analysis (Panama UTC)
###############################################################################

library(dplyr)
library(lubridate)
library(ggplot2)
library(tidyr)
library(patchwork)
library(grid)

#########
# INPUT #
#########

input_model_C13 <- readRDS(
  "Output/timeseries/data/v1.0/timeseries_C13_v1.0_calibrated.rds"
)

input_model_R11 <- readRDS(
  "Output/timeseries/data/v1.0/timeseries_R11_v1.0_calibrated.rds"
)

input_model_C13_v2 <- readRDS(
  "Output/timeseries/data/timeseries_C13_30_gens_5_params_calibrated.rds"
)

input_model_R11_v2 <- readRDS(
  "Output/timeseries/data/timeseries_R11_30_gens_5_params_calibrated.rds"
)

output_path <- "Output/liana-effect/"

###############################################################################
# CLEAN DATETIME
###############################################################################

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

###############################################################################
# DAY/NIGHT CLASSIFICATION (Panama, UTC)
###############################################################################

add_period <- function(df) {

  df %>%
    mutate(
      hour_utc = hour(datetime),

      period = case_when(
        hour_utc >= 11 & hour_utc < 23 ~ "Day",
        TRUE                           ~ "Night"
      )
    )

}

input_model_C13    <- add_period(input_model_C13)
input_model_R11    <- add_period(input_model_R11)
input_model_C13_v2 <- add_period(input_model_C13_v2)
input_model_R11_v2 <- add_period(input_model_R11_v2)

###############################################################################
# AGGREGATION (TEMPERATURE)
###############################################################################

C13_plot <- input_model_C13 %>%
  mutate(date = as.Date(datetime)) %>%
  group_by(date, period) %>%
  summarise(
    obs_C13 = mean(value_mean, na.rm = TRUE),
    mod_C13 = mean(temperature_model, na.rm = TRUE),
    .groups = "drop"
  )

R11_plot <- input_model_R11 %>%
  mutate(date = as.Date(datetime)) %>%
  group_by(date, period) %>%
  summarise(
    obs_R11 = mean(value_mean, na.rm = TRUE),
    mod_R11 = mean(temperature_model, na.rm = TRUE),
    .groups = "drop"
  )

C13_plot_v2 <- input_model_C13_v2 %>%
  mutate(date = as.Date(datetime)) %>%
  group_by(date, period) %>%
  summarise(
    mod_C13_v2 = mean(temperature_model, na.rm = TRUE),
    .groups = "drop"
  )

R11_plot_v2 <- input_model_R11_v2 %>%
  mutate(date = as.Date(datetime)) %>%
  group_by(date, period) %>%
  summarise(
    mod_R11_v2 = mean(temperature_model, na.rm = TRUE),
    .groups = "drop"
  )

###############################################################################
# EXTRA VARIABLES
###############################################################################

C13_extra_v2 <- input_model_C13_v2 %>%
  mutate(date = as.Date(datetime)) %>%
  group_by(date, period) %>%
  summarise(
    LE_C13    = mean(LE_mean, na.rm = TRUE),
    RH_C13    = mean(RH_mean, na.rm = TRUE),
    VPD_C13   = mean(VPD_mean, na.rm = TRUE),
    Fd_C13    = mean(Fd_down_mean, na.rm = TRUE),
    Fb_C13    = mean(Fb_down_mean, na.rm = TRUE),
    Tsurf_C13 = mean(Tsurf_mean, na.rm = TRUE),
    .groups = "drop"
  )

R11_extra_v2 <- input_model_R11_v2 %>%
  mutate(date = as.Date(datetime)) %>%
  group_by(date, period) %>%
  summarise(
    LE_R11    = mean(LE_mean, na.rm = TRUE),
    RH_R11    = mean(RH_mean, na.rm = TRUE),
    VPD_R11   = mean(VPD_mean, na.rm = TRUE),
    Fd_R11    = mean(Fd_down_mean, na.rm = TRUE),
    Fb_R11    = mean(Fb_down_mean, na.rm = TRUE),
    Tsurf_R11 = mean(Tsurf_mean, na.rm = TRUE),
    .groups = "drop"
  )

###############################################################################
# JOIN
###############################################################################

comparison <- C13_plot %>%
  inner_join(R11_plot, by = c("date", "period")) %>%
  inner_join(C13_plot_v2, by = c("date", "period")) %>%
  inner_join(R11_plot_v2, by = c("date", "period")) %>%
  inner_join(C13_extra_v2, by = c("date", "period")) %>%
  inner_join(R11_extra_v2, by = c("date", "period"))

###############################################################################
# REMOVE PERIOD
###############################################################################

exclude_start <- as.Date("2023-06-21")
exclude_end   <- as.Date("2023-06-27")

comparison <- comparison %>%
  mutate(
    across(
      where(is.numeric),
      ~ ifelse(
        date >= exclude_start &
          date <= exclude_end,
        NA,
        .
      )
    )
  )

###############################################################################
# DELTAS
###############################################################################

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

###############################################################################
# SEASON CLASSIFICATION
###############################################################################

comparison <- comparison %>%
  mutate(
    season = ifelse(
      month(date) >= 5 &
        (
          month(date) < 12 |
            (month(date) == 12 & day(date) <= 15)
        ),
      "Wet season",
      "Dry season"
    )
  )

###############################################################################
# WET SEASON BACKGROUND
###############################################################################

years <- unique(year(comparison$date))

wet_season <- data.frame(
  xmin = as.Date(paste0(years, "-05-01")),
  xmax = as.Date(paste0(years, "-12-15")),
  ymin = -Inf,
  ymax = Inf
)

###############################################################################
# SEASON RANGES
###############################################################################

years <- sort(unique(year(comparison$date)))

wet_ranges <- data.frame(
  season = "Wet season",
  xmin = as.Date(paste0(years, "-05-01")),
  xmax = as.Date(paste0(years, "-12-15"))
)

dry_ranges_1 <- data.frame(
  season = "Dry season",
  xmin = as.Date(paste0(years, "-12-16")),
  xmax = as.Date(paste0(years, "-12-31"))
)

dry_ranges_2 <- data.frame(
  season = "Dry season",
  xmin = as.Date(paste0(years, "-01-01")),
  xmax = as.Date(paste0(years, "-04-30"))
)

season_ranges <- bind_rows(
  wet_ranges,
  dry_ranges_1,
  dry_ranges_2
) %>%
  filter(
    xmin <= max(comparison$date) &
      xmax >= min(comparison$date)
  )

###############################################################################
# TIMESERIES FUNCTION
###############################################################################

plot_timeseries <- function(df, period_filter, title_text) {

  df_sub <- df %>%
    filter(period == period_filter)

  season_means <- df_sub %>%
    group_by(season) %>%
    summarise(
      mean_obs    = mean(delta_obs, na.rm = TRUE),
      mean_mod    = mean(delta_mod, na.rm = TRUE),
      mean_mod_v2 = mean(delta_mod_v2, na.rm = TRUE),
      .groups = "drop"
    )

  line_data <- bind_rows(

    season_ranges %>%
      left_join(
        season_means %>%
          select(season, mean_obs),
        by = "season"
      ) %>%
      mutate(
        value = mean_obs,
        type = "Observed"
      ),

    season_ranges %>%
      left_join(
        season_means %>%
          select(season, mean_mod),
        by = "season"
      ) %>%
      mutate(
        value = mean_mod,
        type = "Modelled v1.0"
      ),

    season_ranges %>%
      left_join(
        season_means %>%
          select(season, mean_mod_v2),
        by = "season"
      ) %>%
      mutate(
        value = mean_mod_v2,
        type = "Modelled v1.1"
      )

  )

  label_data <- line_data %>%
    mutate(
      x = ifelse(
        season == "Wet season",
        as.Date("2023-08-15"),
        as.Date("2023-02-15")
      ),

      label = paste0(
        ifelse(season == "Wet season", "Wet", "Dry"),
        " = ",
        round(value, 2),
        " °C"
      )
    )

  p <- ggplot(df_sub, aes(x = date)) +

    geom_rect(
      data = wet_season,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax
      ),
      inherit.aes = FALSE,
      fill = "lightblue",
      alpha = 0.6
    ) +

    geom_smooth(
      method = "gam",
      formula = y ~ s(x, k = 20),
      aes(
        y = delta_obs,
        colour = "Observed"
      ),
      size = 1,
      alpha = 0.8
    ) +

    geom_smooth(
      method = "gam",
      formula = y ~ s(x, k = 20),
      aes(
        y = delta_mod,
        colour = "Modelled v1.0"
      ),
      size = 1,
      alpha = 0.8
    ) +

    geom_smooth(
      method = "gam",
      formula = y ~ s(x, k = 20),
      aes(
        y = delta_mod_v2,
        colour = "Modelled v1.1"
      ),
      size = 1,
      alpha = 0.8
    ) +

    geom_segment(
      data = line_data,
      aes(
        x = xmin,
        xend = xmax,
        y = value,
        yend = value,
        colour = type
      ),
      linewidth = 0.8,
      inherit.aes = FALSE
    ) +

    geom_label(
      data = label_data,
      aes(
        x = x,
        y = value,
        label = label,
        colour = type
      ),
      fill = "white",
      label.size = 0.3,
      label.r = unit(0.15, "lines"),
      label.padding = unit(0.15, "lines"),
      size = 5,
      fontface = "bold",
      inherit.aes = FALSE,
      show.legend = FALSE
    ) +

    scale_colour_manual(
      values = c(
        "Observed" = "blue",
        "Modelled v1.0" = "orange",
        "Modelled v1.1" = "red"
      )
    ) +

    labs(
      title = title_text,
      x = "Date",
      y = "Air temperature difference (R11 - C13, °C)",
      colour = NULL
    ) +

    theme_bw() +

    theme(
      axis.title = element_text(size = 18),
      axis.text = element_text(size = 15),
      legend.text = element_text(size = 15),
      plot.title = element_text(size = 18, face = "bold"),
      legend.position = "top"
    )

  return(p)

}

###############################################################################
# CREATE TIMESERIES FIGURES
###############################################################################

timeseries_day <- plot_timeseries(
  comparison,
  period_filter = "Day",
  title_text = "Daytime (06-18 EST)"
)

timeseries_night <- plot_timeseries(
  comparison,
  period_filter = "Night",
  title_text = "Nighttime (18-06 EST)"
)

###############################################################################
# PROCESS VARIABLES LONG FORMAT
###############################################################################

comparison_long <- comparison %>%
  select(
    date,
    season,
    period,
    delta_LE,
    delta_RH,
    delta_VPD,
    delta_Fd,
    delta_Fb,
    delta_Tsurf
  ) %>%
  pivot_longer(
    cols = -c(date, season, period),
    names_to = "variable",
    values_to = "value"
  )

comparison_long <- comparison_long %>%
  mutate(
    variable = recode(
      variable,
      delta_LE    = "c) Δ LE (W m⁻²)",
      delta_RH    = "e) Δ RH (%)",
      delta_VPD   = "f) Δ VPD (kPa)",
      delta_Fd    = "a) Δ Diffuse SW (W m⁻²)",
      delta_Fb    = "b) Δ Direct SW (W m⁻²)",
      delta_Tsurf = "d) Δ Surface T (°C)"
    )
  )

###############################################################################
# FACET FUNCTION
###############################################################################

plot_facet <- function(df, vars, period_filter, title = NULL) {

  df_sub <- df %>%
    filter(
      variable %in% vars,
      period == period_filter
    )

  means <- df_sub %>%
    group_by(variable, season) %>%
    summarise(
      mean_val = mean(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    merge(season_ranges, by = "season") %>%
    mutate(
      x_pos = ifelse(
        season == "Wet season",
        as.Date("2023-08-15"),
        as.Date("2023-02-15")
      )
    )

  p <- ggplot(
    df_sub,
    aes(
      x = date,
      y = value
    )
  ) +

    geom_rect(
      data = wet_season,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax
      ),
      inherit.aes = FALSE,
      fill = "lightblue",
      alpha = 0.6
    ) +

    geom_smooth(
      method = "gam",
      formula = y ~ s(x, k = 20),
      colour = "red",
      size = 1,
      alpha = 0.8,
      na.rm = TRUE
    ) +

    geom_segment(
      data = means,
      aes(
        x = xmin,
        xend = xmax,
        y = mean_val,
        yend = mean_val
      ),
      colour = "red",
      linewidth = 0.8
    ) +

    geom_label(
      data = means,
      aes(
        x = x_pos,
        y = mean_val,
        label = paste0(
          ifelse(
            season == "Wet season",
            "Wet",
            "Dry"
          ),
          " = ",
          round(mean_val, 2)
        )
      ),
      fill = "white",
      colour = "red",
      label.size = 0.3,
      label.r = unit(0.15, "lines"),
      label.padding = unit(0.15, "lines"),
      size = 5,
      fontface = "bold"
    ) +

    scale_x_date(
      limits = c(
        min(df_sub$date),
        max(df_sub$date) + days(10)
      )
    ) +

    facet_wrap(
      ~variable,
      scales = "free_y",
      ncol = 1
    ) +

    theme_bw() +

    labs(
      title = title,
      x = "Date",
      y = NULL
    ) +

    theme(
      strip.text = element_text(
        size = 18,
        face = "bold"
      ),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 13),
      plot.title = element_text(size = 18, face = "bold"),
      legend.position = "none"
    )

  return(p)

}

###############################################################################
# DAY FACETS
###############################################################################

p_radiation_day <- plot_facet(
  comparison_long,
  vars = c(
    "a) Δ Diffuse SW (W m⁻²)",
    "b) Δ Direct SW (W m⁻²)"
  ),
  period_filter = "Day",
  title = "Daytime (06-18 EST)"
)

p_humidity_day <- plot_facet(
  comparison_long,
  vars = c(
    "e) Δ RH (%)",
    "f) Δ VPD (kPa)"
  ),
  period_filter = "Day",
  title = NULL
)

p_LE_day <- plot_facet(
  comparison_long,
  vars = c(
    "c) Δ LE (W m⁻²)",
    "d) Δ Surface T (°C)"
  ),
  period_filter = "Day",
  title = NULL
)

combined_day <-
  p_radiation_day +
  p_LE_day +
  p_humidity_day +
  plot_layout(ncol = 3)

###############################################################################
# NIGHT FACETS
###############################################################################

p_radiation_night <- plot_facet(
  comparison_long,
  vars = c(
    "a) Δ Diffuse SW (W m⁻²)",
    "b) Δ Direct SW (W m⁻²)"
  ),
  period_filter = "Night",
  title = "Nighttime (18-06 EST)"
)

p_humidity_night <- plot_facet(
  comparison_long,
  vars = c(
    "e) Δ RH (%)",
    "f) Δ VPD (kPa)"
  ),
  period_filter = "Night",
  title = NULL
)

p_LE_night <- plot_facet(
  comparison_long,
  vars = c(
    "c) Δ LE (W m⁻²)",
    "d) Δ Surface T (°C)"
  ),
  period_filter = "Night",
  title = NULL
)

combined_night <-
  p_radiation_night +
  p_LE_night +
  p_humidity_night +
  plot_layout(ncol = 3)

###############################################################################
# PRINT FIGURES
###############################################################################

print(timeseries_day)
print(timeseries_night)

print(combined_day)
print(combined_night)

###############################################################################
# SAVE FIGURES
###############################################################################

ggsave(
  paste0(
    output_path,
    "liana-effect_timeseries_day.png"
  ),
  plot = timeseries_day,
  width = 12,
  height = 6,
  dpi = 500
)

ggsave(
  paste0(
    output_path,
    "liana-effect_timeseries_night.png"
  ),
  plot = timeseries_night,
  width = 12,
  height = 6,
  dpi = 500
)

ggsave(
  paste0(
    output_path,
    "combined_process_facets_day.png"
  ),
  combined_day,
  width = 14,
  height = 8,
  dpi = 300
)

ggsave(
  paste0(
    output_path,
    "combined_process_facets_night.png"
  ),
  combined_night,
  width = 14,
  height = 8,
  dpi = 300
)

###############################################################################
# METRICS
###############################################################################

metrics <- comparison %>%
  group_by(period) %>%
  summarise(
    rmse_v1.1 = sqrt(
      mean(
        (delta_mod_v2 - delta_obs)^2,
        na.rm = TRUE
      )
    ),

    correlation_v1.1 = cor(
      delta_obs,
      delta_mod_v2,
      use = "complete.obs"
    ),
    .groups = "drop"
  )

print(metrics)
