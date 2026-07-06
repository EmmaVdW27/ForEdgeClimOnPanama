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

temperature_from_model = "temperature_model" # Tair_mean_1m or temperature_model (only sensor voxels)

output_path <- "Output/liana-effect/"

control <- readRDS("Output/timeseries/timeseries_control_calibrated.rds")
removal <- readRDS("Output/timeseries/timeseries_removal_calibrated.rds")

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

###############################################################################
# AGGREGATION
###############################################################################

control <- control %>%
  clean_datetime() %>%
  add_period() %>%
  mutate(date = as.Date(datetime))

removal <- removal %>%
  clean_datetime() %>%
  add_period() %>%
  mutate(date = as.Date(datetime))

control_summary <-
  control %>%
  group_by(
    date,
    period
  ) %>%
  summarise(

    n_control = n_distinct(plot),

    obs_control =
      mean(value_mean, na.rm = TRUE),

    mod_control =
      mean(.data[[temperature_from_model]], na.rm = TRUE),

    LE_control =
      mean(LE_mean, na.rm = TRUE),

    RH_control =
      mean(RH_mean, na.rm = TRUE),

    VPD_control =
      mean(VPD_mean, na.rm = TRUE),

    Fd_control =
      mean(Fd_down_mean, na.rm = TRUE),

    Fb_control =
      mean(Fb_down_mean, na.rm = TRUE),

    Tsurf_control =
      mean(Tsurf_mean, na.rm = TRUE),

    .groups = "drop"
  )

removal_summary <-
  removal %>%
  group_by(
    date,
    period
  ) %>%
  summarise(

    n_removal = n_distinct(plot),

    obs_removal =
      mean(value_mean, na.rm = TRUE),

    mod_removal =
      mean(.data[[temperature_from_model]], na.rm = TRUE),

    LE_removal =
      mean(LE_mean, na.rm = TRUE),

    RH_removal =
      mean(RH_mean, na.rm = TRUE),

    VPD_removal =
      mean(VPD_mean, na.rm = TRUE),

    Fd_removal =
      mean(Fd_down_mean, na.rm = TRUE),

    Fb_removal =
      mean(Fb_down_mean, na.rm = TRUE),

    Tsurf_removal =
      mean(Tsurf_mean, na.rm = TRUE),

    .groups = "drop"
  )

###############################################################################
# JOIN
###############################################################################

comparison <-
  inner_join(
    control_summary,
    removal_summary,
    by = c(
      "date",
      "period"
    )
  )

###############################################################################
# REMOVE PERIOD
###############################################################################

# exclude_start <- as.Date("2023-06-21")
# exclude_end   <- as.Date("2023-06-27")
#
# comparison <- comparison %>%
#   mutate(
#     across(
#       where(is.numeric),
#       ~ ifelse(
#         date >= exclude_start &
#           date <= exclude_end,
#         NA,
#         .
#       )
#     )
#   )

###############################################################################
# DELTAS
###############################################################################

comparison <-
  comparison %>%
  mutate(

    delta_obs =
      obs_removal -
      obs_control,

    delta_mod =
      mod_removal -
      mod_control,

    delta_LE =
      LE_removal -
      LE_control,

    delta_RH =
      RH_removal -
      RH_control,

    delta_VPD =
      VPD_removal -
      VPD_control,

    delta_Fd =
      Fd_removal -
      Fd_control,

    delta_Fb =
      Fb_removal -
      Fb_control,

    delta_Tsurf =
      Tsurf_removal -
      Tsurf_control
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
# ABSOLUTE AIR TEMPERATURE
###############################################################################

absolute_air <- bind_rows(

  comparison %>%
    transmute(
      date,
      season,
      period,
      treatment = "Control",
      type = "Observed",
      value = obs_control
    ),

  comparison %>%
    transmute(
      date,
      season,
      period,
      treatment = "Control",
      type = "Modelled",
      value = mod_control
    ),

  comparison %>%
    transmute(
      date,
      season,
      period,
      treatment = "Removal",
      type = "Observed",
      value = obs_removal
    ),

  comparison %>%
    transmute(
      date,
      season,
      period,
      treatment = "Removal",
      type = "Modelled",
      value = mod_removal
    )

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
  ) %>%
  mutate(type = "Modelled")

air_long <-
  comparison %>%
  transmute(
    date,
    season,
    period,
    variable = "c) Δ Air T (°C)",
    Observed = delta_obs,
    Modelled = delta_mod
  ) %>%
  pivot_longer(
    c(Observed, Modelled),
    names_to = "type",
    values_to = "value"
  )

comparison_long <- bind_rows(
  comparison_long,
  air_long
)

comparison_long_day <- comparison_long %>%
  mutate(
    variable = recode(
      variable,
      delta_Fd        = "a) Δ Diffuse SW (W m⁻²)",
      delta_Fb        = "b) Δ Direct SW (W m⁻²)",
      delta_Tsurf     = "c) Δ Surface T (°C)",
      "c) Δ Air T (°C)" = "d) Δ Air T (°C)",
      delta_LE        = "e) Δ LE (W m⁻²)",
      delta_RH        = "f) Δ RH (%)",
      delta_VPD       = "g) Δ VPD (kPa)"
    )
  )

comparison_long_night <- comparison_long %>%
  mutate(
    variable = recode(
      variable,
      delta_Tsurf     = "a) Δ Surface T (°C)",
      "c) Δ Air T (°C)" = "b) Δ Air T (°C)",
      delta_LE        = "c) Δ LE (W m⁻²)",
      delta_RH        = "d) Δ RH (%)",
      delta_VPD       = "e) Δ VPD (kPa)"
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
    group_by(variable, season, type) %>%
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
      aes(colour = type),
      method = "gam",
      formula = y ~ s(x, k = 20),
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
        yend = mean_val,
        colour = type
      ),
      linewidth = 0.8
    ) +

    geom_label(
      data = means,
      aes(
        x = x_pos,
        y = mean_val,
        label = paste0(
          ifelse(season == "Wet season", "Wet", "Dry"),
          " = ",
          round(mean_val, 2)
        ),
        colour = type
      ),
      fill = "white",
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
    ) +
    scale_colour_manual(
      values = c(
        Observed = "blue",
        Modelled = "red"
      )
    )

  return(p)

}

###############################################################################
# DAY FACETS
###############################################################################

p_radiation_day <- plot_facet(
  comparison_long_day,
  vars = c(
    "a) Δ Diffuse SW (W m⁻²)",
    "b) Δ Direct SW (W m⁻²)"
  ),
  period_filter = "Day",
  title = "Daytime (06-18 EST)"
)

p_middle_day <- plot_facet(
  comparison_long_day,
  vars = c(
    "c) Δ Surface T (°C)",
    "d) Δ Air T (°C)"
  ),
  period_filter = "Day"
)

p_right_day <- plot_facet(
  comparison_long_day,
  vars = c(
    "e) Δ LE (W m⁻²)",
    "f) Δ RH (%)",
    "g) Δ VPD (kPa)"
  ),
  period_filter = "Day"
)

combined_day <-
  p_radiation_day +
  p_middle_day +
  p_right_day +
  plot_layout(ncol = 3)

###############################################################################
# NIGHT FACETS
###############################################################################

p_middle_night <- plot_facet(
  comparison_long_night,
  vars = c(
    "a) Δ Surface T (°C)",
    "b) Δ Air T (°C)"
  ),
  period_filter = "Night",
  title = "Nighttime (18-06 EST)"
)

p_right_night <- plot_facet(
  comparison_long_night,
  vars = c(
    "c) Δ LE (W m⁻²)",
    "d) Δ RH (%)",
    "e) Δ VPD (kPa)"
  ),
  period_filter = "Night"
)

combined_night <-
  p_middle_night +
  p_right_night +
  plot_layout(ncol = 2)



###############################################################################
# absolute air plot
###############################################################################

plot_absolute_air <- function(df, period_filter, title = NULL){

  df_sub <- df %>%
    filter(period == period_filter)

  means <- df_sub %>%
    group_by(treatment, season, type) %>%
    summarise(
      mean_val = mean(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      xmin = if_else(
        season == "Wet season",
        as.Date("2023-05-01"),
        as.Date("2023-01-01")
      ),
      xmax = if_else(
        season == "Wet season",
        as.Date("2023-12-15"),
        as.Date("2023-04-30")
      ),
      x_pos = if_else(
        season == "Wet season",
        as.Date("2023-08-15"),
        as.Date("2023-02-15")
      )
    )

  ggplot(df_sub,
         aes(date,
             value,
             colour = type)) +

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
      formula = y ~ s(x, k = 30),
      linewidth = 1
    ) +

    geom_segment(
      data = means,
      aes(
        x = xmin,
        xend = xmax,
        y = mean_val,
        yend = mean_val,
        colour = type
      ),
      linewidth = 0.8
    ) +

    geom_label(
      data = means,
      aes(
        x = x_pos,
        y = mean_val,
        label = round(mean_val,1),
        colour = type
      ),
      fill = "white",
      size = 4
    ) +

    facet_wrap(
      ~treatment,
      ncol = 2,
      scales = "fixed"
    ) +

    scale_colour_manual(
      values = c(
        Observed = "blue",
        Modelled = "red"
      )
    ) +

    labs(
      title = title,
      x = "Date",
      y = "Air temperature (°C)",
      colour = NULL
    ) +

    theme_bw() +

    theme(
      strip.text = element_text(size = 18, face = "bold"),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 13),
      plot.title = element_text(
        size = 18,
        face = "bold"
      ),
      legend.position = "none",
      panel.spacing.x = unit(2, "lines")
    ) +
    scale_x_date(
      limits = c(
        min(df_sub$date),
        max(df_sub$date) + days(10)
      )
    )
}

###############################################################################
# CREATE FIGURES
###############################################################################

absolute_day <- plot_absolute_air(
  absolute_air,
  period_filter = "Day",
  title = "Daytime (06–18 EST)"
)

absolute_night <- plot_absolute_air(
  absolute_air,
  period_filter = "Night",
  title = "Nighttime (18–06 EST)"

)

absolute_combined <-
  wrap_plots(
    absolute_day,
    absolute_night,
    ncol = 1
  )

###############################################################################
# SAVE FIGURES
###############################################################################

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
  width = 9.33,
  height = 8,
  dpi = 300
)


ggsave(
  paste0(
    output_path,
    "absolute_air_temperature_combined.png"
  ),
  plot = absolute_combined,
  width = 12,
  height = 10,
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
        (delta_mod - delta_obs)^2,
        na.rm = TRUE
      )
    ),

    correlation_v1.1 = cor(
      delta_obs,
      delta_mod,
      use = "complete.obs"
    ),
    .groups = "drop"
  )

print(metrics)
