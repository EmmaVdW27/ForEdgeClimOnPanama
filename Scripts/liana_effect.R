###############################################################################
# Liana effect analysis (temperature + process variables)
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
  "Output/timeseries/data/timeseries_C13_30-gens_calibrated.rds"
)

input_model_R11_v2 <- readRDS(
  "Output/timeseries/data/timeseries_R11_30-gens_calibrated.rds"
)

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
  mutate(date = as.Date(datetime)) %>%
  group_by(date) %>%
  summarise(
    obs_C13 = mean(value_mean, na.rm = TRUE),
    mod_C13 = mean(temperature_model, na.rm = TRUE)
  )

R11_plot <- input_model_R11 %>%
  mutate(date = as.Date(datetime)) %>%
  group_by(date) %>%
  summarise(
    obs_R11 = mean(value_mean, na.rm = TRUE),
    mod_R11 = mean(temperature_model, na.rm = TRUE)
  )

C13_plot_v2 <- input_model_C13_v2 %>%
  mutate(date = as.Date(datetime)) %>%
  group_by(date) %>%
  summarise(
    mod_C13_v2 = mean(temperature_model, na.rm = TRUE)
  )

R11_plot_v2 <- input_model_R11_v2 %>%
  mutate(date = as.Date(datetime)) %>%
  group_by(date) %>%
  summarise(
    mod_R11_v2 = mean(temperature_model, na.rm = TRUE)
  )

#########
# EXTRA VARIABLES
#########

C13_extra_v2 <- input_model_C13_v2 %>%
  mutate(date = as.Date(datetime)) %>%
  group_by(date) %>%
  summarise(
    LE_C13    = mean(LE_mean, na.rm = TRUE),
    RH_C13    = mean(RH_mean, na.rm = TRUE),
    VPD_C13   = mean(VPD_mean, na.rm = TRUE),
    Fd_C13    = mean(Fd_down_mean, na.rm = TRUE),
    Fb_C13    = mean(Fb_down_mean, na.rm = TRUE),
    Tsurf_C13 = mean(Tsurf_mean, na.rm = TRUE)
  )

R11_extra_v2 <- input_model_R11_v2 %>%
  mutate(date = as.Date(datetime)) %>%
  group_by(date) %>%
  summarise(
    LE_R11    = mean(LE_mean, na.rm = TRUE),
    RH_R11    = mean(RH_mean, na.rm = TRUE),
    VPD_R11   = mean(VPD_mean, na.rm = TRUE),
    Fd_R11    = mean(Fd_down_mean, na.rm = TRUE),
    Fb_R11    = mean(Fb_down_mean, na.rm = TRUE),
    Tsurf_R11 = mean(Tsurf_mean, na.rm = TRUE)
  )

#########
# JOIN
#########

comparison <- C13_plot %>%
  inner_join(R11_plot, by = "date") %>%
  inner_join(C13_plot_v2, by = "date") %>%
  inner_join(R11_plot_v2, by = "date") %>%
  inner_join(C13_extra_v2, by = "date") %>%
  inner_join(R11_extra_v2, by = "date")

#########
# REMOVE PERIOD
#########

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
# SEASON CLASSIFICATION
#########

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

#########
# WET SEASON BACKGROUND
#########

years <- unique(year(comparison$date))

wet_season <- data.frame(
  xmin = as.Date(paste0(years, "-05-01")),
  xmax = as.Date(paste0(years, "-12-15")),
  ymin = -Inf,
  ymax = Inf
)

###################
# SEASON RANGES
###################

years <- sort(unique(year(comparison$date)))

# Wet season ranges
wet_ranges <- data.frame(
  season = "Wet season",
  xmin = as.Date(paste0(years, "-05-01")),
  xmax = as.Date(paste0(years, "-12-15"))
)

# Dry season ranges
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

###################
# MEANS
###################

season_means <- comparison %>%
  group_by(season) %>%
  summarise(
    mean_obs    = mean(delta_obs, na.rm = TRUE),
    mean_mod    = mean(delta_mod, na.rm = TRUE),
    mean_mod_v2 = mean(delta_mod_v2, na.rm = TRUE),
    .groups = "drop"
  )

###################
# LINE DATA
###################

line_data <- bind_rows(

  season_ranges %>%
    left_join(
      season_means %>%
        select(season, mean_obs),
      by = "season"
    ) %>%
    mutate(
      value = mean_obs,
      type = "Observed",
      colour = "blue"
    ),

  season_ranges %>%
    left_join(
      season_means %>%
        select(season, mean_mod),
      by = "season"
    ) %>%
    mutate(
      value = mean_mod,
      type = "Modelled v1.0",
      colour = "orange"
    ),

  season_ranges %>%
    left_join(
      season_means %>%
        select(season, mean_mod_v2),
      by = "season"
    ) %>%
    mutate(
      value = mean_mod_v2,
      type = "Modelled v1.1",
      colour = "red"
    )

)

###################
# LABEL POSITIONS
###################

label_data <- line_data %>%
  mutate(
    x = ifelse(
      season == "Wet season",
      as.Date("2023-08-15"),
      as.Date("2023-02-15")
    ),

    # Small vertical offsets to avoid overlap
    y = case_when(
      season == "Wet season" &
        type == "Modelled v1.0" ~ value - 0.01,

      season == "Wet season" &
        type == "Modelled v1.1" ~ value + 0.01,

      TRUE ~ value
    ),

    label = paste0(
      ifelse(season == "Wet season", "Wet", "Dry"),
      " = ",
      round(value, 2),
      " °C"
    ),

    linetype = ifelse(
      season == "Wet season",
      "solid",
      "solid"
    )
  )

###################
# TIMESERIES PLOT
###################

timeseries <- ggplot(comparison, aes(x = date)) +

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

  #########
# SEASON LINES
#########

geom_segment(
  data = line_data,
  aes(
    x = xmin,
    xend = xmax,
    y = value,
    yend = value,
    colour = type,
    linetype = season
  ),
  linewidth = 0.8,
  inherit.aes = FALSE,
  show.legend = FALSE
) +

  #########
# LABELS
#########

geom_label(
  data = label_data,
  aes(
    x = x,
    y = y,
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
  guides(
    colour = guide_legend(
      override.aes = list(size = 4)
    ),
    linetype = "none"
  ) +

  scale_linetype_manual(
    values = c(
      "Wet season" = "solid",
      "Dry season" = "solid"
    )
  ) +

  labs(
    x = "Date",
    y = "Air temperature difference (R11 - C13, °C)",
    colour = NULL,
    linetype = NULL
  ) +

  theme_bw() +

  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 15),
    legend.text = element_text(size = 15),
    legend.position = "top",
    legend.key.width = unit(2, "cm")
  )

print(timeseries)

###############################################################################
# FACET PLOTS
###############################################################################

comparison_long <- comparison %>%
  select(
    date,
    season,
    delta_LE,
    delta_RH,
    delta_VPD,
    delta_Fd,
    delta_Fb,
    delta_Tsurf
  ) %>%
  pivot_longer(
    cols = -c(date, season),
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

###################
# FACET FUNCTION
###################

plot_facet <- function(df, vars, title = NULL) {

  df_sub <- df %>%
    filter(variable %in% vars)

  means <- df_sub %>%
    group_by(variable, season) %>%
    summarise(
      mean_val = mean(value, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    left_join(season_ranges, by = "season") %>%
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

    #########
  # SEASON LINES
  #########

  geom_segment(
    data = means,
    aes(
      x = xmin,
      xend = xmax,
      y = mean_val,
      yend = mean_val,
      linetype = season
    ),
    colour = "red",
    linewidth = 0.8
  ) +

    #########
  # LABELS
  #########

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

    scale_linetype_manual(
      values = c(
        "Wet season" = "solid",
        "Dry season" = "solid"
      )
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
      y = NULL,
      linetype = NULL
    ) +

    theme(
      strip.text = element_text(
        size = 18,
        face = "bold"
      ),
      axis.title = element_text(size = 16),
      axis.text = element_text(size = 13),
      plot.title = element_text(size = 18),
      legend.position = "none"
    )

  return(p)

}

###################
# FACET FIGURES
###################

p_radiation <- plot_facet(
  comparison_long,
  vars = c(
    "a) Δ Diffuse SW (W m⁻²)",
    "b) Δ Direct SW (W m⁻²)"
  )
)

p_humidity <- plot_facet(
  comparison_long,
  vars = c(
    "e) Δ RH (%)",
    "f) Δ VPD (kPa)"
  )
)

p_LE <- plot_facet(
  comparison_long,
  vars = c(
    "c) Δ LE (W m⁻²)",
    "d) Δ Surface T (°C)"
  )
)

combined_plot <-
  p_radiation +
  p_LE +
  p_humidity +
  plot_layout(ncol = 3)

print(combined_plot)

###############################################################################
# SAVE FIGURES
###############################################################################

ggsave(
  paste0(
    output_path,
    "liana-effect_timeseries_v1_v2.png"
  ),
  plot = timeseries,
  width = 12,
  height = 6,
  dpi = 500
)

ggsave(
  paste0(
    output_path,
    "combined_process_facets.png"
  ),
  combined_plot,
  width = 14,
  height = 8,
  dpi = 300
)

###############################################################################
# METRICS
###############################################################################

rmse_v2 <- sqrt(
  mean(
    (
      comparison$delta_mod_v2 -
        comparison$delta_obs
    )^2,
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
