################################################################################
#
# ForEdgeClim — Relative spatial light distribution validation
#
# Compares observed PAR and modelled SW-derived PAR at sensor locations,
# focusing exclusively on relative spatial light distributions.
#
# Plots:
#   C8 (liana-present)
#   R7 (liana-removal)
#
# Configurations:
#   Calibrated
#   Uncalibrated
#
# Analysis:
#   - 08:00–16:00 Panama local time (13:00–21:00 UTC)
#   - Daily mean light availability per sensor
#   - Normalisation by the daily spatial mean
#   - Comparison of relative spatial distributions
#   - Wasserstein distance, P10, P50, P90 and heterogeneity
#   - Paired comparison of calibrated and uncalibrated simulations
#
# No absolute PAR performance metrics are calculated.
#
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
})

# ==============================================================================
# 1. SETTINGS
# ==============================================================================

input_dir <- "Output/timeseries/data"
output_dir <- "Output/validation/relative_light_distributions"

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

plots <- c("C8", "R7")
configurations <- c("uncalibrated", "calibrated")

# Minimum number of valid hourly observations per sensor-day
min_hours <- 6

# Minimum number of sensors per plot-day
min_sensors <- 10

# Avoid normalisation of days with negligible observed light.
# Used only as a quality-control threshold, not as a validation metric.
min_mean_PAR <- 5

# ==============================================================================
# 2. FUNCTIONS
# ==============================================================================

# Exact empirical 1D Wasserstein distance.
# Returns a dimensionless distance for normalised light distributions.

wasserstein_1d <- function(x, y) {

  x <- sort(x[is.finite(x)])
  y <- sort(y[is.finite(y)])

  if (length(x) < 2 || length(y) < 2) {
    return(NA_real_)
  }

  z <- sort(unique(c(x, y)))

  if (length(z) < 2) {
    return(0)
  }

  mid <- (head(z, -1) + tail(z, -1)) / 2

  Fx <- findInterval(mid, x) / length(x)
  Fy <- findInterval(mid, y) / length(y)

  sum(abs(Fx - Fy) * diff(z))
}

quant <- function(x, p) {
  as.numeric(quantile(x, p, na.rm = TRUE))
}

# ==============================================================================
# 3. READ HOURLY TIMESERIES
# ==============================================================================

read_configuration <- function(plot_id, config) {

  filename <- file.path(
    input_dir,
    paste0("timeseries_", plot_id, "_", config, ".rds")
  )

  if (!file.exists(filename)) {
    stop("Missing file: ", filename)
  }

  dat <- readRDS(filename)

  required <- c(
    "datetime",
    "location_ID_number",
    "radiation_mean",
    "PAR_model"
  )

  missing <- setdiff(required, names(dat))

  if (length(missing) > 0) {
    stop(
      "Missing columns in ", filename, ": ",
      paste(missing, collapse = ", ")
    )
  }

  dat %>%
    transmute(
      plot = plot_id,
      configuration = config,
      datetime = as.POSIXct(datetime, tz = "UTC"),
      sensor_id = as.character(location_ID_number),
      obs = as.numeric(radiation_mean),

      # The modelled PAR variable is derived from total SW.
      # Its absolute magnitude is not evaluated.
      mod = as.numeric(PAR_model)
    ) %>%
    mutate(
      # Panama uses UTC-5 throughout the year.
      local_datetime = datetime - 5 * 3600,
      date = as.Date(local_datetime, tz = "UTC"),
      hour = as.integer(
        format(local_datetime, "%H", tz = "UTC")
      )
    ) %>%
    filter(
      hour >= 8,
      hour <= 16,
      is.finite(obs),
      is.finite(mod),
      obs >= 0,
      mod >= 0
    ) %>%
    group_by(plot, configuration, date, sensor_id) %>%
    summarise(
      obs = mean(obs),
      mod = mean(mod),
      n_hours = n_distinct(datetime),
      .groups = "drop"
    ) %>%
    filter(n_hours >= min_hours)
}

# ==============================================================================
# 4. COMBINE CONFIGURATIONS
# ==============================================================================

daily_raw <- bind_rows(
  lapply(plots, function(p) {
    bind_rows(
      lapply(configurations, function(cfg) {
        read_configuration(p, cfg)
      })
    )
  })
)

if (nrow(daily_raw) == 0) {
  stop("No valid observations available.")
}

# ==============================================================================
# 5. ENSURE IDENTICAL SENSOR-DAYS
# ==============================================================================

# Only sensor-days present in BOTH configurations are retained.
# This makes calibrated vs uncalibrated comparisons directly comparable.

daily_common <- daily_raw %>%
  group_by(plot, date, sensor_id) %>%
  filter(
    n_distinct(configuration) == length(configurations)
  ) %>%
  ungroup()

# Count common sensors per day.
valid_days <- daily_common %>%
  group_by(plot, date) %>%
  summarise(
    n_sensors = n_distinct(sensor_id),
    .groups = "drop"
  ) %>%
  filter(n_sensors >= min_sensors)

daily_common <- daily_common %>%
  semi_join(
    valid_days,
    by = c("plot", "date")
  )

# ==============================================================================
# 6. NORMALISE RELATIVE LIGHT AVAILABILITY
# ==============================================================================

# For each plot, date and configuration:
#
# relative_light = sensor daily mean / spatial daily mean
#
# The mean relative light availability across sensors is 1.
#
# Observations and model predictions are normalised independently.
# Therefore, only their relative spatial distributions are compared.

daily_relative <- daily_common %>%
  group_by(plot, configuration, date) %>%
  mutate(
    obs_spatial_mean = mean(obs),
    mod_spatial_mean = mean(mod)
  ) %>%
  ungroup() %>%
  filter(
    obs_spatial_mean >= min_mean_PAR,
    mod_spatial_mean > 0
  ) %>%
  mutate(
    obs_relative = obs / obs_spatial_mean,
    mod_relative = mod / mod_spatial_mean
  )

# Ensure the same days are retained for both configurations.
common_dates <- daily_relative %>%
  distinct(plot, date, configuration) %>%
  count(plot, date, name = "n_configurations") %>%
  filter(
    n_configurations == length(configurations)
  ) %>%
  select(plot, date)

daily_relative <- daily_relative %>%
  semi_join(
    common_dates,
    by = c("plot", "date")
  )

if (nrow(daily_relative) == 0) {
  stop("No valid days remain after normalisation.")
}

# ==============================================================================
# 7. DAILY RELATIVE DISTRIBUTION METRICS
# ==============================================================================

daily_metrics <- daily_relative %>%
  group_by(plot, configuration, date) %>%
  summarise(
    n_sensors = n(),

    obs_p10 = quant(obs_relative, 0.10),
    mod_p10 = quant(mod_relative, 0.10),

    obs_p50 = median(obs_relative),
    mod_p50 = median(mod_relative),

    obs_p90 = quant(obs_relative, 0.90),
    mod_p90 = quant(mod_relative, 0.90),

    obs_sd = sd(obs_relative),
    mod_sd = sd(mod_relative),

    wasserstein = wasserstein_1d(
      obs_relative,
      mod_relative
    ),

    .groups = "drop"
  ) %>%
  mutate(
    obs_heterogeneity = obs_p90 - obs_p10,
    mod_heterogeneity = mod_p90 - mod_p10,

    p10_bias = mod_p10 - obs_p10,
    p50_bias = mod_p50 - obs_p50,
    p90_bias = mod_p90 - obs_p90,

    heterogeneity_bias =
      mod_heterogeneity - obs_heterogeneity
  )

write.csv(
  daily_metrics,
  file.path(output_dir, "daily_relative_metrics.csv"),
  row.names = FALSE
)

# ==============================================================================
# 8. SUMMARY STATISTICS
# ==============================================================================

summary_metrics <- daily_metrics %>%
  group_by(plot, configuration) %>%
  summarise(
    n_days = n(),
    mean_sensors = mean(n_sensors),

    median_wasserstein = median(wasserstein),

    median_p10_bias = median(p10_bias),
    median_p50_bias = median(p50_bias),
    median_p90_bias = median(p90_bias),

    median_heterogeneity_bias =
      median(heterogeneity_bias),

    .groups = "drop"
  )

print(summary_metrics, width = Inf)

write.csv(
  summary_metrics,
  file.path(output_dir, "relative_distribution_summary.csv"),
  row.names = FALSE
)

# ==============================================================================
# 9. CALIBRATION IMPROVEMENT
# ==============================================================================

# Positive improvement means calibrated distributions
# are closer to observations than uncalibrated distributions.

calibration_comparison <- daily_metrics %>%
  select(
    plot,
    date,
    configuration,
    wasserstein
  ) %>%
  pivot_wider(
    names_from = configuration,
    values_from = wasserstein
  ) %>%
  mutate(
    improvement = uncalibrated - calibrated
  )

calibration_summary <- calibration_comparison %>%
  group_by(plot) %>%
  summarise(
    n_days = n(),
    median_improvement = median(improvement),
    mean_improvement = mean(improvement),
    fraction_days_improved = mean(improvement > 0),
    .groups = "drop"
  )

print(calibration_summary)

write.csv(
  calibration_comparison,
  file.path(output_dir, "daily_calibration_comparison.csv"),
  row.names = FALSE
)

write.csv(
  calibration_summary,
  file.path(output_dir, "calibration_improvement_summary.csv"),
  row.names = FALSE
)

# ==============================================================================
# 10. PREPARE PLOTTING DATA
# ==============================================================================

plot_data <- daily_relative %>%
  select(
    plot,
    configuration,
    date,
    sensor_id,
    obs_relative,
    mod_relative
  ) %>%
  pivot_longer(
    cols = c(obs_relative, mod_relative),
    names_to = "source",
    values_to = "relative_light"
  ) %>%
  mutate(
    source = recode(
      source,
      obs_relative = "Observed",
      mod_relative = "Modelled"
    ),
    source = factor(
      source,
      levels = c("Observed", "Modelled")
    ),
    configuration = factor(
      configuration,
      levels = c("uncalibrated", "calibrated")
    )
  )

# ==============================================================================
# COLOURS AND FIGURE THEME
# ==============================================================================

# Observed = blue
# Modelled = red

colours <- c(
  "Observed" = "blue",
  "Modelled" = "red"
)

# Consistent figure styling across all plots.

base_theme <- theme_bw(base_size = 17) +
  theme(
    plot.title = element_text(
      size = 21,
      face = "bold",
      margin = margin(b = 14)
    ),

    plot.subtitle = element_blank(),

    axis.title.x = element_text(
      size = 20,
      margin = margin(t = 12)
    ),

    axis.title.y = element_text(
      size = 20,
      margin = margin(r = 12)
    ),

    axis.text.x = element_text(size = 16),
    axis.text.y = element_text(size = 16),

    strip.text = element_text(
      size = 18
    ),

    strip.text.x = element_text(
      size = 18,
      margin = margin(t = 8, b = 8)
    ),

    strip.text.y = element_text(
      size = 18,
      margin = margin(l = 8, r = 8)
    ),

    legend.position = "bottom",

    legend.text = element_text(size = 17),
    legend.title = element_text(size = 17),

    legend.key.size = grid::unit(1.2, "cm"),

    panel.grid.minor = element_blank(),

    plot.margin = margin(18, 22, 18, 22)
  )

# ==============================================================================
# 11. OVERALL RELATIVE DISTRIBUTIONS — DENSITY
# ==============================================================================

p_density <- ggplot(
  plot_data,
  aes(
    x = relative_light,
    colour = source,
    fill = source
  )
) +
  geom_density(
    alpha = 0.20,
    linewidth = 1.1,
    adjust = 1.2
  ) +
  facet_grid(plot ~ configuration) +
  scale_colour_manual(values = colours) +
  scale_fill_manual(values = colours) +
  labs(
    title = "Relative spatial light distributions",
    x = "Relative light availability (spatial mean = 1)",
    y = "Density",
    colour = NULL,
    fill = NULL
  ) +
  base_theme

ggsave(
  file.path(output_dir, "relative_light_density.png"),
  p_density,
  width = 12,
  height = 8,
  dpi = 300
)

# ==============================================================================
# 12. OVERALL RELATIVE DISTRIBUTIONS — ECDF
# ==============================================================================

p_ecdf <- ggplot(
  plot_data,
  aes(
    x = relative_light,
    colour = source,
    linetype = source
  )
) +
  stat_ecdf(linewidth = 1.1) +
  facet_grid(plot ~ configuration) +
  scale_colour_manual(values = colours) +
  labs(
    title = "Relative spatial light distributions",
    x = "Relative light availability (spatial mean = 1)",
    y = "Cumulative proportion",
    colour = NULL,
    linetype = NULL
  ) +
  base_theme

ggsave(
  file.path(output_dir, "relative_light_ECDF.png"),
  p_ecdf,
  width = 12,
  height = 8,
  dpi = 300
)

# ==============================================================================
# 13. DAILY WASSERSTEIN DISTANCE
# ==============================================================================

p_wasserstein <- ggplot(
  daily_metrics,
  aes(
    x = date,
    y = wasserstein,
    colour = configuration
  )
) +
  geom_line(linewidth = 0.7, alpha = 0.8) +
  geom_point(size = 1.5, alpha = 0.8) +
  facet_wrap(~plot, ncol = 1) +
  labs(
    title = "Daily relative distribution mismatch",
    x = "Date",
    y = "Wasserstein distance (dimensionless)",
    colour = "Configuration"
  ) +
  base_theme

ggsave(
  file.path(output_dir, "relative_light_wasserstein.png"),
  p_wasserstein,
  width = 12,
  height = 7,
  dpi = 300
)

# ==============================================================================
# 14. DAILY QUANTILE COMPARISON
# ==============================================================================

quantile_data <- daily_metrics %>%
  select(
    plot,
    configuration,
    date,
    obs_p10, mod_p10,
    obs_p50, mod_p50,
    obs_p90, mod_p90
  ) %>%
  pivot_longer(
    cols = matches("^(obs|mod)_p"),
    names_to = c("source", "quantile"),
    names_pattern = "(obs|mod)_(p[0-9]+)",
    values_to = "relative_light"
  ) %>%
  mutate(
    source = recode(
      source,
      obs = "Observed",
      mod = "Modelled"
    ),
    quantile = factor(
      quantile,
      levels = c("p10", "p50", "p90"),
      labels = c("P10", "P50", "P90")
    )
  )

p_quantiles <- ggplot(
  quantile_data,
  aes(
    x = date,
    y = relative_light,
    colour = source
  )
) +
  geom_line(alpha = 0.55) +
  geom_point(size = 1.5, alpha = 0.7) +
  facet_grid(
    quantile ~ plot + configuration,
    scales = "free_y"
  ) +
  scale_colour_manual(values = colours) +
  labs(
    title = "Daily quantiles of relative light availability",
    x = "Date",
    y = "Relative light availability",
    colour = NULL
  ) +
  base_theme

ggsave(
  file.path(output_dir, "relative_light_quantiles.png"),
  p_quantiles,
  width = 16,
  height = 10,
  dpi = 300
)

# ==============================================================================
# 15. CALIBRATION IMPROVEMENT FIGURE
# ==============================================================================

p_improvement <- ggplot(
  calibration_comparison,
  aes(
    x = date,
    y = improvement
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    colour = "grey40"
  ) +
  geom_col(
    fill = "#0072B2",
    alpha = 0.8
  ) +
  facet_wrap(~plot, ncol = 1) +
  labs(
    title = "Improvement in relative light distribution",
    x = "Date",
    y = "Reduction in Wasserstein distance"
  ) +
  base_theme

ggsave(
  file.path(output_dir, "calibration_improvement.png"),
  p_improvement,
  width = 12,
  height = 7,
  dpi = 300
)

# ==============================================================================
# 16. DATA AVAILABILITY
# ==============================================================================

availability <- daily_relative %>%
  group_by(plot, configuration) %>%
  summarise(
    n_days = n_distinct(date),
    first_date = min(date),
    last_date = max(date),
    mean_sensors = n() / n_distinct(date),
    .groups = "drop"
  )

print(availability)

write.csv(
  availability,
  file.path(output_dir, "data_availability.csv"),
  row.names = FALSE
)

# ==============================================================================
# 17. SAVE PROCESSED DATA
# ==============================================================================

saveRDS(
  daily_relative,
  file.path(output_dir, "daily_relative_light.rds")
)

cat("\nRelative light distribution validation completed.\n")
cat("Output directory:", output_dir, "\n")
