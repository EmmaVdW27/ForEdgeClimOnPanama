################################################################################
# ForEdgeClim — Canopy structure and relative light validation
#
# Plots: C8 and R7
# Year: 2023
# Configurations: uncalibrated and calibrated
#
# STRUCTURAL PROXY:
#   Mean voxel density above 1 m within a 3 x 3 neighbourhood:
#   the central voxel column plus eight neighbouring columns.
#
# FIVE COMMON CANOPY STRUCTURE CLASSES:
#   Defined from pooled C8 and R7 sensor locations using quintiles.
#   Identical structural boundaries are applied to both plots.
#
# LIGHT AGGREGATION:
#   Annual MEAN relative light per sensor.
#   Class MEAN across sensors.
#   Error bars: +/- 1 SD between sensor-level annual means.
#
# FIGURE STYLE:
#   Observed: blue
#   Modelled: red
#   Transparent fill with coloured outlines
#   C8 above R7
#   Uncalibrated left, calibrated right
################################################################################

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(readxl)
  library(RANN)
})

# ==============================================================================
# 1. SETTINGS
# ==============================================================================

plots <- c("C8", "R7")

configurations <- c(
  "uncalibrated",
  "calibrated"
)

relative_file <-
  "Output/validation/relative_light_distributions/daily_relative_light.rds"

sensor_file <- "Data/sensor_locations_all_plots.xlsx"

voxel_files <- c(
  C8 = "Data/TLS_scaled_DTM_and_grid_c8_Panama.rds",
  R7 = "Data/TLS_scaled_DTM_and_grid_r7_Panama.rds"
)

output_dir <-
  "Output/validation/canopy_openness_C8_R7_3x3"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

sensor_height <- 1

# Maximum matching distance for diagnostics only.
# Sensors are not automatically excluded.
match_distance_warning <- 2

# Neighbourhood radius in grid cells:
# 1 = central column plus eight neighbours (3 x 3).
neighbourhood_radius <- 1L

# Five common canopy structure classes.
class_labels <- c(
  "Very closed",
  "Closed",
  "Intermediate",
  "Open",
  "Very open"
)

colours <- c(
  "Observed" = "blue",
  "Modelled" = "red"
)

# ==============================================================================
# 2. FUNCTIONS
# ==============================================================================

check_columns <- function(dat, required, label) {

  missing <- setdiff(required, names(dat))

  if (length(missing) > 0) {
    stop(
      label,
      " missing columns: ",
      paste(missing, collapse = ", ")
    )
  }
}

normalise_id <- function(x) {

  x <- trimws(as.character(x))

  numeric_id <- grepl("^[0-9]+(\\.0+)?$", x)

  x[numeric_id] <- as.character(
    as.integer(as.numeric(x[numeric_id]))
  )

  x
}

safe_cor <- function(x, y) {

  valid <- is.finite(x) & is.finite(y)

  x <- x[valid]
  y <- y[valid]

  if (
    length(x) < 3 ||
    sd(x) == 0 ||
    sd(y) == 0
  ) {
    return(NA_real_)
  }

  cor(x, y)
}

# ==============================================================================
# 3. READ DAILY RELATIVE LIGHT DATA
# ==============================================================================

if (!file.exists(relative_file)) {
  stop("Missing file: ", relative_file)
}

daily <- readRDS(relative_file)

check_columns(
  daily,
  c(
    "plot",
    "configuration",
    "date",
    "sensor_id",
    "obs_relative",
    "mod_relative"
  ),
  "Daily relative light"
)

daily <- daily %>%
  mutate(
    plot = as.character(plot),
    configuration = tolower(as.character(configuration)),
    date = as.Date(date),
    sensor_id = normalise_id(sensor_id)
  ) %>%
  filter(
    plot %in% plots,
    configuration %in% configurations,
    format(date, "%Y") == "2023",
    is.finite(obs_relative),
    is.finite(mod_relative)
  )

if (nrow(daily) == 0) {
  stop("No valid daily relative light data for 2023.")
}

# ==============================================================================
# 4. READ SENSOR LOCATIONS
# ==============================================================================

sensors_raw <- read_excel(sensor_file)

check_columns(
  sensors_raw,
  c("x", "y", "SN", "sensor", "plot"),
  "Sensor locations"
)

sensors <- sensors_raw %>%
  mutate(
    plot = case_when(
      tolower(trimws(as.character(plot))) == "cplot08" ~ "C8",
      tolower(trimws(as.character(plot))) == "rplot07" ~ "R7",
      TRUE ~ NA_character_
    ),
    sensor_type = toupper(trimws(as.character(sensor))),
    sensor_id = normalise_id(SN),
    x = as.numeric(x),
    y = as.numeric(y),
    z = sensor_height
  ) %>%
  filter(
    !is.na(plot),
    sensor_type == "LS",
    !grepl("vertical", sensor_id, ignore.case = TRUE),
    grepl("^[0-9]+$", sensor_id),
    is.finite(x),
    is.finite(y)
  ) %>%
  select(plot, sensor_id, x, y, z)

if (anyDuplicated(sensors[c("plot", "sensor_id")])) {
  stop("Duplicate sensor IDs within a plot.")
}

cat("\nNumber of LS sensors:\n")
print(sensors %>% count(plot))

# ==============================================================================
# 5. CALCULATE 3 x 3 STRUCTURAL OPENNESS PROXY
# ==============================================================================

calculate_structure <- function(plot_id) {

  cat("\nProcessing structure:", plot_id, "\n")

  voxel_TLS <- readRDS(voxel_files[[plot_id]])

  if (is.null(voxel_TLS$grid)) {
    stop("No $grid found for ", plot_id)
  }

  grid <- voxel_TLS$grid

  check_columns(
    grid,
    c("X", "Y", "Z", "density"),
    paste("Voxel grid", plot_id)
  )

  grid <- grid %>%
    transmute(
      X = as.numeric(X),
      Y = as.numeric(Y),
      Z = as.numeric(Z),
      density = as.numeric(density)
    )

  if (any(!is.finite(as.matrix(grid)))) {
    stop("Non-finite voxel values in ", plot_id)
  }

  if (any(grid$density < 0 | grid$density > 1)) {
    stop("Density outside [0,1] in ", plot_id)
  }

  if (anyDuplicated(grid[c("X", "Y", "Z")])) {
    stop("Duplicate voxel coordinates in ", plot_id)
  }

  plot_sensors <- sensors %>%
    filter(plot == plot_id)

  # --------------------------------------------------------------------------
  # 5.1. Identify grid spacing and XY column indices
  # --------------------------------------------------------------------------

  x_coords <- sort(unique(grid$X))
  y_coords <- sort(unique(grid$Y))

  if (length(x_coords) < 3 || length(y_coords) < 3) {
    stop("Voxel grid is too small for 3 x 3 analysis.")
  }

  dx <- median(diff(x_coords))
  dy <- median(diff(y_coords))

  cat("\nGrid spacing:\n")
  cat("dx =", dx, "m\n")
  cat("dy =", dy, "m\n")

  # Verify that the XY axes are regular.
  if (
    any(abs(diff(x_coords) - dx) > abs(dx) * 1e-5) ||
    any(abs(diff(y_coords) - dy) > abs(dy) * 1e-5)
  ) {
    stop(
      "Irregular voxel grid spacing in ", plot_id,
      ". The 3 x 3 neighbourhood requires a regular grid."
    )
  }

  grid <- grid %>%
    mutate(
      ix = match(X, x_coords),
      iy = match(Y, y_coords)
    )

  # --------------------------------------------------------------------------
  # 5.2. Match sensors to nearest voxel XY column
  # --------------------------------------------------------------------------

  grid_xy <- grid %>%
    distinct(X, Y, ix, iy)

  nn <- RANN::nn2(
    data = as.matrix(grid_xy[, c("X", "Y")]),
    query = as.matrix(plot_sensors[, c("x", "y")]),
    k = 1
  )

  plot_sensors <- plot_sensors %>%
    mutate(
      voxel_x = grid_xy$X[nn$nn.idx[, 1]],
      voxel_y = grid_xy$Y[nn$nn.idx[, 1]],
      centre_ix = grid_xy$ix[nn$nn.idx[, 1]],
      centre_iy = grid_xy$iy[nn$nn.idx[, 1]],
      match_distance_m = as.numeric(nn$nn.dists[, 1])
    )

  matching <- plot_sensors %>%
    transmute(
      plot,
      sensor_id,
      sensor_x = x,
      sensor_y = y,
      voxel_x,
      voxel_y,
      match_distance_m
    )

  write.csv(
    matching,
    file.path(
      output_dir,
      paste0("sensor_voxel_matching_", plot_id, ".csv")
    ),
    row.names = FALSE
  )

  cat("\nMatching distances:\n")
  print(summary(plot_sensors$match_distance_m))

  if (
    any(plot_sensors$match_distance_m >
        match_distance_warning)
  ) {
    warning(
      plot_id,
      ": Some sensors are more than ",
      match_distance_warning,
      " m from their nearest voxel column."
    )
  }

  # --------------------------------------------------------------------------
  # 5.3. Calculate structural proxy from 3 x 3 columns
  # --------------------------------------------------------------------------

  structure <- lapply(
    seq_len(nrow(plot_sensors)),
    function(i) {

      s <- plot_sensors[i, ]

      # Select the central column and its eight neighbours
      # using grid indices, avoiding floating-point XY issues.

      neighbourhood <- grid %>%
        filter(
          abs(ix - s$centre_ix) <= neighbourhood_radius,
          abs(iy - s$centre_iy) <= neighbourhood_radius,
          Z > s$z
        )

      if (nrow(neighbourhood) == 0) {
        return(NULL)
      }

      n_columns <- neighbourhood %>%
        distinct(ix, iy) %>%
        nrow()

      n_voxels <- nrow(neighbourhood)

      # Verify completeness of the vertical columns.
      # Missing voxels are NOT silently treated as empty.

      n_vertical_levels <- neighbourhood %>%
        distinct(Z) %>%
        nrow()

      expected_voxels <- n_columns * n_vertical_levels

      if (n_voxels != expected_voxels) {
        stop(
          "Incomplete voxel neighbourhood in ",
          plot_id,
          " at sensor ",
          s$sensor_id,
          ". Expected ",
          expected_voxels,
          " voxels, found ",
          n_voxels,
          ". Check whether zero-density voxels are present."
        )
      }

      # Mean density across all voxels in the neighbourhood.
      # The gap proxy is one minus mean density.

      mean_density <- mean(neighbourhood$density)

      data.frame(
        plot = plot_id,
        sensor_id = s$sensor_id,
        sensor_x = s$x,
        sensor_y = s$y,
        voxel_x = s$voxel_x,
        voxel_y = s$voxel_y,
        match_distance_m = s$match_distance_m,
        n_columns = n_columns,
        n_vertical_levels = n_vertical_levels,
        n_voxels = n_voxels,
        mean_density = mean_density,
        column_gap_proxy = 1 - mean_density
      )
    }
  )

  bind_rows(structure)
}

structure <- bind_rows(
  lapply(plots, calculate_structure)
)

if (nrow(structure) == 0) {
  stop("No valid structural data.")
}

# ==============================================================================
# 6. CHECK 3 x 3 NEIGHBOURHOOD COMPLETENESS
# ==============================================================================

cat("\nNumber of columns per sensor:\n")

print(
  table(
    structure$plot,
    structure$n_columns
  )
)

cat("\nNumber of voxels per sensor neighbourhood:\n")

print(
  structure %>%
    group_by(plot) %>%
    summarise(
      min_voxels = min(n_voxels),
      max_voxels = max(n_voxels),
      .groups = "drop"
    )
)

if (any(structure$n_columns < 9)) {
  warning(
    "Some sensors have fewer than nine voxel columns. ",
    "Check their positions relative to the grid boundaries."
  )
}

# ==============================================================================
# 7. FIVE COMMON CANOPY STRUCTURE CLASSES
# ==============================================================================

# Only sensors represented in the light validation are included.

structure <- structure %>%
  semi_join(
    daily %>% distinct(plot, sensor_id),
    by = c("plot", "sensor_id")
  )

# Common quintile boundaries calculated from BOTH plots together.

class_breaks <- as.numeric(
  quantile(
    structure$column_gap_proxy,
    probs = seq(0, 1, 0.2),
    na.rm = TRUE,
    names = FALSE
  )
)

if (anyDuplicated(class_breaks)) {
  stop(
    "Identical quintile boundaries were found. ",
    "Five distinct structural classes cannot be defined ",
    "using these pooled quantiles."
  )
}

class_table <- data.frame(
  canopy_class = class_labels,
  lower_bound = head(class_breaks, -1),
  upper_bound = tail(class_breaks, -1)
)

cat("\nCommon canopy structure class boundaries:\n")
print(class_table)

write.csv(
  class_table,
  file.path(output_dir, "canopy_class_boundaries.csv"),
  row.names = FALSE
)

structure <- structure %>%
  mutate(
    canopy_class = cut(
      column_gap_proxy,
      breaks = class_breaks,
      labels = class_labels,
      include.lowest = TRUE
    ),
    canopy_class = factor(
      canopy_class,
      levels = class_labels
    )
  )

class_inventory <- structure %>%
  count(
    plot,
    canopy_class,
    name = "n_sensors",
    .drop = FALSE
  )

cat("\nSensor counts per class:\n")
print(class_inventory)

write.csv(
  class_inventory,
  file.path(output_dir, "canopy_class_inventory.csv"),
  row.names = FALSE
)

write.csv(
  structure,
  file.path(output_dir, "sensor_canopy_structure.csv"),
  row.names = FALSE
)

# ==============================================================================
# 8. JOIN STRUCTURE AND RELATIVE LIGHT
# ==============================================================================

missing_sensors <- daily %>%
  distinct(plot, sensor_id) %>%
  anti_join(
    structure %>% distinct(plot, sensor_id),
    by = c("plot", "sensor_id")
  )

write.csv(
  missing_sensors,
  file.path(output_dir, "unmatched_sensors.csv"),
  row.names = FALSE
)

joined <- daily %>%
  inner_join(
    structure,
    by = c("plot", "sensor_id")
  )

if (nrow(joined) == 0) {
  stop("No matching sensors between structure and light data.")
}

# Re-normalise using the retained sensor subset.
# When all original sensors are retained, this has no effect.

joined <- joined %>%
  group_by(plot, configuration, date) %>%
  mutate(
    obs_relative = obs_relative / mean(obs_relative),
    mod_relative = mod_relative / mean(mod_relative)
  ) %>%
  ungroup()

write.csv(
  joined,
  file.path(output_dir, "daily_sensor_relative_light.csv"),
  row.names = FALSE
)

# ==============================================================================
# 9. ANNUAL MEAN RELATIVE LIGHT PER SENSOR
# ==============================================================================

annual <- joined %>%
  group_by(
    plot,
    configuration,
    sensor_id,
    column_gap_proxy,
    canopy_class
  ) %>%
  summarise(
    obs = mean(obs_relative),
    mod = mean(mod_relative),
    n_days = n_distinct(date),
    .groups = "drop"
  )

write.csv(
  annual,
  file.path(output_dir, "annual_sensor_means.csv"),
  row.names = FALSE
)

# ==============================================================================
# 10. DAILY CLASS-LEVEL LIGHT
# ==============================================================================

daily_class <- joined %>%
  group_by(
    plot,
    configuration,
    date,
    canopy_class
  ) %>%
  summarise(
    obs = mean(obs_relative),
    mod = mean(mod_relative),
    n_sensors = n_distinct(sensor_id),
    .groups = "drop"
  )

write.csv(
  daily_class,
  file.path(output_dir, "daily_class_relative_light.csv"),
  row.names = FALSE
)

# ==============================================================================
# 11. CLASS-WISE MODEL PERFORMANCE
# ==============================================================================

class_performance <- daily_class %>%
  mutate(
    error = mod - obs
  ) %>%
  group_by(
    plot,
    configuration,
    canopy_class
  ) %>%
  summarise(
    n_days = n_distinct(date),
    mean_obs = mean(obs),
    mean_mod = mean(mod),
    mean_bias = mean(error),
    rmse = sqrt(mean(error^2)),
    .groups = "drop"
  )

write.csv(
  class_performance,
  file.path(output_dir, "class_performance_summary.csv"),
  row.names = FALSE
)

cat("\nClass-wise performance:\n")
print(class_performance)

# ==============================================================================
# 12. OPEN VERSUS CLOSED CONTRAST
# ==============================================================================

contrasts_wide <- daily_class %>%
  filter(
    canopy_class %in% c("Very closed", "Very open")
  ) %>%
  mutate(
    class_key = gsub(
      " ",
      "_",
      as.character(canopy_class)
    )
  ) %>%
  select(
    plot,
    configuration,
    date,
    class_key,
    obs,
    mod
  ) %>%
  pivot_wider(
    names_from = class_key,
    values_from = c(obs, mod),
    names_sep = "_"
  )

needed <- c(
  "obs_Very_open",
  "obs_Very_closed",
  "mod_Very_open",
  "mod_Very_closed"
)

for (nm in setdiff(needed, names(contrasts_wide))) {
  contrasts_wide[[nm]] <- NA_real_
}

contrasts <- contrasts_wide %>%
  filter(
    if_all(all_of(needed), is.finite)
  ) %>%
  mutate(
    obs_contrast = obs_Very_open - obs_Very_closed,
    mod_contrast = mod_Very_open - mod_Very_closed,
    contrast_bias = mod_contrast - obs_contrast
  )

write.csv(
  contrasts,
  file.path(output_dir, "daily_extreme_class_contrast.csv"),
  row.names = FALSE
)

if (nrow(contrasts) > 0) {

  contrast_summary <- contrasts %>%
    group_by(plot, configuration) %>%
    summarise(
      n_days = n(),
      median_obs_contrast = median(obs_contrast),
      median_mod_contrast = median(mod_contrast),
      median_contrast_bias = median(contrast_bias),
      fraction_correct_sign =
        mean(sign(obs_contrast) == sign(mod_contrast)),
      .groups = "drop"
    )

} else {

  contrast_summary <- data.frame(
    plot = character(),
    configuration = character(),
    n_days = integer(),
    median_obs_contrast = numeric(),
    median_mod_contrast = numeric(),
    median_contrast_bias = numeric(),
    fraction_correct_sign = numeric()
  )

  message(
    "No days contain both extreme classes. ",
    "Contrast summary will be empty."
  )
}

write.csv(
  contrast_summary,
  file.path(output_dir, "canopy_contrast_summary.csv"),
  row.names = FALSE
)

# ==============================================================================
# 13. DAILY SPATIAL CORRELATIONS
# ==============================================================================

spatial <- joined %>%
  group_by(plot, configuration, date) %>%
  summarise(
    n_sensors = n_distinct(sensor_id),
    spatial_pearson = safe_cor(
      obs_relative,
      mod_relative
    ),
    .groups = "drop"
  )

write.csv(
  spatial,
  file.path(output_dir, "daily_spatial_correlations.csv"),
  row.names = FALSE
)

# ==============================================================================
# 14. FIGURE SETTINGS
# ==============================================================================

base_theme <- theme_bw(base_size = 17) +
  theme(
    strip.text = element_text(size = 18),
    axis.title = element_text(size = 20),
    axis.text = element_text(size = 15),
    legend.position = "bottom",
    legend.text = element_text(size = 17),
    panel.grid.minor = element_blank()
  )

# ==============================================================================
# 15. FIGURE 1 — ANNUAL MEAN LIGHT VS STRUCTURAL OPENNESS
# ==============================================================================

annual_long <- annual %>%
  pivot_longer(
    cols = c(obs, mod),
    names_to = "source",
    values_to = "relative_light"
  ) %>%
  mutate(
    source = recode(
      source,
      obs = "Observed",
      mod = "Modelled"
    ),
    source = factor(
      source,
      levels = c("Observed", "Modelled")
    ),
    plot = factor(
      plot,
      levels = plots
    ),
    configuration = factor(
      configuration,
      levels = configurations
    )
  )

p1 <- ggplot(
  annual_long,
  aes(
    x = column_gap_proxy,
    y = relative_light,
    colour = source
  )
) +
  geom_point(
    alpha = 0.65,
    size = 2.5
  ) +
  geom_smooth(
    method = "lm",
    se = TRUE,
    linewidth = 1
  ) +
  facet_grid(plot ~ configuration) +
  scale_colour_manual(values = colours) +
  labs(
    x = "Structural openness proxy (3 x 3 voxel columns)",
    y = "Annual mean relative light availability",
    colour = NULL
  ) +
  base_theme

ggsave(
  file.path(output_dir, "annual_light_vs_column_gap.png"),
  p1,
  width = 13,
  height = 9,
  dpi = 300
)

# ==============================================================================
# 16. FIGURE 2 — ANNUAL MEAN LIGHT BY CANOPY CLASS
# ==============================================================================

# Annual mean per sensor, followed by class mean.
# Error bars: +/- 1 SD across annual sensor means.

annual_class <- annual_long %>%
  group_by(
    plot,
    configuration,
    canopy_class,
    source
  ) %>%
  summarise(
    mean_light = mean(relative_light),
    sd_light = if (n() > 1) {
      sd(relative_light)
    } else {
      NA_real_
    },
    n_sensors = n(),
    .groups = "drop"
  )

write.csv(
  annual_class,
  file.path(output_dir, "annual_class_light_summary.csv"),
  row.names = FALSE
)

dodge <- position_dodge(width = 0.8)

p2 <- ggplot(
  annual_class,
  aes(
    x = canopy_class,
    y = mean_light,
    fill = source,
    colour = source,
    group = source
  )
) +
  geom_col(
    position = dodge,
    width = 0.7,
    linewidth = 1.1,
    alpha = 0.20
  ) +
  geom_errorbar(
    aes(
      ymin = pmax(0, mean_light - sd_light),
      ymax = mean_light + sd_light
    ),
    position = dodge,
    width = 0.18,
    linewidth = 0.8,
    na.rm = TRUE
  ) +
  facet_grid(plot ~ configuration) +
  scale_fill_manual(
    values = colours,
    breaks = c("Observed", "Modelled")
  ) +
  scale_colour_manual(
    values = colours,
    breaks = c("Observed", "Modelled")
  ) +
  scale_x_discrete(
    drop = FALSE,
    labels = c(
      "Very\nclosed",
      "Closed",
      "Intermediate",
      "Open",
      "Very\nopen"
    )
  ) +
  labs(
    title = "Relative understorey light availability across canopy structure classes",
    x = "Canopy structure class",
    y = "Annual mean relative light availability",
    fill = NULL,
    colour = NULL
  ) +
  guides(
    colour = "none",
    fill = guide_legend(
      override.aes = list(
        alpha = 0.20,
        linewidth = 1.0,
        colour = c("blue", "red")
      )
    )
  ) +
  base_theme +
  theme(
    # Main title
    plot.title = element_text(
      size = 26,
      face = "bold",
      hjust = 0,
      margin = margin(b = 18)
    ),

    # Axis titles
    axis.title.x = element_text(
      size = 25,
      margin = margin(t = 15)
    ),
    axis.title.y = element_text(
      size = 25,
      margin = margin(r = 15)
    ),

    # Axis tick labels
    axis.text.x = element_text(
      size = 19,
      angle = 0,
      hjust = 0.5,
      colour = "black"
    ),
    axis.text.y = element_text(
      size = 20,
      colour = "black"
    ),

    # Facet labels
    strip.text = element_text(
      size = 23,
      face = "plain",
      colour = "black"
    ),

    # Legend
    legend.position = "bottom",
    legend.text = element_text(size = 22),
    legend.key.width = grid::unit(1.1, "cm"),
    legend.key.height = grid::unit(1.1, "cm"),
    legend.spacing.x = grid::unit(0.5, "cm"),

    # Additional spacing
    panel.spacing = grid::unit(0.9, "cm"),
    plot.margin = margin(20, 25, 20, 20)
  )

ggsave(
  file.path(output_dir, "relative_light_by_canopy_class.png"),
  p2,
  width = 16,
  height = 10,
  dpi = 300
)

# ==============================================================================
# 17. FIGURE 3 — DAILY OPEN-CLOSED CONTRAST
# ==============================================================================

if (nrow(contrasts) > 0) {

  contrast_long <- contrasts %>%
    pivot_longer(
      cols = c(obs_contrast, mod_contrast),
      names_to = "source",
      values_to = "contrast"
    ) %>%
    mutate(
      source = recode(
        source,
        obs_contrast = "Observed",
        mod_contrast = "Modelled"
      ),
      source = factor(
        source,
        levels = c("Observed", "Modelled")
      ),
      plot = factor(
        plot,
        levels = plots
      ),
      configuration = factor(
        configuration,
        levels = configurations
      )
    )

  p3 <- ggplot(
    contrast_long,
    aes(
      x = date,
      y = contrast,
      colour = source
    )
  ) +
    geom_line(
      alpha = 0.7,
      linewidth = 0.7
    ) +
    facet_grid(plot ~ configuration) +
    scale_colour_manual(values = colours) +
    labs(
      x = "Date",
      y = "Very open - very closed relative light",
      colour = NULL
    ) +
    base_theme

  ggsave(
    file.path(output_dir, "daily_extreme_class_contrast.png"),
    p3,
    width = 13,
    height = 9,
    dpi = 300
  )
}

# ==============================================================================
# 18. COMPLETION
# ==============================================================================

cat("\n========================================\n")
cat("CANOPY STRUCTURE VALIDATION COMPLETED\n")
cat("========================================\n")

cat("\nOutput directory:\n")
cat(output_dir, "\n")

cat("\nCommon class boundaries:\n")
print(class_table)

cat("\nSensor distribution across classes:\n")
print(class_inventory)

cat("\nNeighbourhood column counts:\n")
print(table(structure$plot, structure$n_columns))

cat("\nIMPORTANT:\n")
cat(
  "Structural openness is based on mean voxel density ",
  "in a 3 x 3 neighbourhood above sensor height.\n",
  "It is not a true hemispherical canopy gap fraction.\n",
  "Annual MEAN relative light availability is retained.\n",
  "Error bars show between-sensor SD.\n",
  "Inspect voxel completeness and sensor matching diagnostics.\n",
  sep = ""
)
