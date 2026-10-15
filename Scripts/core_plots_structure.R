###############################################################################
# Forest structural profile figure
#
# Creates vertical forest structural profiles for the four intensively
# monitored plots in Panama (C8, C13, R7, and R11).
#
# Each panel contains:
#   - An X-Z cross-section of mean voxel density along Y
#   - A smoothed canopy-height profile
#   - A grey vertical vegetation-density distribution on the right
#
# All four panels share:
#   - The same structural density colour scale
#   - The same vertical density profile scale
#   - One common legend
#
# Author: Emma Van de Walle - Q-ForestLab
###############################################################################

# ------------------------------------------------------------
# Load packages
# ------------------------------------------------------------

library(patchwork)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cowplot)
library(zoo)

# ------------------------------------------------------------
# Load TLS-derived structural grids
# ------------------------------------------------------------

tls_C8  <- readRDS("Data/TLS_scaled_DTM_and_grid_C8_Panama.rds")
tls_C13 <- readRDS("Data/TLS_scaled_DTM_and_grid_C13_Panama.rds")
tls_R7  <- readRDS("Data/TLS_scaled_DTM_and_grid_R7_Panama.rds")
tls_R11 <- readRDS("Data/TLS_scaled_DTM_and_grid_R11_Panama.rds")

grids <- list(
  C8  = tls_C8$grid,
  C13 = tls_C13$grid,
  R7  = tls_R7$grid,
  R11 = tls_R11$grid
)

# ------------------------------------------------------------
# Settings
# ------------------------------------------------------------

# Original structural profile dimensions
x_min <- 0
x_max <- 60

z_min <- 0
z_max <- 40

# Position and width of grey density distributions
# These are expressed in the same coordinate system as the X-Z plots.

# Grey density profile directly adjacent to structural profile
profile_x_start <- 60.5
profile_width <- 8

panel_x_max <- profile_x_start + profile_width

# ------------------------------------------------------------
# Calculate vertical density distributions for all plots
# ------------------------------------------------------------

vertical_profiles <- lapply(grids, function(grid) {

  grid |>
    group_by(Z) |>
    summarise(
      mean_density = mean(density, na.rm = TRUE),
      .groups = "drop"
    ) |>
    complete(
      Z = seq(z_min, z_max, by = 1),
      fill = list(mean_density = 0)
    ) |>
    filter(
      Z >= z_min,
      Z <= z_max
    )
})

# Common maximum density across all four plots.
# This ensures that profile widths are directly comparable.

max_vertical_density <- max(
  vapply(
    vertical_profiles,
    function(d) {
      max(d$mean_density, na.rm = TRUE)
    },
    numeric(1)
  ),
  na.rm = TRUE
)

if (!is.finite(max_vertical_density) ||
    max_vertical_density <= 0) {
  stop("No positive vertical vegetation density found.")
}

cat(
  "Maximum vertical density across plots:",
  max_vertical_density,
  "\n"
)

# ------------------------------------------------------------
# Function to create a structural profile
# ------------------------------------------------------------

create_structure_plot <- function(grid,
                                  site_name,
                                  show_x = FALSE,
                                  show_y = FALSE) {

  # ----------------------------------------------------------
  # Average voxel density along Y
  # ----------------------------------------------------------

  average_density <- grid |>
    group_by(X, Z) |>
    summarise(
      mean_density = mean(density, na.rm = TRUE),
      .groups = "drop"
    ) |>
    complete(
      X = seq(x_min, x_max, by = 1),
      Z = seq(z_min, z_max, by = 1),
      fill = list(mean_density = 0)
    )

  average_density$plot_density <- ifelse(
    average_density$mean_density == 0,
    NA_real_,
    average_density$mean_density
  )

  # ----------------------------------------------------------
  # Determine canopy-height profile
  # ----------------------------------------------------------

  canopy_line <- average_density |>
    filter(!is.na(plot_density)) |>
    group_by(X) |>
    summarise(
      canopy_height = max(Z),
      .groups = "drop"
    ) |>
    arrange(X)

  # Smooth canopy-height profile
  if (nrow(canopy_line) >= 6) {

    canopy_line <- canopy_line |>
      mutate(
        canopy_height = zoo::rollmean(
          canopy_height,
          k = 6,
          fill = "extend"
        )
      )
  }

  # ----------------------------------------------------------
  # Vertical vegetation-density distribution
  # ----------------------------------------------------------

  vertical_density <- vertical_profiles[[site_name]]

  vertical_density <- vertical_density |>
    mutate(
      profile_x = profile_x_start +
        (mean_density / max_vertical_density) *
        profile_width
    )

  # ----------------------------------------------------------
  # Create plot
  # ----------------------------------------------------------

  p <- ggplot() +

    # Original X-Z structural profile
    geom_tile(
      data = average_density,
      aes(
        x = X,
        y = Z,
        fill = plot_density
      ),
      width = 1,
      height = 1
    ) +

    # Smoothed canopy-height line
    geom_line(
      data = canopy_line,
      aes(
        x = X,
        y = canopy_height
      ),
      colour = "black",
      linewidth = 0.8,
      lineend = "round"
    ) +

    # Grey vertical vegetation-density distribution
    geom_ribbon(
      data = vertical_density,
      aes(
        y = Z,
        xmin = profile_x_start,
        xmax = profile_x
      ),
      orientation = "y",
      fill = "grey65",
      colour = NA
    ) +

    # Subtle separation between structure and density distribution
    geom_segment(
      x = profile_x_start,
      xend = profile_x_start,
      y = z_min,
      yend = z_max,
      colour = "grey50",
      linewidth = 0.25
    ) +

    # Original shared density colour scale
    scale_fill_viridis_c(
      option = "viridis",
      direction = -1,
      limits = c(0, 1),
      na.value = "#CFE8FF",
      name = "Normalised\ndensity",
      guide = guide_colorbar(
        barwidth = 1,
        barheight = 15,
        frame.colour = "black",
        ticks.colour = "black"
      )
    ) +

    # Keep original X and Y axis ticks
    scale_x_continuous(
      limits = c(x_min, panel_x_max),
      breaks = seq(0, 60, 20),
      expand = c(0, 0)
    ) +

    scale_y_continuous(
      limits = c(z_min, z_max),
      breaks = seq(0, 40, 10),
      expand = c(0, 0)
    ) +

    # Fixed aspect ratio for structural profile
    coord_fixed(
      ratio = 1,
      expand = FALSE
    ) +

    labs(
      title = site_name,
      x = NULL,
      y = NULL
    ) +

    theme_bw(base_size = 18) +
    theme(
      plot.title = element_text(
        size = 22,
        face = "bold",
        hjust = 0.5
      ),

      axis.title = element_text(
        size = 22
      ),

      axis.text = element_text(
        size = 20
      ),

      strip.text = element_text(
        size = 20,
        face = "bold"
      ),

      legend.title = element_text(
        size = 20,
        margin = margin(b = 25)
      ),

      legend.text = element_text(
        size = 22
      )
    )

  # ----------------------------------------------------------
  # Show X-axis labels only for bottom panels
  # ----------------------------------------------------------

  if (!show_x) {

    p <- p +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }

  # ----------------------------------------------------------
  # Show Y-axis labels only for left panels
  # ----------------------------------------------------------

  if (!show_y) {

    p <- p +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  }

  return(p)
}

# ------------------------------------------------------------
# Create four structural profiles
# ------------------------------------------------------------

p1 <- create_structure_plot(
  grids$C8,
  "C8",
  show_y = TRUE
)

p2 <- create_structure_plot(
  grids$C13,
  "C13"
)

p3 <- create_structure_plot(
  grids$R7,
  "R7",
  show_x = TRUE,
  show_y = TRUE
)

p4 <- create_structure_plot(
  grids$R11,
  "R11",
  show_x = TRUE
)

# ------------------------------------------------------------
# Combine the four structural profiles
# ------------------------------------------------------------

# Original 2 x 2 arrangement:
#
# C8   C13
# R7   R11
#
# The grey vertical profiles are now part of each plot,
# so no additional patchwork panels are needed.

final_plot <- (
  (p1 + p2) /
    (p3 + p4)
) +
  plot_layout(
    guides = "collect"
  ) +
  plot_annotation(
    theme = theme(
      plot.margin = margin(
        t = 20,
        r = 20,
        b = 30,
        l = 30
      )
    )
  ) &
  theme(
    legend.position = "right"
  )

# ------------------------------------------------------------
# Add common X and Y axis titles
# ------------------------------------------------------------

final_plot_axes <- ggdraw() +

  draw_plot(
    final_plot,
    x = 0.025,
    y = 0.045,
    width = 0.965,
    height = 0.94
  ) +

  draw_label(
    "Height (m)",
    x = 0.035,
    y = 0.52,
    angle = 90,
    size = 22
  ) +

  draw_label(
    "Easting (m)",
    x = 0.45,
    y = 0.09,
    size = 22
  )

# ------------------------------------------------------------
# Display figure
# ------------------------------------------------------------

print(final_plot_axes)

# ------------------------------------------------------------
# Save figure
# ------------------------------------------------------------

dir.create(
  "Output",
  showWarnings = FALSE,
  recursive = TRUE
)

ggsave(
  "Output/Figure_forest_structure.png",
  plot = final_plot_axes,
  width = 14,
  height = 10,
  dpi = 600
)


