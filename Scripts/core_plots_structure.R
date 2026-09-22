###############################################################################
# Forest structural profile figure
#
# Creates vertical forest structural profiles for the four intensively
# monitored plots in Panama (C8, C13, R7, and R11) from 1-m TLS-derived
# structural grids. Voxel density is averaged along the northing (Y) direction
# to obtain an X-Z cross-sectional representation of forest structure.
# A smoothed canopy-height profile is overlaid on each structural profile.
# The four plots are combined into a single figure with a shared density scale.
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

# ------------------------------------------------------------
# Function to create a vertical forest structural profile
# ------------------------------------------------------------

create_structure_plot <- function(grid,
                                  site_name,
                                  show_x = FALSE,
                                  show_y = FALSE){

  average_density <- grid |>
    dplyr::group_by(X, Z) |>
    dplyr::summarise(
      mean_density = mean(density, na.rm = TRUE),
      .groups = "drop"
    ) |>
    tidyr::complete(
      X = seq(min(X), max(X), by = 1),
      Z = seq(0, 40, by = 1),
      fill = list(mean_density = 0)
    )

  average_density$plot_density <-
    ifelse(
      average_density$mean_density == 0,
      NA,
      average_density$mean_density
    )

  # Determine canopy surface (highest voxel containing vegetation)
  canopy_line <- average_density %>%
    filter(!is.na(plot_density)) %>%
    group_by(X) %>%
    summarise(
      canopy_height = max(Z),
      .groups = "drop"
    ) %>%
    mutate(
      canopy_height = zoo::rollmean(canopy_height,
                                    k = 6,
                                    fill = "extend")
    )

  p <- ggplot(
    average_density,
    aes(X, Z, fill = plot_density)
  ) +
    geom_tile() +
    geom_line(
      data = canopy_line,
      aes(x = X, y = canopy_height),
      inherit.aes = FALSE,
      colour = "black",
      linewidth = 0.8,
      lineend = "round"
    ) +
    scale_fill_viridis_c(
      option = "viridis",
      direction = -1,
      #trans = "sqrt",
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
    coord_fixed() +
    scale_x_continuous(
      limits = c(0, 60),
      breaks = seq(0, 60, 20)
    ) +
    scale_y_continuous(
      limits = c(0, 40),
      breaks = seq(0, 40, 10)
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

  if(!show_x){
    p <- p +
      theme(
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank()
      )
  }

  if(!show_y){
    p <- p +
      theme(
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
      )
  }

  p
}

# ------------------------------------------------------------
# Load TLS-derived structural grids
# ------------------------------------------------------------

tls_C8  <- readRDS("Data/TLS_scaled_DTM_and_grid_C8_Panama.rds")
tls_C13 <- readRDS("Data/TLS_scaled_DTM_and_grid_C13_Panama.rds")
tls_R7  <- readRDS("Data/TLS_scaled_DTM_and_grid_R7_Panama.rds")
tls_R11 <- readRDS("Data/TLS_scaled_DTM_and_grid_R11_Panama.rds")

# ------------------------------------------------------------
# Create vertical forest structural profiles
# ------------------------------------------------------------

p1 <- create_structure_plot(
  tls_C8$grid,
  "C8",
  show_y = TRUE
)

p2 <- create_structure_plot(
  tls_C13$grid,
  "C13"
)

p3 <- create_structure_plot(
  tls_R7$grid,
  "R7",
  show_x = TRUE,
  show_y = TRUE
)

p4 <- create_structure_plot(
  tls_R11$grid,
  "R11",
  show_x = TRUE
)

# ------------------------------------------------------------
# Combine the four structural profiles
# ------------------------------------------------------------

final_plot <-
  (p1 + p2) /
  (p3 + p4) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = NULL, #"Forest structural profiles of liana-present control plots (C) \nand liana-removal plots (R) in Panama."
    theme = theme(
      plot.title = element_text(
        size = 24,
        face = "bold",
        hjust = 0.5
      ),
      plot.margin = margin(
        t = 20,
        r = 20,
        b = 30,
        l = 30
      )
    )
  )

# ------------------------------------------------------------
# Add common x- and y-axis titles
# ------------------------------------------------------------

final_plot_axes <- ggdraw() +
  draw_plot(final_plot) +
  draw_label(
    "Height (m)",
    x = 0.012,
    y = 0.5,
    angle = 90,
    size = 22
  ) +
  draw_label(
    "Easting (m)",
    x = 0.46,
    y = 0.04,
    size = 22
  )

# ------------------------------------------------------------
# Display figure
# ------------------------------------------------------------

final_plot_axes

# ------------------------------------------------------------
# Save figure
# ------------------------------------------------------------

ggsave(
  "Output/Figure_forest_structure.png",
  final_plot_axes,
  width = 14,
  height = 10,
  dpi = 600
)


