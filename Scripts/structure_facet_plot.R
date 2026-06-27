library(patchwork)
library(dplyr)
library(ggplot2)

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
      average_density$mean_density^2
    )

  p <- ggplot(
    average_density,
    aes(X, Z, fill = plot_density)
  ) +
    geom_tile() +
    # scale_fill_gradient(
    #   low = "lightgreen",
    #   high = "darkgreen",
    #   limits = c(0, 1),
    #   trans = "sqrt",
    #   na.value = "#CFE8FF",
    #   name = "Normalised\ndensity",
    #   guide = guide_colorbar(
    #     barwidth = 1,
    #     barheight = 12,
    #     frame.colour = "black",
    #     ticks.colour = "black"
    #   )
    # ) +
    scale_fill_viridis_c(
      option = "viridis",
      direction = -1,
      trans = "sqrt",
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

tls_C8  <- readRDS("Data/TLS_scaled_DTM_and_grid_C8_Panama.rds")
tls_C13 <- readRDS("Data/TLS_scaled_DTM_and_grid_C13_Panama.rds")
tls_R7  <- readRDS("Data/TLS_scaled_DTM_and_grid_R7_Panama.rds")
tls_R11 <- readRDS("Data/TLS_scaled_DTM_and_grid_R11_Panama.rds")

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


final_plot <-
  (p1 + p2) /
  (p3 + p4) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title =
      "Forest structural profiles of liana-present control plots (C) \nand liana-removal plots (R) in Panama.",
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

library(cowplot)

final_plot_axes <- ggdraw() +
  draw_plot(final_plot) +
  draw_label(
    "Height (m)",
    x = 0.03,
    y = 0.5,
    angle = 90,
    size = 22
  ) +
  draw_label(
    "Easting (m)",
    x = 0.5,
    y = 0.04,
    size = 22
  )

final_plot_axes

ggsave(
  "Output/Figure_forest_structure.png",
  final_plot_axes,
  width = 14,
  height = 10,
  dpi = 600
)


