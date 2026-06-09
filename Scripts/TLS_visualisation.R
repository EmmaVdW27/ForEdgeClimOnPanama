# =========================================================
# BEAUTIFUL 3D VOXEL VISUALISATION OF TROPICAL TLS PLOTS
# Panama liana removal experiment
# =========================================================
#
# This script:
# 1. Reads four LAS files
# 2. Creates voxel density grids
# 3. Creates beautiful interactive 3D renders
# 4. Uses identical colour scales for all plots
# 5. Uses identical camera angle for all plots
# 6. Saves HTML files for manual PNG export
#
# Recommended workflow:
# - Open HTML files
# - Rotate slightly if desired
# - Use screenshot tool
# - Combine PNGs afterwards
#
# =========================================================

# -----------------------------
# LOAD LIBRARIES
# -----------------------------

library(lidR)
library(data.table)
library(plotly)
library(viridis)
library(htmlwidgets)

# -----------------------------
# SETTINGS
# -----------------------------

voxel_size <- 1

density_threshold <- 0.08

# -----------------------------
# CREATE OUTPUT FOLDER
# -----------------------------

dir.create(
  "Output/TLS_Figures",
  recursive = TRUE,
  showWarnings = FALSE
)

# =========================================================
# FUNCTION TO CREATE VOXEL GRID
# =========================================================

generate_voxel_grid <- function(
    las_file,
    voxel_size = 1
) {

  cat("\n")
  cat("Processing:\n")
  cat(las_file, "\n")

  # -----------------------------
  # READ LAS
  # -----------------------------

  las <- readLAS(las_file)

  # -----------------------------
  # CLASSIFY GROUND
  # -----------------------------

  las_ground <- classify_ground(
    las,
    csf()
  )

  # -----------------------------
  # CREATE DTM
  # -----------------------------

  dtm <- rasterize_terrain(
    las_ground,
    res = 1,
    algorithm = knnidw(k = 10L)
  )

  # -----------------------------
  # NORMALISE HEIGHT
  # -----------------------------

  las <- normalize_height(
    las,
    dtm
  )

  # Remove points below ground
  las <- filter_poi(
    las,
    Z >= 0
  )

  # -----------------------------
  # CREATE VOXEL GRID
  # -----------------------------

  voxels <- voxel_metrics(
    las,
    ~length(Z),
    res = voxel_size
  )

  voxel_grid <- as.data.table(voxels)

  setnames(
    voxel_grid,
    c("X", "Y", "Z", "n")
  )

  # -----------------------------
  # REMOVE OUTER BOUNDARY
  # -----------------------------

  voxel_grid <- voxel_grid[
    X > min(X) &
      X < max(X)
  ]

  voxel_grid <- voxel_grid[
    Y > min(Y) &
      Y < max(Y)
  ]

  # -----------------------------
  # RESCALE COORDINATES
  # -----------------------------

  voxel_grid$X <- voxel_grid$X -
    min(voxel_grid$X) + 1

  voxel_grid$Y <- voxel_grid$Y -
    min(voxel_grid$Y) + 1

  voxel_grid$Z <- voxel_grid$Z -
    min(voxel_grid$Z) + 1

  # -----------------------------
  # FILL MISSING VOXELS
  # -----------------------------

  full_grid <- CJ(
    X = seq(
      min(voxel_grid$X),
      max(voxel_grid$X),
      by = voxel_size
    ),

    Y = seq(
      min(voxel_grid$Y),
      max(voxel_grid$Y),
      by = voxel_size
    ),

    Z = seq(
      min(voxel_grid$Z),
      max(voxel_grid$Z),
      by = voxel_size
    )
  )

  full_grid <- merge(
    full_grid,
    voxel_grid,
    by = c("X", "Y", "Z"),
    all.x = TRUE
  )

  full_grid[
    is.na(n),
    n := 0
  ]

  # -----------------------------
  # DENSITY TRANSFORMATION
  # -----------------------------

  scaling <- log1p(full_grid$n) *
    exp(
      0.5 * full_grid$Z /
        max(full_grid$Z)
    )

  full_grid$density <- scaling /
    max(scaling)

  # -----------------------------
  # REMOVE LOW DENSITY VOXELS
  # -----------------------------

  full_grid <- full_grid[
    density > density_threshold
  ]

  return(full_grid)
}

# =========================================================
# READ LAS FILES
# =========================================================

plot07 <- generate_voxel_grid(
  "Data/Data_Kasper/plot07_trimmed_1cm.las",
  voxel_size
)

plot08 <- generate_voxel_grid(
  "Data/Data_Kasper/plot08_trimmed_1cm.las",
  voxel_size
)

plot11 <- generate_voxel_grid(
  "Data/Data_Kasper/plot11_trimmed_5cm.las",
  voxel_size
)

plot13 <- generate_voxel_grid(
  "Data/Data_Kasper/plot13_trimmed_5cm.las",
  voxel_size
)

# =========================================================
# GLOBAL COLOUR SCALE
# SAME FOR ALL PLOTS
# =========================================================

global_zmin <- min(
  plot07$Z,
  plot08$Z,
  plot11$Z,
  plot13$Z
)

global_zmax <- max(
  plot07$Z,
  plot08$Z,
  plot11$Z,
  plot13$Z
)

# =========================================================
# FIXED CAMERA
# =========================================================

camera <- list(

  eye = list(
    x = 1.4,
    y = 2.1,
    z = 0.8
  ),

  up = list(
    x = 0,
    y = 0,
    z = 1
  )
)

# =========================================================
# FUNCTION TO CREATE BEAUTIFUL RENDER
# =========================================================

make_voxel_render <- function(
    grid,
    title,
    output_html
) {

  fig <- plot_ly()

  fig <- fig %>%

    add_markers(

      data = grid,

      x = ~X,
      y = ~Y,
      z = ~Z,

      # Colour by height
      color = ~Z,

      # Same scale for all plots
      cmin = global_zmin,
      cmax = global_zmax,

      colors = viridis(100),

      # Transparency
      opacity = ~pmax(density, 0.35),

      # Remove legend
      showscale = FALSE,

      marker = list(

        # Bigger voxels
        size = 5
      )
    )

  # -----------------------------
  # LAYOUT
  # -----------------------------

  fig <- fig %>%

    layout(

      title = list(

        width = 1200,
        height = 1200,

        text = title,

        x = 0.5,
        y = 0.96,

        xanchor = "center",

        font = list(
          size = 26,
          color = "white"
        )
      ),

      paper_bgcolor = "black",
      plot_bgcolor = "black",

      margin = list(
        l = 0,
        r = 0,
        b = 0,
        t = 70
      ),

      scene = list(

        bgcolor = "black",

        camera = camera,

        aspectmode = "data",

        xaxis = list(
          visible = FALSE
        ),

        yaxis = list(
          visible = FALSE
        ),

        zaxis = list(
          visible = FALSE
        )
      )
    )

  # -----------------------------
  # SAVE HTML
  # -----------------------------

  saveWidget(
    fig,
    output_html,
    selfcontained = TRUE
  )

  return(fig)
}

# =========================================================
# CREATE FIGURES
# =========================================================

fig07 <- make_voxel_render(
  plot07,
  "Plot 07 - Removal",
  "Output/TLS_Figures/plot07_voxel.html"
)

fig08 <- make_voxel_render(
  plot08,
  "Plot 08 - Control",
  "Output/TLS_Figures/plot08_voxel.html"
)

fig11 <- make_voxel_render(
  plot11,
  "Plot 11 - Removal",
  "Output/TLS_Figures/plot11_voxel.html"
)

fig13 <- make_voxel_render(
  plot13,
  "Plot 13 - Control",
  "Output/TLS_Figures/plot13_voxel.html"
)

# =========================================================
# BEAUTIFUL SEPARATE LEGEND
# White background + black ticks + border ONLY around scale
# =========================================================

library(ggplot2)
library(viridis)

# -----------------------------
# Dummy data
# -----------------------------

legend_data <- data.frame(
  x = 1,
  y = seq(
    global_zmin,
    global_zmax,
    length.out = 100
  )
)

# -----------------------------
# Create legend plot
# -----------------------------

legend_plot <- ggplot(
  legend_data,
  aes(
    x = x,
    y = y,
    fill = y
  )
) +

  geom_tile() +

  scale_fill_viridis_c(

    option = "viridis",

    name = "Height (m)",

    limits = c(
      global_zmin,
      global_zmax
    )
  ) +

  theme_void() +

  theme(

    # White background
    plot.background = element_rect(
      fill = "black",
      colour = "black"
    ),

    panel.background = element_rect(
      fill = "black",
      colour = "black"
    ),

    # Legend title
    legend.title = element_text(
      colour = "white",
      size = 18
    ),

    # Legend tick labels
    legend.text = element_text(
      colour = "white",
      size = 14
    ),

    # Legend position
    legend.position = "right"

  ) +

  guides(

    fill = guide_colourbar(

      # Black ticks
      ticks.colour = "white",

      ticks.linewidth = 1,

      # Black border ONLY around scale
      frame.colour = "white",

      frame.linewidth = 1,

      # Scale dimensions
      barwidth = 1.2,
      barheight = 20
    )
  )

# -----------------------------
# Save legend
# -----------------------------

ggsave(

  filename = "Output/TLS_Figures/legend.png",

  plot = legend_plot,

  width = 2.5,
  height = 5,

  dpi = 400,

  bg = "black"
)
# =========================================================
# DISPLAY FIGURES
# =========================================================

fig07
fig08
fig11
fig13

# =========================================================
# FINISHED
# =========================================================

cat("\n")
cat("Finished!\n")
cat("Figures saved in:\n")
cat("Output/TLS_Figures/\n")
