############################################################################################
# ForEdgeClim - 3D thermal heterogeneity
#
# Main figure:
#   (a) Vertical distribution of spatial thermal range
#   (b) Treatment difference across forest height
#
# Appendix figure:
#   Distribution of spatial thermal range across vertical forest zones
#
# Author: Emma Van de Walle - Q-ForestLab
############################################################################################


# ---- Packages ----
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(gridExtra)
  library(grid)
})


# ==========================================================================================
# INPUT
# ==========================================================================================

control <- readRDS(
  "Output/heterogeneity/heterogeneity_control_calibrated.rds"
)

removal <- readRDS(
  "Output/heterogeneity/heterogeneity_removal_calibrated.rds"
)


# ---- Add treatment ----

control <- control %>%
  mutate(
    treatment = "Control"
  )

removal <- removal %>%
  mutate(
    treatment = "Liana removal"
  )


# ==========================================================================================
# COMBINE DATA
# ==========================================================================================

dat <- bind_rows(
  control,
  removal
) %>%
  mutate(

    treatment = factor(
      treatment,
      levels = c(
        "Control",
        "Liana removal"
      )
    ),

    # Seasons
    season = case_when(
      doy >= 350 | doy <= 120 ~ "Dry season",
      TRUE                    ~ "Wet season"
    ),

    # Day/night
    period = case_when(
      hour >= 11 & hour <= 22 ~ "Day",
      TRUE                    ~ "Night"
    ),

    season = factor(
      season,
      levels = c(
        "Dry season",
        "Wet season"
      )
    ),

    period = factor(
      period,
      levels = c(
        "Day",
        "Night"
      )
    )
  )


# ==========================================================================================
# COLOURS
# ==========================================================================================

treatment_colours <- c(
  "Control"       = "blue",
  "Liana removal" = "red"
)


# ==========================================================================================
# OUTPUT
# ==========================================================================================

output_path <- "Output/heterogeneity/figures/"

dir.create(
  output_path,
  recursive = TRUE,
  showWarnings = FALSE
)



############################################################################################
############################################################################################
#
# MAIN FIGURE
#
############################################################################################
############################################################################################


# ==========================================================================================
# STEP 1
# SUMMARIZE THROUGH TIME WITHIN EACH INDIVIDUAL PLOT
#
# Result:
#
# treatment x plot x season x period x height
#
# Spatial thermal range = P95(Tair) - P05(Tair)
#
# Every plot therefore receives equal weight.
# ==========================================================================================

plot_summary <- dat %>%
  group_by(
    treatment,
    plot,
    season,
    period,
    z
  ) %>%
  summarise(

    Tair_range90 = median(
      Tair_range90,
      na.rm = TRUE
    ),

    .groups = "drop"
  )


# ---- Remove non-finite values ----

plot_summary <- plot_summary[
  is.finite(plot_summary$Tair_range90),
]


# ==========================================================================================
# STEP 2
# SUMMARIZE ACROSS PLOTS WITHIN EACH TREATMENT
#
# Median = treatment-level central tendency
# Q25-Q75 = between-plot variability
# ==========================================================================================

vertical_summary <- plot_summary %>%
  group_by(
    treatment,
    season,
    period,
    z
  ) %>%
  summarise(

    median_range90 = median(
      Tair_range90,
      na.rm = TRUE
    ),

    q25_range90 = quantile(
      Tair_range90,
      probs = 0.25,
      na.rm = TRUE,
      names = FALSE
    ),

    q75_range90 = quantile(
      Tair_range90,
      probs = 0.75,
      na.rm = TRUE,
      names = FALSE
    ),

    .groups = "drop"
  )


# ---- Remove possible non-finite values ----

vertical_summary <- vertical_summary[
  is.finite(vertical_summary$median_range90) &
    is.finite(vertical_summary$q25_range90) &
    is.finite(vertical_summary$q75_range90),
]


# ==========================================================================================
# FIGURE A
#
# VERTICAL DISTRIBUTION OF SPATIAL THERMAL HETEROGENEITY
#
# Separate figures are created for daytime and nighttime.
#
# Broad coloured marks:
#   treatment median
#
# Transparent coloured intervals:
#   Q25-Q75 among plots
# ==========================================================================================


# ---- Function for Figure A ----

create_plot_A <- function(data, selected_period) {

  plot_dat <- data %>%
    filter(
      period == selected_period
    )


  ggplot(
    plot_dat,
    aes(
      y = z,
      colour = treatment
    )
  ) +

    # ---- Interquartile range among plots ----

  geom_ribbon(
    aes(
      xmin = q25_range90,
      xmax = q75_range90,
      group = treatment,
      fill = treatment
    ),
    alpha = 0.20,
    colour = NA,
    orientation = "y"
  ) +

    # ---- Treatment median: connected vertically across height ----

  geom_path(
    aes(
      x = median_range90,
      group = treatment
    ),
    linewidth = 1.5
  ) +

    # ---- Facets: seasons only ----

  facet_grid(
    . ~ season
  ) +

  scale_y_continuous(
    limits = c(0, 40),
    breaks = seq(0, 40, by = 10),
    expand = expansion(mult = c(0, 0))
  ) +

    # ---- Colours ----

  scale_colour_manual(
    values = treatment_colours,
    name = NULL
  ) +

    scale_fill_manual(
      values = treatment_colours,
      name = NULL
    ) +

    guides(
      colour = guide_legend(
        override.aes = list(
          linewidth = 1.5,
          alpha = 1
        )
      ),
      fill = "none"
    ) +

    # ---- Labels ----

  labs(

    x = expression(
      "Spatial thermal range (" *
        P[95] - P[5] *
        ", " *
        degree * C *
        ")"
    ),

    y = "Height (m)",

    title = "(a) Vertical distribution of spatial thermal heterogeneity"
  ) +

    # ---- Theme ----

  theme_bw() +

    theme(

      plot.title = element_text(
        size = 34
      ),

      axis.title = element_text(
        size = 32
      ),

      axis.text = element_text(
        size = 28
      ),

      strip.text = element_text(
        size = 35,
        face = "bold"
      ),

      legend.text = element_text(
        size = 28
      ),

      legend.position = "top",

      legend.justification = "left",

      legend.key.width = unit(
        2.5,
        "cm"
      ),

      legend.margin = margin(
        t = 0,
        r = 0,
        b = 5,
        l = 0
      ),

      plot.margin = margin(
        10,
        15,
        10,
        10
      )
    )
}


# ==========================================================================================
# TREATMENT DIFFERENCE
#
# Liana removal - Control
#
# Difference calculated for all pairwise combinations of
# removal and control plots at each height.
#
# Line   = median pairwise treatment difference
# Ribbon = Q25-Q75 among pairwise treatment differences
# ==========================================================================================

control_diff <- plot_summary %>%
  filter(treatment == "Control") %>%
  select(
    plot_control = plot,
    season,
    period,
    z,
    range_control = Tair_range90
  )

removal_diff <- plot_summary %>%
  filter(treatment == "Liana removal") %>%
  select(
    plot_removal = plot,
    season,
    period,
    z,
    range_removal = Tair_range90
  )


# ---- All 4 x 4 treatment combinations at each height ----

pairwise_differences <- removal_diff %>%
  inner_join(
    control_diff,
    by = c("season", "period", "z"),
    relationship = "many-to-many"
  ) %>%
  mutate(
    delta_range90 = range_removal - range_control
  )


treatment_difference <- pairwise_differences %>%
  group_by(
    season,
    period,
    z
  ) %>%
  summarise(

    median_delta = median(
      delta_range90,
      na.rm = TRUE
    ),

    q25_delta = quantile(
      delta_range90,
      probs = 0.25,
      na.rm = TRUE,
      names = FALSE
    ),

    q75_delta = quantile(
      delta_range90,
      probs = 0.75,
      na.rm = TRUE,
      names = FALSE
    ),

    .groups = "drop"
  ) %>%

  # Rename only AFTER all summaries have been calculated
  rename(
    delta_range90 = median_delta
  )

# ==========================================================================================
# FIGURE B
#
# TREATMENT DIFFERENCE ACROSS FOREST HEIGHT
# ==========================================================================================


# ---- Function for Figure B ----

create_plot_B <- function(data, selected_period) {

  plot_dat <- data %>%

    filter(
      period == selected_period
    )


  ggplot(
    plot_dat,
    aes(
      x = delta_range90,
      y = z
    )
  ) +

    # ---- Zero = no treatment difference ----

  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    colour = "grey30",
    linewidth = 1.3
  ) +

    # ---- Difference at each 1-m layer ----

  # ---- Q25-Q75 of pairwise treatment differences ----

  geom_ribbon(
    aes(
      y = z,
      xmin = q25_delta,
      xmax = q75_delta,
      group = season
    ),
    orientation = "y",
    inherit.aes = FALSE,
    fill = "grey70",
    alpha = 0.5,
    colour = NA
  ) +

    # ---- Median pairwise treatment difference ----

  geom_path(
    aes(
      x = delta_range90,
      y = z,
      group = season
    ),
    inherit.aes = FALSE,
    linewidth = 1.5,
    colour = "grey20"
  ) +

    # ---- Facets: seasons only ----

  facet_grid(
    . ~ season
  ) +

  scale_y_continuous(
    limits = c(0, 40),
    breaks = seq(0, 40, by = 10),
    expand = expansion(mult = c(0, 0))
  ) +

    # ---- Labels ----

  labs(

    x = expression(
      Delta *
        " spatial thermal range (removal - control, " *
        degree * C *
        ")"
    ),

    y = "Height (m)",

    title = "(b) Liana removal effect on spatial thermal heterogeneity"
  ) +

    # ---- Theme ----

  theme_bw() +

    theme(

      plot.title = element_text(
        size = 34
      ),

      axis.title = element_text(
        size = 32
      ),

      axis.text = element_text(
        size = 28
      ),

      strip.text = element_text(
        size = 35,
        face = "bold"
      ),

      legend.position = "none",

      plot.margin = margin(
        10,
        15,
        10,
        10
      )
    )
}


# ==========================================================================================
# CREATE DAY FIGURE
# ==========================================================================================

plot_A_day <- create_plot_A(
  vertical_summary,
  "Day"
)

plot_B_day <- create_plot_B(
  treatment_difference,
  "Day"
)


main_heterogeneity_plot_day <- gridExtra::arrangeGrob(

  plot_A_day,

  plot_B_day,

  ncol = 1,

  heights = c(
    1,
    1
  )
)


# ==========================================================================================
# CREATE NIGHT FIGURE
# ==========================================================================================

plot_A_night <- create_plot_A(
  vertical_summary,
  "Night"
)

plot_B_night <- create_plot_B(
  treatment_difference,
  "Night"
)


main_heterogeneity_plot_night <- gridExtra::arrangeGrob(

  plot_A_night,

  plot_B_night,

  ncol = 1,

  heights = c(
    1,
    1
  )
)


# ==========================================================================================
# SAVE DAY FIGURE
# ==========================================================================================

ggsave(

  filename = file.path(
    output_path,
    "vertical_thermal_heterogeneity_control_vs_removal_day.png"
  ),

  plot = main_heterogeneity_plot_day,

  width = 13,
  height = 13,

  dpi = 300
)


# ==========================================================================================
# SAVE NIGHT FIGURE
# ==========================================================================================

ggsave(

  filename = file.path(
    output_path,
    "vertical_thermal_heterogeneity_control_vs_removal_night.png"
  ),

  plot = main_heterogeneity_plot_night,

  width = 13,
  height = 13,

  dpi = 300
)


############################################################################################
############################################################################################
#
# APPENDIX FIGURE
#
# DISTRIBUTION OF SPATIAL THERMAL HETEROGENEITY
# ACROSS VERTICAL FOREST ZONES
#
############################################################################################
############################################################################################


# ==========================================================================================
# DEFINE VERTICAL FOREST ZONES
#
# 0-2 m
# 2-10 m
# 10-20 m
# 20-30 m
# >30 m
# ==========================================================================================

vertical_zone_data <- dat %>%

  mutate(

    height_class = cut(

      z,

      breaks = c(
        -Inf,
        2,
        10,
        20,
        30,
        Inf
      ),

      labels = c(
        "0-2 m",
        "2-10 m",
        "10-20 m",
        "20-30 m",
        ">30 m"
      ),

      right = FALSE
    )
  )

# ---- Reverse vertical order for plotting ----

vertical_zone_data <- vertical_zone_data %>%
  mutate(
    height_class = factor(
      height_class,
      levels = c(
        ">30 m",
        "20-30 m",
        "10-20 m",
        "2-10 m",
        "0-2 m"
      )
    )
  )

# ---- Remove missing height classes ----

vertical_zone_data <- vertical_zone_data[
  !is.na(vertical_zone_data$height_class),
]


# ==========================================================================================
# STEP 1
#
# Average hourly heterogeneity within each:
#
# plot x day x period x 1-m layer
# ==========================================================================================

vertical_zone_data <- vertical_zone_data %>%

  group_by(
    treatment,
    plot,
    season,
    period,
    date,
    height_class,
    z
  ) %>%

  summarise(

    Tair_range90 = mean(
      Tair_range90,
      na.rm = TRUE
    ),

    .groups = "drop"
  )


# ==========================================================================================
# STEP 2
#
# Average 1-m layers within each vertical forest zone.
#
# This prevents wider zones from receiving more weight simply because
# they contain more individual 1-m layers.
# ==========================================================================================

vertical_zone_data <- vertical_zone_data %>%

  group_by(
    treatment,
    plot,
    season,
    period,
    date,
    height_class
  ) %>%

  summarise(

    Tair_range90 = mean(
      Tair_range90,
      na.rm = TRUE
    ),

    .groups = "drop"
  )


# ---- Remove non-finite observations ----

vertical_zone_data <- vertical_zone_data[
  is.finite(vertical_zone_data$Tair_range90),
]


# ==========================================================================================
# CREATE TREATMENT POSITIONS
#
# Control       = position 1
# Liana removal = position 2
#
# Treatment names are deliberately not shown on the x-axis.
# The colours and legend identify the treatments.
# ==========================================================================================

vertical_zone_data <- vertical_zone_data %>%

  mutate(

    treatment_position = case_when(

      treatment == "Control" ~ 1,

      treatment == "Liana removal" ~ 2
    )
  )


# ==========================================================================================
# APPENDIX VIOLIN PLOT
# ==========================================================================================

vertical_zone_plot <- ggplot(
  vertical_zone_data,
  aes(
    x = treatment_position,
    y = Tair_range90,
    fill = treatment,
    colour = treatment,
    group = treatment
  )
) +

  # ---- Violin distributions ----

geom_violin(
  trim = TRUE,
  alpha = 0.25,
  scale = "width",
  linewidth = 1.3,
  width = 0.8
) +

  # ---- Boxplots ----

geom_boxplot(
  width = 0.16,
  outlier.shape = NA,
  alpha = 0.65,
  linewidth = 1
) +

  # ---- Facets ----

facet_grid(

  height_class ~

    interaction(
      period,
      season,
      sep = " - "
    ),

  scales = "free_y"
) +

  # ---- Colours ----

scale_fill_manual(
  values = treatment_colours,
  name = NULL
) +

  scale_colour_manual(
    values = treatment_colours,
    name = NULL
  ) +

  # ---- X-axis ----
#
# Treatments remain separated spatially,
# but treatment names are not printed underneath.
# ----

scale_x_continuous(
  limits = c(
    0.4,
    2.6
  ),
  breaks = NULL
) +

  # ---- Labels ----

labs(

  x = NULL,

  y = expression(
    "Spatial thermal range (" *
      P[95] - P[5] *
      ", " *
      degree * C *
      ")"
  )
) +

  # ---- Theme ----

theme_bw() +

  theme(

    axis.title = element_text(
      size = 32
    ),

    axis.text.y = element_text(
      size = 28
    ),

    # No treatment labels underneath violins
    axis.text.x = element_blank(),

    axis.ticks.x = element_blank(),

    strip.text = element_text(
      size = 30,
      face = "bold"
    ),

    legend.text = element_text(
      size = 28,
      margin = margin(l = 12, r = 25)
    ),

    legend.position = "top",

    legend.key.width = unit(
      2.5,
      "cm"
    ),

    plot.margin = margin(
      10,
      15,
      10,
      10
    )
  )

# ==========================================================================================
# SAVE APPENDIX FIGURE
# ==========================================================================================

ggsave(

  filename = file.path(
    output_path,
    "vertical_zone_thermal_heterogeneity_distributions.png"
  ),

  plot = vertical_zone_plot,

  width = 18,
  height = 16,

  dpi = 300
)
