###############################################################################
# ForEdgeClim - Visualisation of Sobol sensitivity analysis
#
# Processes and visualises total-order Sobol sensitivity indices for control
# and liana-removal treatments under dry- and wet-season daytime and nighttime
# scenarios.
#
# Total-order indices are truncated at zero and normalised within each
# treatment x season x period x metric combination. Parameters are ranked by
# their mean normalised total-order sensitivity, and the smallest common set
# explaining at least 50% of the normalised sensitivity in every corresponding
# scenario is retained separately for temperature and shortwave metrics.
# Remaining parameters are grouped as "Other parameters".
#
# The script creates a stacked bar plot of normalised total-order Sobol indices
# and calculates treatment differences and sensitivity coverage of the retained
# parameters.
#
# Author: Emma Van de Walle - Q-ForestLab
###############################################################################

library(dplyr)
library(tidyr)
library(purrr)
library(stringr)
library(ggplot2)
library(forcats)
library(viridis)

###############################################################################
# INPUT
###############################################################################

input_dir <- "Output/sensitivity_analysis/data/all_scenarios"

output_dir <- "Output/sensitivity_analysis/plots"

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

###############################################################################
# METRICS TO USE
###############################################################################

metrics_to_use <- c(
  "mean_understorey_T",
  "sd_understorey_T",
  "vertical_T_gradient",
  "mean_understorey_SW",
  "sd_understorey_SW",
  "vertical_SW_extinction"
)

###############################################################################
# SCENARIOS
###############################################################################

scenario_dirs <- list.dirs(
  input_dir,
  recursive = FALSE,
  full.names = TRUE
)

###############################################################################
# CHECK FILES
###############################################################################

missing_files <- c()

for(dir in scenario_dirs){

  scenario_name <- basename(dir)

  parts <- str_split(
    scenario_name,
    "_",
    simplify = TRUE
  )

  period <- parts[3]

  metrics_expected <- c(
    "mean_understorey_T",
    "sd_understorey_T",
    "vertical_T_gradient"
  )

  if(period == "noon"){

    metrics_expected <- c(
      metrics_expected,
      "mean_understorey_SW",
      "sd_understorey_SW",
      "vertical_SW_extinction"
    )

  }

  for(metric in metrics_expected){

    f <- file.path(
      dir,
      paste0(metric, ".rds")
    )

    if(!file.exists(f)){

      missing_files <- c(
        missing_files,
        f
      )
    }
  }
}

if(length(missing_files) > 0){

  stop(
    paste(
      "Missing files:\n",
      paste(
        missing_files,
        collapse = "\n"
      )
    )
  )
}

###############################################################################
# READ SOBOL RESULTS
###############################################################################

sobol_list <- list()

for(dir in scenario_dirs){

  scenario_name <- basename(dir)

  parts <- str_split(
    scenario_name,
    "_",
    simplify = TRUE
  )

  treatment <- parts[1]
  season    <- parts[2]
  period    <- parts[3]

  metrics_available <- c(
    "mean_understorey_T",
    "sd_understorey_T",
    "vertical_T_gradient"
  )

  if(period == "noon"){

    metrics_available <- c(
      metrics_available,
      "mean_understorey_SW",
      "sd_understorey_SW",
      "vertical_SW_extinction"
    )

  }

  for(metric in metrics_available){

    f <- file.path(
      dir,
      paste0(metric, ".rds")
    )

    sobol_res <- readRDS(f)

    tmp <- data.frame(

      treatment = treatment,
      season = season,
      period = period,
      metric = metric,

      parameter = rownames(
        sobol_res$T
      ),

      ST = sobol_res$T$original,

      stringsAsFactors = FALSE

    )

    sobol_list[[length(sobol_list)+1]] <- tmp
  }
}

sobol_df <- bind_rows(
  sobol_list
)

###############################################################################
# CLEAN
###############################################################################

sobol_df <- sobol_df %>%
  mutate(
    ST = pmax(ST,0)
  )

sobol_df <- sobol_df %>%
  group_by(
    treatment,
    season,
    period,
    metric
  ) %>%
  mutate(
    ST_norm = ST / sum(ST)
  ) %>%
  ungroup()

###############################################################################
# FACTOR ORDERING
###############################################################################

param_order <- c(
  "betad",
  "beta0",
  "omega",
  "Kd_v",
  "Kb_v",
  "omega_g_v",

  "e_forest",

  "beta_lw",
  "omega_lw",
  "Kd_lw_v",
  "omega_g_lw_v",

  "h",

  "g_macro_ref",
  "infl_macro_ref",

  "infl_soil_ref",
  "infl_forest_ref",

  "g_forest_ref",

  "p_ground",

  "g_soil_ref",
  "k_soil",

  "gs_max",
  "stomata_closure",

  "k_u",
  "u_ref"
)

metric_order <- c(
  "mean_understorey_T",
  "sd_understorey_T",
  "vertical_T_gradient",

  "mean_understorey_SW",
  "sd_understorey_SW",
  "vertical_SW_extinction"
)

metric_labels <- c(
  mean_understorey_T = "Mean understorey\ntemperature",

  sd_understorey_T = "SD understorey\ntemperature",

  `vertical_T_gradient` =
    "Vertical temperature\ngradient",

  mean_understorey_SW =
    "Mean understorey\nshortwave radiation",

  sd_understorey_SW =
    "SD understorey\nshortwave radiation",

  `vertical_SW_extinction` =
    "Vertical shortwave\nextinction"
)

sobol_df$metric_label <- factor(
  metric_labels[
    as.character(sobol_df$metric)
  ],
  levels = metric_labels[metric_order]
)

sobol_df$parameter <- factor(
  sobol_df$parameter,
  levels = param_order
)

sobol_df$metric <- factor(
  sobol_df$metric,
  levels = metric_order
)

sobol_df$condition <- factor(
  paste(
    sobol_df$season,
    sobol_df$period
  ),
  levels = c(
    "dry noon",
    "dry night",
    "wet noon",
    "wet night"
  )
)

###############################################################################
# TOP PARAMETERS EXPLAINING X% OF VARIANCE
###############################################################################

variance_threshold <- 0.50

temperature_metrics <- c(
  "mean_understorey_T",
  "sd_understorey_T",
  "vertical_T_gradient"
)

sw_metrics <- c(
  "mean_understorey_SW",
  "sd_understorey_SW",
  "vertical_SW_extinction"
)

###############################################################################
# TEMPERATURE PARAMETERS
###############################################################################

temp_importance <- sobol_df %>%

  filter(
    metric %in% temperature_metrics
  ) %>%

  group_by(parameter) %>%

  summarise(
    mean_ST = mean(ST_norm),
    .groups = "drop"
  ) %>%

  arrange(
    desc(mean_ST)
  )

temp_top_parameters <- c()

for(par in temp_importance$parameter){

  temp_top_parameters <- c(
    temp_top_parameters,
    par
  )

  coverage <- sobol_df %>%

    filter(
      metric %in% temperature_metrics
    ) %>%

    mutate(
      selected =
        parameter %in% temp_top_parameters
    ) %>%

    group_by(
      treatment,
      season,
      period,
      metric
    ) %>%

    summarise(
      coverage =
        sum(ST_norm[selected]),
      .groups = "drop"
    )

  if(min(coverage$coverage) >= variance_threshold){

    break

  }
}
###############################################################################
# SHORTWAVE PARAMETERS
###############################################################################

sw_importance <- sobol_df %>%

  filter(
    metric %in% sw_metrics
  ) %>%

  group_by(parameter) %>%

  summarise(
    mean_ST = mean(ST_norm),
    .groups = "drop"
  ) %>%

  arrange(
    desc(mean_ST)
  )

sw_top_parameters <- c()

for(par in sw_importance$parameter){

  sw_top_parameters <- c(
    sw_top_parameters,
    par
  )

  coverage <- sobol_df %>%

    filter(
      metric %in% sw_metrics
    ) %>%

    mutate(
      selected =
        parameter %in% sw_top_parameters
    ) %>%

    group_by(
      treatment,
      season,
      period,
      metric
    ) %>%

    summarise(
      coverage =
        sum(ST_norm[selected]),
      .groups = "drop"
    )

  if(min(coverage$coverage) >= variance_threshold){

    break

  }
}

write.csv(
  temp_importance,
  file.path(
    output_dir,
    "top_parameters_temperature.csv"
  ),
  row.names = FALSE
)

write.csv(
  sw_importance,
  file.path(
    output_dir,
    "top_parameters_shortwave.csv"
  ),
  row.names = FALSE
)



###############################################################################
# CREATE GROUPED PARAMETER COLUMN
###############################################################################

sobol_df <- sobol_df %>%

  mutate(

    parameter_plot = case_when(

      metric %in% temperature_metrics &
        parameter %in% temp_top_parameters
      ~ as.character(parameter),

      metric %in% sw_metrics &
        parameter %in% sw_top_parameters
      ~ as.character(parameter),

      TRUE
      ~ "Other parameters"

    )

  )

###############################################################################
# RECOMPUTE ST AFTER GROUPING
###############################################################################

sobol_df_plot <- sobol_df %>%

  group_by(
    treatment,
    season,
    period,
    metric,
    metric_label,
    condition,
    parameter_plot
  ) %>%

  summarise(
    ST_norm = sum(ST_norm),
    .groups = "drop"
  )

###############################################################################
# LEGEND ORDER
###############################################################################

legend_order <- c(

  temp_top_parameters,

  setdiff(
    sw_top_parameters,
    temp_top_parameters
  ),

  "Other parameters"

)

sobol_df_plot$parameter_plot <- factor(
  sobol_df_plot$parameter_plot,
  levels = legend_order
)

###############################################################################
# COLOURS
###############################################################################

temperature_only <- temp_top_parameters

shortwave_only <- sw_top_parameters

# temperature: blue to green
temp_cols <- colorRampPalette(
  c(
    "cornflowerblue",
    "royalblue",
    "chartreuse4",
    "darkgreen"

  )
)(
  max(1, length(temperature_only))
)

names(temp_cols) <- temperature_only

# shortwave:  yellow to orange
sw_cols <- colorRampPalette(
  c(    "gold2",
    "orange"
  )
)(
  max(1, length(shortwave_only))
)

names(sw_cols) <- shortwave_only

cols <- c(
  temp_cols,
  sw_cols,
  "Other parameters" = "grey80"
)

# keep legend order
cols <- cols[
  match(
    legend_order,
    names(cols)
  )
]

###############################################################################
# PARAMETER LABELS
###############################################################################

parameter_labels <- c(

  beta0 = expression(beta[0]),

  infl_macro_ref =
    expression(i[m]),

  k_soil =
    expression(k[s]),

  k_u =
    expression(k[u]),

  infl_soil_ref =
    expression(i[s]),

  u_ref =
    expression(u[r]),

  Kb_v =
    expression(K[b]),

  Kd_v =
    expression(K[d]),

  omega =
    expression(omega),

  omega_g_v =
    expression(omega[g]),

  betad =
    expression(beta[d]),

  h =
    expression(h),

  g_macro_ref =
    expression(g[m]),

  g_soil_ref =
    expression(g[s]),

  g_forest_ref =
    expression(g[f]),

  infl_forest_ref =
    expression(i[f]),

  gs_max =
    expression(g[s-max]),

  stomata_closure =
    expression(s[c]),

  p_ground =
    expression(p),

  e_forest =
    expression(e[f]),

  beta_lw =
    expression(beta[lw]),

  omega_lw =
    expression(omega[lw]),

  omega_g_lw_v =
    expression(omega[g-lw]),

  Kd_lw_v =
    expression(K[d-lw]),

  "Other parameters" =
    "Other parameters"
)

###############################################################################
# FIGURE 1
# STACKED BARPLOT
###############################################################################

plot_df <- sobol_df_plot

# treatment labels
plot_df <- plot_df %>%
  mutate(
    treatment = recode(
      treatment,
      control = "Control",
      removal = "Liana removal"
    )
  )

# drop unused factor levels
plot_df <- plot_df %>%
  mutate(
    condition = forcats::fct_drop(condition)
  )


p1 <- ggplot(
  plot_df,
  aes(
    x = ST_norm,
    y = condition,
    fill = parameter_plot
  )
) +

  geom_col(
    width = 0.60,
    colour = "white",
    linewidth = 0.3
  ) +

  facet_grid(
    metric_label ~ treatment,
    scales = "free_y",
    space = "free_y",
    switch = "y"
  ) +

  # y scale
  scale_y_discrete(
    position = "right",
    drop = TRUE
  ) +

  scale_fill_manual(
    values = cols,
    labels = parameter_labels,
    drop = FALSE
  ) +

  labs(
    x = "Normalised total-order Sobol index",
    y = NULL,
    fill = "Parameter"
  ) +

  theme_bw(base_size = 14) +

  theme(
    strip.text =
      element_text(
        size = 13,
        lineheight = 0.9
      ),

    strip.background =
      element_rect(
        fill = "grey90",
        colour = "black"
      ),

    legend.position = "right"
  ) +
  theme(

    axis.text.y =
      element_text(
        size = 11
      )

  ) +
  theme(

    axis.text.x =
      element_text(
        size = 12
      ),

    axis.title.x =
      element_text(
        size = 14
      )

  ) +


  theme(

    strip.placement = "outside",

    strip.text.y.left =
      element_text(
        angle = 0,
        size = 13
      )

  ) +

  theme(

    legend.title =
      element_text(
        size = 16
      ),

    legend.text =
      element_text(
        size = 14
      )

  ) +
  scale_x_continuous(
    breaks = c(0, 0.5, 1)
  ) +
  coord_cartesian(
    xlim = c(0, 1)
  )

ggsave(
  file.path(
    output_dir,
    "stacked_barplot_all_parameters.png"
  ),
  p1,
  width = 10,
  height = 5,
  dpi = 300
)


###############################################################################
# SAVE MASTER DATAFRAME
###############################################################################

saveRDS(
  sobol_df,
  file.path(
    output_dir,
    "sobol_master_dataframe.rds"
  )
)

###############################################################################
# CREATE TABLE WITH MEAN ABSOLUTE DELTA (CONTROL - REMOVAL) OVER ALL SCENARIOS
###############################################################################

delta_table <- sobol_df %>%
  select(treatment, metric, season, period, parameter, ST_norm) %>%
  pivot_wider(
    names_from = treatment,
    values_from = ST_norm
  ) %>%
  mutate(
    delta = control - removal
  ) %>%
  group_by(parameter) %>%
  summarise(
    Mean_Abs_Delta = mean(abs(delta), na.rm = TRUE),
    .groups = "drop"
  )%>%
  arrange(desc(Mean_Abs_Delta))

print(delta_table, n = 24)


###############################################################################
# median and emperical ranges on top parameters
###############################################################################

# Seven key parameters across temperature and shortwave metrics
top_7_parameters <- union(
  temp_top_parameters,
  sw_top_parameters
)

# Variance explained by these parameters in each scenario
top7_coverage <- sobol_df %>%

  filter(
    metric %in% c(temperature_metrics, sw_metrics)
  ) %>%

  group_by(
    treatment,
    season,
    period,
    metric
  ) %>%

  summarise(
    coverage = sum(
      ST_norm[parameter %in% top_7_parameters],
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# Summary across all scenarios
top7_summary <- top7_coverage %>%
  summarise(
    mean = mean(coverage),
    median = median(coverage),
    CI95_low = quantile(coverage, 0.025),
    CI95_high = quantile(coverage, 0.975),
    min = min(coverage),
    max = max(coverage)
  )

print(top7_summary)

###############################################################################
# FINISHED
###############################################################################

cat(
  "\nFinished.\n",
  nrow(sobol_df),
  "Sobol entries loaded.\n"
)
