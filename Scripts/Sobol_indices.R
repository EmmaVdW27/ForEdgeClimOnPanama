###############################################################################
# In this script all 'Sobol.R' rds outputfiles for the several
# time points are merged together to get an overview of the
# Quantities of Interest (3 time points x 3 metrics).
# This script analyzes the air temperature. Sobol_indices_forest_surfaceT.R is similar,
# but for the forest surface temperature.
#
# Author: Emma Van de Walle - Q-ForestLab
#
###############################################################################


library(dplyr)
library(ggplot2)
library(readr)
library(purrr)
library(tidyr)
library(forcats)
library(latex2exp)
library(glue)

#########
# INPUT #
#########

plots = c("C13", "R11")        # C13 or R11 ~ control or removal plot
temps = c("a", "f")      # a or f ~ air or forest surface temperature

#########

# paths:

input_path <- "Output/sensitivity_analysis/data/"
output_path_plots <- "Output/sensitivity_analysis/plots_output/"
output_path_numbers <- "Output/sensitivity_analysis/numbers_output/"

for(plotnr in plots){
  for(temp_type in temps){
    output_plot_focused_parameters_by_condition_normalized <- glue("Sobol_indices_{plotnr}_T{temp_type}.png")

    if (temp_type == "a"){
      plot_title = "air"
    } else{
      plot_title = "forest surface"
    }


    # File list with labels
    # date: 26042023 or 26112023 ~ end of dry or wet season in Panama
    # hour of the day in UTC: 19 or 7 ~ noon or night in Panama (EST)
    files_info <- tribble(
      ~file, ~label,
      ## Dry season
      ### Noon
      glue("400samples_24parameters_19h_26042023_{plotnr}_avT{temp_type}.rds"),     glue("Dr_No_≈T{temp_type}"),
      glue("400samples_24parameters_19h_26042023_{plotnr}_SDT{temp_type}.rds"),     glue("Dr_No_σT{temp_type}"),
      glue("400samples_24parameters_19h_26042023_{plotnr}_avT{temp_type}_v.rds"),     glue("Dr_No_≈T{temp_type}_v"),
      glue("400samples_24parameters_19h_26042023_{plotnr}_SDT{temp_type}_v.rds"),     glue("Dr_No_σT{temp_type}_v"),
      glue("400samples_24parameters_19h_26042023_{plotnr}_gradT{temp_type}_v.rds"),   glue("Dr_No_∇T{temp_type}_v"),
      ### Night
      glue("400samples_24parameters_07h_26042023_{plotnr}_avT{temp_type}.rds"),     glue("Dr_Ni_≈T{temp_type}"),
      glue("400samples_24parameters_07h_26042023_{plotnr}_SDT{temp_type}.rds"),     glue("Dr_Ni_σT{temp_type}"),
      glue("400samples_24parameters_07h_26042023_{plotnr}_avT{temp_type}_v.rds"),     glue("Dr_Ni_≈T{temp_type}_v"),
      glue("400samples_24parameters_07h_26042023_{plotnr}_SDT{temp_type}_v.rds"),     glue("Dr_Ni_σT{temp_type}_v"),
      glue("400samples_24parameters_07h_26042023_{plotnr}_gradT{temp_type}_v.rds"),   glue("Dr_Ni_∇T{temp_type}_v"),
      ## Wet season
      ### Noon
      glue("400samples_24parameters_19h_26112023_{plotnr}_avT{temp_type}.rds"),     glue("We_No_≈T{temp_type}"),
      glue("400samples_24parameters_19h_26112023_{plotnr}_SDT{temp_type}.rds"),     glue("We_No_σT{temp_type}"),
      glue("400samples_24parameters_19h_26112023_{plotnr}_avT{temp_type}_v.rds"),     glue("We_No_≈T{temp_type}_v"),
      glue("400samples_24parameters_19h_26112023_{plotnr}_SDT{temp_type}_v.rds"),     glue("We_No_σT{temp_type}_v"),
      glue("400samples_24parameters_19h_26112023_{plotnr}_gradT{temp_type}_v.rds"),   glue("We_No_∇T{temp_type}_v"),
      ### Night
      glue("400samples_24parameters_07h_26112023_{plotnr}_avT{temp_type}.rds"),     glue("We_Ni_≈T{temp_type}"),
      glue("400samples_24parameters_07h_26112023_{plotnr}_SDT{temp_type}.rds"),     glue("We_Ni_σT{temp_type}"),
      glue("400samples_24parameters_07h_26112023_{plotnr}_avT{temp_type}_v.rds"),     glue("We_Ni_≈T{temp_type}_v"),
      glue("400samples_24parameters_07h_26112023_{plotnr}_SDT{temp_type}_v.rds"),     glue("We_Ni_σT{temp_type}_v"),
      glue("400samples_24parameters_07h_26112023_{plotnr}_gradT{temp_type}_v.rds"),   glue("We_Ni_∇T{temp_type}_v"),
    )

    ##########
    # OUTPUT #
    ##########

    # Parameter sequence
    param_order <- c(
      "betad", "beta0", "omega", "Kd_v", "Kb_v", "omega_g_v",
      "e_forest", "beta_lw", "omega_lw", "Kd_lw_v", "omega_g_lw_v",
      "h", "g_macro_ref", "infl_macro_ref", "infl_soil_ref", "infl_forest_ref", "g_forest_ref", "p_ground", "g_soil_ref", "k_soil",
      "gs_max", "stomata_closure", "k_u", "u_ref"
    )
    # Focused parameter sequence in order from most influential to least influential (to be plot in this order)
    focus_params <- c("infl_macro_ref", "infl_soil_ref", "g_macro_ref",
                      "g_soil", "infl_forest_ref", "k_soil",
                       "h", "p_ground", "g_forest_ref",
                      "gs_max", "stomata_closure", "k_u", "u_ref")
    # Shortwave parameter sequence
    SW_params <- c("betad", "beta0", "omega", "Kd_v", "Kb_v", "omega_g_v")
    # Longwave parameter sequence
    LW_params <- c("e_forest", "beta_lw", "omega_lw", "Kd_lw_v", "omega_g_lw_v")

    focus_colors <- c(
      "#3B4CC0",  # deep blue
      "#5E7FDB",  # blue
      "#F4A259",  # orange
      "#F7C58A",  # light orange
      "#87AFDE",  # light blue
      "#D62828",  # red
      "#6A040F",  # very dark red
      "#A4161A",  # dark red
      "#E98B1E",  # deep orange
      "darkgreen", # dark green
      "lightgreen",
      "purple",   # purple
      "pink"     # pink
    )



    # Color palette
    all_colors <- c(
      setNames(focus_colors, focus_params),
      "SW" = "darkgray",
      "LW" = "black"
    )

    # legend labels
    labels_expr <- c(
      infl_macro_ref  = expression(italic(i)[m]),
      infl_soil_ref   = expression(italic(i)[s]),
      k_soil      = expression(italic(k)[s]),
      g_macro_ref     = expression(italic(g)[m]),
      infl_forest_ref = expression(italic(i)[f]),
      g_soil_ref      = expression(italic(g)[s]),
      h           = expression(italic(h)),
      p_ground    = expression(italic(p)),
      g_forest_ref    = expression(italic(g)[f]),
      gs_max    = expression(italic(g)[s * "," * max]),
      stomata_closure     = expression(italic(s)[c]),
      k_u     = expression(italic(k)[u]),
      u_ref     = expression(italic(u)[r]),
      `SW parameters` = "SW parameters",
      `LW parameters` = "LW parameters"
    )

    scale_fill_manual(values = all_colors, labels = labels_expr)


    # Load, combine and normalize
    sobol_df <- files_info %>%
      mutate(data = map(file, ~ readRDS(file.path(input_path, .x)))) %>%
      mutate(S = map(data, ~ .x$T$original)) %>%
      select(label, S) %>%
      unnest(S) %>%
      mutate(parameter = rep(param_order, times = nrow(.) / length(param_order))) %>%
      group_by(label) %>%
      mutate(norm_value = pmax(0, S) / sum(pmax(0, S))) %>%
      ungroup() %>%
      mutate(parameter = factor(parameter, levels = param_order),
             label = fct_rev(fct_inorder(label))) %>%
      # define conditions
      mutate(
        season = substr(label, 1, 2),
        moment = substr(label, 4, 5),
        metric = case_when(
          grepl("≈T.*_v", label) ~ "avT_v",
          grepl("σT.*_v", label) ~ "sigmaT_v",
          grepl("∇T.*_v", label) ~ "gradT_v",
          grepl("≈T", label) ~ "avT",
          grepl("σT", label) ~ "sigmaT",
          TRUE ~ NA_character_
        )
      )


    ###############################################################################
    # PLOTTING RESULTS
    ###############################################################################

    # Plot mean parameter contribution per condition (metric, moment, season)

    # Dataframe per condition
    # Mean per metric
    mean_metric_n <- sobol_df %>%
      filter(!is.na(metric)) %>%
      group_by(metric, parameter) %>%
      summarise(mean_index = mean(norm_value), .groups = "drop") %>%
      rename(condition = metric)

    # Mean per moment
    mean_moment_n <- sobol_df %>%
      filter(!is.na(moment)) %>%
      group_by(moment, parameter) %>%
      summarise(mean_index = mean(norm_value), .groups = "drop") %>%
      rename(condition = moment)

    # Mean per season
    mean_season_n <- sobol_df %>%
      filter(!is.na(season)) %>%
      group_by(season, parameter) %>%
      summarise(mean_index = mean(norm_value), .groups = "drop") %>%
      rename(condition = season)

    # Combine dataframes
    mean_combined_n <- bind_rows(
      mutate(mean_metric_n, group = "Metric"),
      mutate(mean_moment_n, group = "Moment"),
      mutate(mean_season_n, group = "Season")
    ) %>%
      mutate(
        condition = as.character(condition),
        condition = recode(condition,
                           "avT" = "≈T",
                           "sigmaT" = "σT",
                           "gradT_v" = "∇T_v",
                           "avT_v" = "≈T_v",
                           "sigmaT_v" = "σT_v"

        ),
        # moments en seasons voluit
        condition = case_when(
          group == "Moment" ~ recode(condition,
                                     "Ni" = "Night",
                                     "No" = "Noon"),
          group == "Season" ~ recode(condition,
                                     "Dr" = "Dry",
                                     "We" = "Wet"),
          TRUE ~ condition
        ),
        condition = factor(condition, levels = rev(c("≈T", "σT", "≈T_v", "σT_v", "∇T_v",
                                                     "Night", "Noon",
                                                     "Dry","Wet"))),
        group = factor(group, levels = c("Metric", "Moment", "Season")),
        parameter = factor(parameter, levels = param_order)
      )



    mean_combined_focusedplus_n <- mean_combined_n %>%
      mutate(grouped_param = case_when(
        parameter %in% focus_params ~ as.character(parameter),
        parameter %in% SW_params   ~ "SW",
        parameter %in% LW_params   ~ "LW",
        TRUE ~ NA_character_
      )) %>%
      filter(!is.na(grouped_param)) %>%
      group_by(group, condition, grouped_param) %>%
      summarise(mean_index = sum(mean_index), .groups = "drop") %>%
      mutate(
        grouped_param = factor(grouped_param, levels = rev(c(focus_params, "SW", "LW")))
      )

    p_condition_focused_n <- ggplot(mean_combined_focusedplus_n, aes(x = mean_index, y = condition, fill = grouped_param)) +
      geom_col(width = 0.5) +
      facet_wrap(~ group, scales = "free_y", ncol = 1) +
      scale_fill_manual(values = all_colors, labels = labels_expr) +
      labs(
        title = glue("Influence on {plot_title} temperature in {plotnr}"),
        x = "Average normalised total-order Sobol-index",
        y = "Condition",
        fill = "Parameter"
      ) +
      theme_bw(base_size = 14) +
      guides(fill = guide_legend(ncol = 1, reverse = TRUE)) +
      theme(
        plot.title = element_text(size = 40),
        axis.title    = element_text(size = 40),
        axis.text     = element_text(size = 40),
        legend.position = "right",
        legend.direction = "vertical",
        legend.title    = element_text(size = 40),
        legend.text     = element_text(size = 40),
        legend.key.size = unit(1.0, "cm"),
        legend.spacing.x = unit(0.5, "cm"),
        strip.text = element_text(size = 40),
        plot.caption = element_text(hjust = 0)
      )

    ggsave(paste0(output_path_plots, output_plot_focused_parameters_by_condition_normalized), plot = p_condition_focused_n, width = 16, height = 14, dpi = 500)


    ###############################################################################
    # EXTRA OUTPUT: per conditie % variance + bereik per parameter + top-set tot 65%
    ###############################################################################

    # A) Per conditie: % variance per parameter (aflopend)
    sobol_percent_per_condition <- mean_combined_n %>%
      group_by(group, condition) %>%
      mutate(mean_index_pos = pmax(0, mean_index)) %>%
      mutate(
        norm_in_condition = mean_index_pos / sum(mean_index_pos),
        percent = 100 * norm_in_condition
      ) %>%
      arrange(desc(percent), .by_group = TRUE) %>%
      mutate(
        cum_percent = cumsum(percent),
        rank = row_number()
      ) %>%
      ungroup() %>%
      select(group, condition, rank, parameter, mean_index, percent, cum_percent)

    write_csv(
      sobol_percent_per_condition,
      file.path(output_path_numbers, glue("Sobol_percent_variance_per_condition_{plotnr}_T{temp_type}.csv"))
    )

    # B) Bereik per parameter over alle condities (min–max + optioneel mean/median)
    sobol_range_per_parameter <- sobol_percent_per_condition %>%
      group_by(parameter) %>%
      summarise(
        min_percent = min(percent, na.rm = TRUE),
        max_percent = max(percent, na.rm = TRUE),
        range_percent = max_percent - min_percent,
        median_percent = median(percent, na.rm = TRUE),
        mean_percent = mean(percent, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(range_percent))

    write_csv(
      sobol_range_per_parameter,
      file.path(output_path_numbers, glue("Sobol_percent_range_per_parameter_over_conditions_{plotnr}_T{temp_type}.csv"))
    )

    # C) Per conditie: minimale parameter-set die cumulatief ≥ 65% verklaart
    sobol_top65_per_condition_long <- sobol_percent_per_condition %>%
      group_by(group, condition) %>%
      mutate(cut_rank_65 = which(cum_percent >= 65)[1]) %>%
      mutate(cut_rank_65 = ifelse(is.na(cut_rank_65), n(), cut_rank_65)) %>%
      filter(rank <= cut_rank_65) %>%
      ungroup()

    # Samenvatting per conditie (één rij per conditie)
    sobol_top65_per_condition_summary <- sobol_top65_per_condition_long %>%
      group_by(group, condition) %>%
      summarise(
        n_params_65 = max(rank),
        explained_percent_65 = max(cum_percent),
        params_65 = paste(as.character(parameter), collapse = ", "),
        .groups = "drop"
      )

    write_csv(
      sobol_top65_per_condition_summary,
      file.path(output_path_numbers, glue("Sobol_parameters_explaining_ge65pct_per_condition_{plotnr}_T{temp_type}.csv"))
    )

    # D) Over alle condities: hoe "stabiel" zit een parameter in de 65%-topset?
    #    (hoe vaak komt hij voor, en rank-bereik als hij voorkomt)
    n_conditions <- sobol_percent_per_condition %>%
      distinct(group, condition) %>%
      nrow()

    sobol_top65_parameter_stability <- sobol_top65_per_condition_long %>%
      group_by(parameter) %>%
      summarise(
        n_conditions_in_top65 = n_distinct(paste(group, condition, sep = " / ")),
        frac_conditions_in_top65 = n_conditions_in_top65 / n_conditions,
        min_rank_in_top65 = min(rank, na.rm = TRUE),
        max_rank_in_top65 = max(rank, na.rm = TRUE),
        median_rank_in_top65 = median(rank, na.rm = TRUE),
        .groups = "drop"
      ) %>%
      arrange(desc(n_conditions_in_top65), median_rank_in_top65)

    write_csv(
      sobol_top65_parameter_stability,
      file.path(output_path_numbers, glue("Sobol_top65_parameter_stability_over_conditions_{plotnr}_T{temp_type}.csv"))
    )

    message(glue(
      "Saved 4 files: per-condition %, range per parameter, per-condition ≥65% set, and stability over conditions."
    ))


    # ###############################################################################
    # # QUANTITATIVE RESULTS (for normalized Sobol analysis)
    # ###############################################################################
    #
    # # Top 3 parameter contributors per condition
    #
    # top3_per_label <- sobol_df %>%
    #   group_by(label) %>%
    #   arrange(desc(norm_value), .by_group = TRUE) %>%
    #   mutate(rank = row_number()) %>%
    #   filter(rank <= 3) %>%
    #   summarise(
    #     top1_param = parameter[1], top1_val = norm_value[1],
    #     top2_param = parameter[2], top2_val = norm_value[2],
    #     top3_param = parameter[3], top3_val = norm_value[3],
    #     top3_sum   = sum(norm_value[1:3], na.rm = TRUE),
    #     .groups = "drop"
    #   )
    #
    # # Herfindahl–Hirschman (concentration) and entropy as measure for spread of influence
    #
    # concentration_per_label <- sobol_df %>%
    #   group_by(label) %>%
    #   summarise(
    #     # HHI: square each normalized value and sum them
    #     # If all influence is concentrated in one category (monoculture), norm_value = 1 for one, 0 for rest => HHI = 1
    #     # If influence is evenly distributed across n categories: HHI = 1/n (so the lower the value, the greater the spread)
    #     HHI = sum(norm_value^2),
    #     # Shannon entropy measures uncertainty/spread (-sum(norm_value * log(norm_value)))
    #     # Maximum entropy occurs when all categories are equal in size (maximum spread)
    #     # Entropy = 0 when one category holds 100% of the influence (no spread)
    #     entropy = -sum(ifelse(norm_value>0, norm_value*log(norm_value), 0)),
    #     .groups = "drop"
    #   )
    #
    # # Define to each parameter a process group
    # param_sets <- tribble(
    #   ~parameter,        ~process,
    #   "betad","SW","beta0","SW","omega","SW","Kd_v","SW","Kb_v","SW","omega_g_v","SW",
    #   "e_forest","LW","beta_lw","LW","omega_lw","LW","Kd_lw_v","LW","omega_g_lw_v","LW",
    #   "h","HEAT","g_macro","HEAT","infl_macro","HEAT","infl_soil","HEAT","infl_forest","HEAT","g_forest","HEAT","p_ground","HEAT","g_soil","HEAT","k_soil","HEAT"
    # ) %>% mutate(process = factor(process, levels = c("SW","LW","HEAT")))
    #
    # sobol_df2 <- sobol_df %>%
    #   left_join(param_sets, by = "parameter")
    #
    # # Share per process per condition
    #
    # process_share_per_label <- sobol_df2 %>%
    #   group_by(label, process) %>%
    #   summarise(process_share = sum(norm_value), .groups = "drop")
    #
    #
    # #------------------------------------------------------------------------------
    #
    # # Mean and sd parameter contribution per metric, moment and season
    # avg_param_by_metric  <- sobol_df %>% filter(!is.na(metric)) %>%
    #   group_by(metric, parameter) %>%
    #   summarise(mean_index = mean(norm_value), sd_index = sd(norm_value), .groups = "drop")
    #
    # avg_param_by_moment  <- sobol_df %>% filter(!is.na(moment)) %>%
    #   group_by(moment, parameter) %>%
    #   summarise(mean_index = mean(norm_value), sd_index = sd(norm_value), .groups = "drop")
    #
    # avg_param_by_season  <- sobol_df %>% filter(!is.na(season)) %>%
    #   group_by(season, parameter) %>%
    #   summarise(mean_index = mean(norm_value), sd_index = sd(norm_value), .groups = "drop")
    #
    # # get ranges on the contributions for top 3 parameters across conditions
    # target_params <- c("infl_soil" , "infl_macro", "g_macro")
    # # SW parameters:
    # #target_params <- c("betad", "beta0", "omega", "Kd_v", "Kb_v", "omega_g_v")
    # # LW parameters:
    # #target_params <- c("e_forest", "beta_lw", "omega_lw", "Kd_lw_v", "omega_g_lw_v")
    # # rad params
    # # target_params <- c("betad", "beta0", "omega", "Kd_v", "Kb_v", "omega_g_v",
    # #                    "e_forest", "beta_lw", "omega_lw", "Kd_lw_v", "omega_g_lw_v")
    #
    # range_metric <- avg_param_by_metric %>%
    #   filter(parameter %in% target_params) %>%
    #   group_by(metric) %>%
    #   summarise(sum_share = sum(mean_index), .groups = "drop") %>%
    #   summarise(min_sum = min(sum_share),
    #             max_sum = max(sum_share))
    #
    # range_moment <- avg_param_by_moment %>%
    #   filter(parameter %in% target_params) %>%
    #   group_by(moment) %>%
    #   summarise(sum_share = sum(mean_index), .groups = "drop") %>%
    #   summarise(min_sum = min(sum_share),
    #             max_sum = max(sum_share))
    #
    # range_season <- avg_param_by_season %>%
    #   filter(parameter %in% target_params) %>%
    #   group_by(season) %>%
    #   summarise(sum_share = sum(mean_index), .groups = "drop") %>%
    #   summarise(min_sum = min(sum_share),
    #             max_sum = max(sum_share))
    #
    # print(range_metric)
    # print(range_moment)
    # print(range_season)
    #
    #
    # # Mean and sd process contribution per metric, moment and season
    # process_by_label <- sobol_df2 %>%
    #   group_by(label, metric, moment, season, process) %>%
    #   summarise(process_share = sum(norm_value), .groups = "drop")
    #
    # avg_process_by_metric <- process_by_label %>% filter(!is.na(metric)) %>%
    #   group_by(metric, process) %>%
    #   summarise(mean_share = mean(process_share), sd_share   = sd(process_share), .groups = "drop")
    #
    # avg_process_by_moment <- process_by_label %>% filter(!is.na(moment)) %>%
    #   group_by(moment, process) %>%
    #   summarise(mean_share = mean(process_share), sd_share   = sd(process_share), .groups = "drop")
    #
    # avg_process_by_season <- process_by_label %>% filter(!is.na(season)) %>%
    #   group_by(season, process) %>%
    #   summarise(mean_share = mean(process_share), sd_share   = sd(process_share), .groups = "drop")
    #
    # #------------------------------------------------------------------------------
    #
    # # Stability : rank correlation and variation
    # # Rank per condition (1 = most important)
    # ranks_long <- sobol_df %>%
    #   group_by(label) %>%
    #   mutate(rank = dense_rank(desc(norm_value))) %>%
    #   ungroup()
    #
    # # Consistency: mean rang + sd per parameter for every condition
    # rank_consistency_moment <- ranks_long %>%
    #   group_by(moment, parameter) %>%
    #   summarise(mean_rank = mean(rank), sd_rank = sd(rank), .groups = "drop")
    #
    # rank_consistency_season <- ranks_long %>%
    #   group_by(season, parameter) %>%
    #   summarise(mean_rank = mean(rank), sd_rank = sd(rank), .groups = "drop")
    #
    # rank_consistency_metric <- ranks_long %>%
    #   group_by(metric, parameter) %>%
    #   summarise(mean_rank = mean(rank), sd_rank = sd(rank), .groups = "drop")
    #
    # # get ranges on the rank and sd on rank for top 3 parameters across conditions
    # # Moment
    # range_moment_rank <- rank_consistency_moment %>%
    #   filter(parameter %in% target_params) %>%
    #   summarise(
    #     min_mean_rank = min(mean_rank),
    #     max_mean_rank = max(mean_rank),
    #     min_sd_rank   = min(sd_rank),
    #     max_sd_rank   = max(sd_rank)
    #   )
    #
    # # Season
    # range_season_rank <- rank_consistency_season %>%
    #   filter(parameter %in% target_params) %>%
    #   summarise(
    #     min_mean_rank = min(mean_rank),
    #     max_mean_rank = max(mean_rank),
    #     min_sd_rank   = min(sd_rank),
    #     max_sd_rank   = max(sd_rank)
    #   )
    #
    # # Metric
    # range_metric_rank <- rank_consistency_metric %>%
    #   filter(parameter %in% target_params) %>%
    #   summarise(
    #     min_mean_rank = min(mean_rank),
    #     max_mean_rank = max(mean_rank),
    #     min_sd_rank   = min(sd_rank),
    #     max_sd_rank   = max(sd_rank)
    #   )
    #
    # # print(range_metric_rank)
    # # print(range_moment_rank)
    # # print(range_season_rank)
    #
    # #------------------------------------------------------------------------------
    #
    # # export tables
    # write_csv(top3_per_label,          file.path(output_path_numbers, "top3_per_condition_ah.csv"))
    # write_csv(concentration_per_label, file.path(output_path_numbers, "concentration_per_condition_ah.csv"))
    # write_csv(process_share_per_label, file.path(output_path_numbers, "process_share_per_condition_ah.csv"))
    #
    # write_csv(avg_param_by_metric,     file.path(output_path_numbers, "avg_param_by_metric_ah.csv"))
    # write_csv(avg_param_by_moment,     file.path(output_path_numbers, "avg_param_by_moment_ah.csv"))
    # write_csv(avg_param_by_season,     file.path(output_path_numbers, "avg_param_by_season_ah.csv"))
    #
    # write_csv(avg_process_by_metric,   file.path(output_path_numbers, "avg_process_by_metric_ah.csv"))
    # write_csv(avg_process_by_moment,   file.path(output_path_numbers, "avg_process_by_moment_ah.csv"))
    # write_csv(avg_process_by_season,   file.path(output_path_numbers, "avg_process_by_season_ah.csv"))
    #
    # write_csv(rank_consistency_metric, file.path(output_path_numbers, "rank_consistency_metric_ah.csv"))
    # write_csv(rank_consistency_moment, file.path(output_path_numbers, "rank_consistency_moment_ah.csv"))
    # write_csv(rank_consistency_season, file.path(output_path_numbers, "rank_consistency_season_ah.csv"))
  }
}




