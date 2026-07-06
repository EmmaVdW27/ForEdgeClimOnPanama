#' Function to create several longwave radiation plots
#'
#' @param lw_rad Dataframe with longwave radiative values
#' @param output_path The output_path to store the plots
#' @return Several radiation plots (plotted and saved)
#' @importFrom dplyr filter group_by summarise
#' @importFrom scales rescale
#' @import ggplot2
#' @export
#'
plots_lw <- function(lw_rad, output_path){

  # 2D RTM plots
  ##############

  # Calculate averages across Y-slices
  final_avg_results <- lw_rad |>
    group_by(X, Z) |>
    summarise(avg_F_d_down = mean(F_d_down, na.rm = TRUE),
              avg_F_d_up = mean(F_d_up, na.rm = TRUE), .groups = 'keep')


  # F_d_down plot:
  F_d_down = ggplot(final_avg_results, aes(x = X, y = Z, fill = avg_F_d_down)) +
    geom_tile() +

    scale_fill_viridis_c(option = "inferno",#,  limits = c(300, 850), breaks = c(400, 600, 800),
                         guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = paste("Longwave downward radiation"),
         x = "Easting (m)", y = "Height (m)", fill = bquote("Flux"~(W~m^{-2})) ) +
    coord_fixed(ratio = 1) +
    theme_bw() #+
    # theme(
    #   plot.title = element_text(size = 20),
    #   plot.caption = element_text(hjust = 0, color = "black"),
    #   axis.title = element_text(size = 20),
    #   axis.text = element_text(size = 16),
    #   strip.text = element_text(size = 18),
    #   legend.title = element_text(size = 18),
    #   legend.text  = element_text(size = 16)
    # )

  print(F_d_down)

  ggsave(paste0(output_path, '/LW_F_d_down.png'), plot = F_d_down, width = 9, height = 3, dpi = 500)


  # F_d_up plot:
  F_d_up = ggplot(final_avg_results, aes(x = X, y = Z, fill = avg_F_d_up)) +
    geom_tile() +
    scale_x_continuous(
      breaks = seq(0, 150, by = 50),
      labels = rev(seq(0, 150, by = 50))
    ) +
    scale_fill_viridis_c(option = "inferno",# limits = c(300, 850), breaks = c(400, 600, 800),
                         guide = guide_colorbar(barwidth = 1, barheight = 10, frame.colour = "black", ticks.colour = "black")) +
    labs(title = paste("Longwave upward radiation"),
         x = "Easting (m)", y = "Height (m)", fill = bquote("Flux"~(W~m^{-2})) )+
    coord_fixed(ratio = 1) +
    theme_bw() #+
    # theme(
    #   plot.title = element_text(size = 20),
    #   plot.caption = element_text(hjust = 0, color = "black"),
    #   axis.title = element_text(size = 20),
    #   axis.text = element_text(size = 16),
    #   strip.text = element_text(size = 18),
    #   legend.title = element_text(size = 18),
    #   legend.text  = element_text(size = 16)
    # )

  print(F_d_up)

  ggsave(paste0(output_path, '/LW_F_d_up.png'), plot = F_d_up, width = 9, height = 3, dpi = 500)




}
