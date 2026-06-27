############################################################################################
# This script combines plots.
#
# Author: Emma Van de Walle - Q-ForestLab
############################################################################################


library(magick)


# validation plots
#-----------------

# Read images
img1 <- image_read("Output/validation/plots/model_vs_obs_over_time_and_position_R11_both.png")
img2 <- image_read("Output/validation/plots/residual_by_hour_R11_both.png")

img3 <- image_read("Output/validation/plots/temporal_validation_R7_both.png")
img4 <- image_read("Output/validation/plots/spatial_validation_R7_both.png")

# Horizontal spacers
spacer_h <- image_blank(
  width = 500,
  height = image_info(img1)$height,
  color = "white"
)

# Create top row (a,b)
row1 <- image_append(
  c(img1, spacer_h, img2),
  stack = FALSE
)

# Create bottom row (c,d)
row2 <- image_append(
  c(img3, spacer_h, img4),
  stack = FALSE
)

# Vertical spacer between rows
spacer_v <- image_blank(
  width = image_info(row1)$width,
  height = 300,
  color = "white"
)

# Combine rows
combined <- image_append(
  c(row1, spacer_v, row2),
  stack = TRUE
)

# Save
image_write(
  combined,
  path = "Output/validation/plots/combined_validation_R11_R7.png",
  format = "png"
)


# calibration plots
#------------------
img_cal_rmse <- image_read("Output/calibration/plots/C13_50_generations_4_params.png")
img_cal_pca <- image_border(image_trim(image_read("Output/calibration/plots/PCA_space_C13_50_generations_4_params.png")), color = "white", geometry = "50x50")

spacer_cal <- image_blank(height = image_info(img_cal_pca)$height, width = 500, color = "white")

combined_cal <- image_append(c(img_cal_rmse, spacer_cal, img_cal_pca), stack = FALSE)

image_write(combined_cal, path = "Output/calibration/plots/combined_calibration_C13.png", format = "png")


# sensitivity plots
#------------------
img_SA_Tair <- image_read("Output/sensitivity_analysis/plots_output/Sobol_indices_C13_Ta.png")
img_SA_Tsurf <- image_read("Output/sensitivity_analysis/plots_output/Sobol_indices_C13_Tf.png")

spacer_SA <- image_blank(height = image_info(img_SA_Tair)$height, width = 500, color = "white")

combined_SA <- image_append(c(img_SA_Tair, spacer_SA, img_SA_Tsurf))

image_write(combined_SA, path = "Output/sensitivity_analysis/plots_output/combined_sensitivity.png", format = "png")


# flux/temperature plots
#-----------------------
img_SW <- image_read("Output/2023-07-08_14/SW_2D_F_b_down.png")
img_Tair <- image_read("Output/2023-07-08_14/temp_air.png")
img_Tsurf <- image_read("Output/2023-07-08_14/temp_surface.png")

spacer_tiles <- image_blank(height = 0, width = image_info(img_Tair)$width, color = "white")

combined_tiles <- image_append(c(img_SW, spacer_tiles, img_Tair, spacer_tiles, img_Tsurf), stack = TRUE)

image_write(combined_tiles, path = "Output/combined_tiles.png", format = "png")

# Shortwave RTM plots
#--------------------
img_direct <- image_read("Output/2023-07-08_14/SW_2D_F_b_down.png")
img_diffd <- image_read("Output/2023-07-08_14/SW_2D_F_d_down.png")
img_diffu <- image_read("Output/2023-07-08_14/SW_2D_F_d_up.png")

spacer_tiles <- image_blank(height = 0, width = image_info(img_direct)$width, color = "white")

combined_SW <- image_append(c(img_direct, spacer_tiles, img_diffd, spacer_tiles, img_diffu), stack = TRUE)

image_write(combined_SW, path = "Output/combined_SW.png", format = "png")

# Longwave RTM plots
#--------------------
img_Ld <- image_read("Output/2023-07-08_14/LW_2D_F_d_down.png")
img_Lu <- image_read("Output/2023-07-08_14/LW_2D_F_d_up.png")

spacer_tiles <- image_blank(height = 0, width = image_info(img_Ld)$width, color = "white")

combined_LW <- image_append(c(img_Ld, spacer_tiles, img_Lu), stack = TRUE)

image_write(combined_LW, path = "Output/combined_LW.png", format = "png")

# time series plots
#------------------
img_h <- image_read("Output/Tair_horizontal.png")
img_v <- image_read("Output/Tair_vertical.png")

spacer_tiles <- image_blank(height = image_info(img_h)$height, width = 150, color = "white")

combined_time <- image_append(c(img_h, spacer_tiles, img_v))

image_write(combined_time, path = "Output/combined_timeseries.png", format = "png")

# perfomance plots
#------------------
img_mo <- image_read("Output/timeseries/plots/Observed_vs_modelled_gradient.png")
img_bin <- image_read("Output/timeseries/plots/binned_performance.png")

spacer_tiles <- image_blank(height = image_info(img_mo)$height, width = 150, color = "white")

combined_time <- image_append(c(img_mo, spacer_tiles, img_bin))

image_write(combined_time, path = "Output/combined_performance.png", format = "png")


