############################################################################################
# Combine model output figures
#
# This script combines individual validation and calibration plots into composite figures.
#
# Author: Emma Van de Walle - Q-ForestLab
############################################################################################


library(magick)


# validation plots
#-----------------

# Read images
img1 <- image_read("Output/validation/plots/model_vs_obs_over_time_and_position_C13_both.png")
img2 <- image_read("Output/validation/plots/residual_by_hour_C13_both.png")

img3 <- image_read("Output/validation/plots/temporal_validation_C8_both.png")
img4 <- image_read("Output/validation/plots/spatial_validation_C8_both.png")

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
  path = "Output/validation/plots/combined_validation_C13_C8.png",
  format = "png"
)


# calibration plots
#------------------
img_cal_rmse <- image_read("Output/calibration/plots/r11_50_generations_5_params.png")
img_cal_pca <- image_border(image_trim(image_read("Output/calibration/plots/PCA_space_R11_50_generations_5_params.png")), color = "white", geometry = "50x50")

spacer_cal <- image_blank(height = image_info(img_cal_pca)$height, width = 500, color = "white")

combined_cal <- image_append(c(img_cal_rmse, spacer_cal, img_cal_pca), stack = FALSE)

image_write(combined_cal, path = "Output/calibration/plots/combined_calibration_R11.png", format = "png")


