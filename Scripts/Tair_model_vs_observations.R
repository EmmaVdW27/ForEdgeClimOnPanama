###############################################################################
#
# In this script modelled and observed air temperature are compared.
# The observations are done by TOMST TMS-4 sensors. 27 such sensors are
# positioned within each investigated Panama plot (C13 or R11).
#
###############################################################################

library(dplyr)
library(ggplot2)
library(lubridate)
library(readxl)
library(stringr)
library(FNN)
library(Metrics)

#########
# INPUT #
#########

plotnr = "C13" #  C13 or R11

time_point <- ymd_h("2023-04-26 14")

req_height <- 1  # (m)

#########

if(plotnr == "C13"){
  model_input_file <- readRDS("Output/2023-04-26_14/C13/model_results.rds")

  tomst_input_file <- read.csv("Data/microtemp_air_obs_20230426_14_C13.csv")
} else { # plotnr == "R11"
  model_input_file <- readRDS("Output/2023-04-26_14/R11/model_results.rds")

  tomst_input_file <- read.csv("Data/microtemp_air_obs_20230426_14_R11.csv")
}


tomst_coordinates_file <- read_excel("Data/Data_Kasper/calculated_locations_all_plots_Emma.xlsx")

####################
# DATA  EXTRACTION #
####################

# Model data
micro_grid <- model_input_file$micro_grid
air_temp   <- model_input_file$air_temperature

model_xy <- micro_grid %>%
  mutate(T_model = air_temp - 273.15) %>%
  filter(z == req_height) %>%
  select(x, y, T_model)


# TOMST data
tomst <- tomst_input_file %>%
transmute(
  SN = str_trim(as.character(SN)),
  T_obs = as.numeric(value_mean)
)
locs <- tomst_coordinates_file %>%
  transmute(
    SN = str_trim(as.character(SN)),
    x  = as.numeric(x),
    y  = as.numeric(y)
  )
tomst_xy <- tomst %>%
  left_join(locs, by = "SN")

#####################################
# COMPARING MODEL WITH OBSERVATIONS #
#####################################
# Nearest neighbour match TOMST -> model cell
nn <- get.knnx(
  data  = model_xy[, c("x","y")],
  query = tomst_xy[, c("x","y")],
  k = 1 # search for exactly 1 closest voxel per TOMST sensor
)

# Combine model with TOMST values
comp <- tomst_xy %>%
  mutate(
    nn_dist  = nn$nn.dist[,1],
    T_model  = model_xy$T_model[nn$nn.index[,1]],
    residual = T_model - T_obs
  )

# Metrics
rmse_val <- rmse(comp$T_obs, comp$T_model)
mae_val  <- mae(comp$T_obs, comp$T_model)
bias_val <- mean(comp$residual)
r_val    <- cor(comp$T_obs, comp$T_model)

cat(
  "RMSE:", round(rmse_val, 2), "°C | ",
  "MAE:",  round(mae_val,  2), "°C | ",
  "Bias:", round(bias_val, 2), "°C | ",
  "R:",    round(r_val,    2), "\n"
)

#########
# PLOTS #
#########

# Spatial map: model temperature field + TOMST points (black)
p_map <- ggplot(model_xy, aes(x, y, fill = T_model)) +
  geom_raster() +
  geom_point(data = comp, aes(x, y, shape = "TOMST sensor"), colour = "black", size = 2) +
  coord_equal() +
  theme_bw() +
  scale_fill_viridis_c(name = "Model (°C)") +
  scale_shape_manual(name   = NULL, values = c("TOMST sensor" = 16)) +
  labs(title = "Modelled air temperature field at 1 m height", x = "Easting (m)", y = "Northing (m)") +
  guides(
    fill = guide_colourbar(
      frame.colour = "black",
      frame.linewidth = 1,
      ticks.colour = "black"
    ),
    shape = guide_legend(override.aes = list(colour = "black"))
  )

print(p_map)

ggsave(paste("Output/modelled_Tair_1m_height_", plotnr, ".png"), p_map, width = 5, height = 5, dpi = 300)



# Residual map
p_resid <- ggplot(comp, aes(x, y, colour = residual)) +
  geom_point(size = 3) +
  coord_equal() +
  theme_bw() +
  scale_colour_viridis_c(name   = "Residual (°C)") +
  labs(title = "Residuals (Model − observations)", x = "Easting (m)", y = "Northing (m)") +
  guides(
    colour = guide_colourbar(
      frame.colour   = "black",
      frame.linewidth= 1,
      ticks.colour   = "black"
    )
  )

print(p_resid)

ggsave(paste("Output/Residu_Tair_model_vs_obs_", plotnr, ".png"), p_resid, width = 5, height = 5, dpi = 300)

# Scatter: model vs obs
p_scatter <- ggplot(comp, aes(T_obs, T_model)) +
  geom_point(size = 3) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  # 1:1 line
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    colour = "cornflowerblue",
    linewidth = 1
  ) +
  theme_bw() +
  labs(title = "Model vs TOMST", x = "TOMST observation (°C)", y = "Model (°C)") +
  annotate("text",
           x = min(comp$T_obs, na.rm=TRUE),
           y = max(comp$T_model, na.rm=TRUE),
           hjust = 0,
           label = paste0("RMSE=", round(rmse_val,2),
                          " | Bias=", round(bias_val,2),
                          " | R=", round(r_val,2))) +
  # 1:1 label
  annotate(
    "text",
    x = min(comp$T_obs, na.rm = TRUE) +
      0.05 * diff(range(comp$T_obs, na.rm = TRUE)),
    y = min(comp$T_obs, na.rm = TRUE) +
      0.05 * diff(range(comp$T_obs, na.rm = TRUE)),
    label = "1:1",
    colour = "cornflowerblue"
  )
print(p_scatter)

ggsave(paste("Output/Scatter_Tair_model_vs_obs_", plotnr, ".png"), p_scatter, width = 5, height = 5, dpi = 300)

