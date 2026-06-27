library(dplyr)
library(ggplot2)

# load data
sensor_metrics <- readRDS(
  "Output/validation/data/sensor_metrics_C8_uncalibrated.rds"
)

val_data <- readRDS(
  "Output/validation/data/val_data_C8_uncalibrated.rds"
)

# unique sensor locations
sensor_locations <- val_data %>%
  group_by(sensor_id) %>%
  summarise(
    sensor_x = first(sensor_x),
    sensor_y = first(sensor_y),
    .groups = "drop"
  )

# combine
sensor_map <- left_join(
  sensor_locations,
  sensor_metrics,
  by = "sensor_id"
)


ggplot(sensor_map,
       aes(sensor_x,
           sensor_y,
           colour = R2)) +
  geom_point(size = 5) +
  geom_text(
    aes(label = sensor_id),
    nudge_y = 1
  ) +
  coord_equal() +
  scale_colour_viridis_c() +
  theme_bw()


ggplot(sensor_map,
       aes(sensor_x,
           sensor_y,
           colour = ME)) +
  geom_point(size = 5) +
  geom_text(
    aes(label = sensor_id),
    nudge_y = 1
  ) +
  coord_equal() +
  scale_colour_gradient2(
    midpoint = 0
  ) +
  theme_bw()


ggplot(sensor_map,
       aes(sensor_x,
           sensor_y,
           colour = ME,
           size = R2)) +
  geom_point() +
  geom_text(
    aes(label = sensor_id),
    nudge_y = 1
  ) +
  coord_equal() +
  scale_colour_gradient2(
    midpoint = 0
  ) +
  theme_bw()
