###############################################################################
#
# In this script observed macro and micro air temperature plotted over time.
# The micro observations are done by TOMST TMS-4 sensors. 27 such sensors are
# positioned within each investigated Panama plot (C13 or R11).
#
###############################################################################

library(dplyr)
library(ggplot2)
library(tidyr)
library(lubridate)
library(dygraphs)
library(xts)
library(htmlwidgets)



#########
# INPUT #
#########

plotnr = "C13" #  C13 or R11

# time window (Panama time = EST, UTC-5)
t_start <- as.POSIXct("2023-01-01 00:00:00", tz = "EST")
t_end   <- as.POSIXct("2023-12-31 00:00:00", tz = "EST")



#########

# input file observations
input_temps <- readRDS("Data/TMS_macro_2023_C13_and_R11.rds")
# input file model
input_model <- readRDS(paste0("Output/timeseries/data/timeseries_", plotnr, ".rds"))

# mutate model outputs
model <- input_model %>%
  mutate(
    datetime_panama = with_tz(datetime, tzone = "EST"),
    datetime_hour   = floor_date(datetime_panama, unit = "hour")
  ) %>%
  group_by(datetime_hour) %>%
  summarise(
    Mean_model = mean(temperature_model, na.rm = TRUE),
    .groups = "drop"
  )

# filter and mutate observations
temps <- input_temps %>%
  filter(sensor_name == "TMS_T3",
         plot == plotnr) %>%
  filter(datetime_EST >= t_start & datetime_EST <= t_end)

df_mean <- temps %>%
  mutate(datetime_hour = floor_date(datetime_EST, unit = "hour")) %>%
  group_by(datetime_hour) %>%
  summarise(
    Mean_micro = mean(value, na.rm = TRUE),
    Macro      = mean(temp_macro, na.rm = TRUE),
    season     = names(which.max(table(season))),
    .groups = "drop"
  )

# merge model output with observations
df_mean <- df_mean %>%
  left_join(model, by = "datetime_hour")


season_rects <- df_mean %>%
  arrange(datetime_hour) %>%
  mutate(change = season != lag(season, default = first(season)),
         grp = cumsum(change)) %>%
  group_by(grp, season) %>%
  summarise(
    from = min(datetime_hour),
    to   = max(datetime_hour) + hours(1),
    .groups = "drop"
  )



# xts object for dygraphs
# xts object for dygraphs
x <- xts(df_mean[, c("Macro", "Mean_micro", "Mean_model")], order.by = df_mean$datetime_hour)
names(x) <- c("macro", "TOMST-mean micro", "model-mean micro")



p <- dygraph(x, main = paste("Timeseries of observed macro- & TOMST-mean microtemperature and model-mean microtemperature for plot", plotnr,
                             "<br><span style='font-size:12px'>
                          <span style='background-color:rgba(255,165,0,0.3);padding:2px 8px;'>dry season</span>
                          &nbsp;
                          <span style='background-color:rgba(173,216,230,0.3);padding:2px 8px;'>wet season</span>
                          </span>")) %>%
  dyRangeSelector() %>%                 # slider at the bottom
  dyOptions(drawPoints = FALSE, strokeWidth = 1.2) %>%
  dyAxis("y", label = "Temperature (°C)")
for(i in seq_len(nrow(season_rects))) {
  col <- if (season_rects$season[i] == "dry") "rgba(255,165,0,0.15)" else "rgba(173,216,230,0.15)"
  p <- p %>% dyShading(from = season_rects$from[i], to = season_rects$to[i], color = col)
}
#print(p)

saveWidget(
  p,
  file = paste0("Output/timeseries/plots/timeseries_slider_", plotnr, ".html"),
  selfcontained = TRUE
)
