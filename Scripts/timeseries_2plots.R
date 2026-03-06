###############################################################################
# Combined dygraph for C13 + R11
# Macro is identical -> 1 curve
# Micro + model differ per plot -> separate curves for C13 and R11
###############################################################################

library(dplyr)
library(lubridate)
library(dygraphs)
library(xts)
library(htmlwidgets)

#########
# INPUT #
#########

plots <- c("C13", "R11")

# time window (Panama time = EST, UTC-5)
t_start <- as.POSIXct("2023-01-01 00:00:00", tz = "EST")
t_end   <- as.POSIXct("2023-12-31 00:00:00", tz = "EST")

# observations (contains both plots)
input_temps <- readRDS("Data/TMS_macro_2023_C13_and_R11.rds")

# model files (one per plot)
input_model_C13 <- readRDS("Output/timeseries/data/timeseries_C13_calibrated.rds")
input_model_R11 <- readRDS("Output/timeseries/data/timeseries_R11_calibrated.rds")

# helper: hourly mean model per plot
prep_model <- function(input_model, plotnr){
  input_model %>%
    mutate(
      plot = plotnr,
      datetime_panama = with_tz(datetime, tzone = "EST"),
      datetime_hour   = floor_date(datetime_panama, unit = "hour")
    ) %>%
    group_by(plot, datetime_hour) %>%
    summarise(
      Mean_model = mean(temperature_model, na.rm = TRUE),
      .groups = "drop"
    )
}

model_all <- bind_rows(
  prep_model(input_model_C13, "C13"),
  prep_model(input_model_R11, "R11")
)

# helper: hourly mean observations per plot
temps_all <- input_temps %>%
  filter(sensor_name == "TMS_T3",
         plot %in% plots) %>%
  filter(datetime_EST >= t_start & datetime_EST <= t_end)

obs_all <- temps_all %>%
  mutate(datetime_hour = floor_date(datetime_EST, unit = "hour")) %>%
  group_by(plot, datetime_hour) %>%
  summarise(
    Mean_micro = mean(value, na.rm = TRUE),
    Macro      = mean(temp_macro, na.rm = TRUE),
    season     = names(which.max(table(season))),
    .groups = "drop"
  )

# merge obs + model per plot
df_all <- obs_all %>%
  left_join(model_all, by = c("plot", "datetime_hour"))

# wide format -> columns become curves
df_wide <- df_all %>%
  select(datetime_hour, plot, Macro, Mean_micro, Mean_model) %>%
  tidyr::pivot_wider(
    names_from  = plot,
    values_from = c(Macro, Mean_micro, Mean_model),
    names_sep   = "_"
  ) %>%
  arrange(datetime_hour)

# Macro is the same: keep only one (use C13 as reference)
df_wide <- df_wide %>%
  mutate(Macro = Macro_C13) %>%
  select(
    datetime_hour,
    Macro,
    Mean_micro_C13, Mean_model_C13,
    Mean_micro_R11, Mean_model_R11
  ) %>%
  rename(
    `micro C13`  = Mean_micro_C13,
    `model C13`  = Mean_model_C13,
    `micro R11`  = Mean_micro_R11,
    `model R11`  = Mean_model_R11
  )

# season shading: macro/season is identical -> just use season from C13 time series
season_rects <- df_all %>%
  filter(plot == "C13") %>%
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
x <- xts(df_wide[, c("Macro", "micro C13", "model C13", "micro R11", "model R11")],
         order.by = df_wide$datetime_hour)

names(x) <- c("macro", "micro TOMST-mean (C13)", "micro model-mean (C13)",
              "micro TOMST-mean (R11)", "micro model-mean (R11)")

p <- dygraph(
  x,
  main = paste(
    "Timeseries of observed macro- & TOMST-mean microtemperature and model-mean microtemperature for plots C13 and R11",
    "<br><span style='font-size:12px'>
       <span style='background-color:rgba(255,165,0,0.3);padding:2px 8px;'>dry season</span>
       &nbsp;
       <span style='background-color:rgba(100,149,237,0.35);padding:2px 8px;'>wet season</span>
     </span>"
  )
) %>%
  dyRangeSelector() %>%
  dyOptions(
    useDataTimezone = TRUE,
    drawPoints = FALSE,
    strokeWidth = 1.2
  ) %>%
  dyLegend(
    show = "always",
    hideOnMouseOut = FALSE,
    width = 280,
    labelsSeparateLines = TRUE
  ) %>%
  dyAxis("y", label = "Temperature (°C)") %>%

  # MACRO (dikke rode lijn)
  dySeries("macro",
           label = "macro",
           color = "#D7191C",
           strokeWidth = 3,
           drawPoints = FALSE) %>%

  # --- R11 ---
  dySeries("micro TOMST-mean (R11)",
           label = "micro TOMST-mean (R11)",
           color = "#006400",
           strokePattern = "dashed",
           strokeWidth = 2,
           drawPoints = FALSE) %>%

  dySeries("micro model-mean (R11)",
           label = "micro model-mean (R11)",
           color = "#66C266",
           strokeWidth = 2,
           drawPoints = FALSE) %>%

  # --- C13 ---
  dySeries("micro TOMST-mean (C13)",
           label = "micro TOMST-mean (C13)",
           color = "#08306B",
           strokePattern = "dashed",
           strokeWidth = 2,
           drawPoints = FALSE) %>%

  dySeries("micro model-mean (C13)",
           label = "micro model-mean (C13)",
           color = "#6BAED6",
           strokeWidth = 2,
           drawPoints = FALSE)


for(i in seq_len(nrow(season_rects))) {
  col <- if (season_rects$season[i] == "dry") "rgba(255,165,0,0.15)" else "rgba(100,149,237,0.25)"
  p <- p %>% dyShading(from = season_rects$from[i], to = season_rects$to[i], color = col)
}

saveWidget(
  p,
  file = "Output/timeseries/plots/timeseries_slider_C13_R11.html",
  selfcontained = TRUE
)
