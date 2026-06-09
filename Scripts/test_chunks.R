library(dplyr)
library(lubridate)
library(readxl)



time <- as.POSIXct("2023-09-28 15:00:00", tz = "UTC")

plotnr <- "C8"

PAR_data_full_2023 <<- readRDS("Data/PAR_data_full_2023.rds")

sensor_locations <- read_excel("Data/sensor_locations_all_plots.xlsx")
sensor_locations <- sensor_locations %>%
  mutate(
    plotnr_clean = gsub("plot0?", "", plot)
  ) %>%
  filter(plotnr_clean == plotnr)

# rename datetime column and create UTC datetime
PAR_obs <- PAR_data_full_2023 %>%
  rename(datetime_EST = datetime) %>%
  mutate(
    datetime_EST = force_tz(datetime_EST, tzone = "EST"),
    datetime_UTC = with_tz(datetime_EST, tzone = "UTC"),
    datetime_UTC = floor_date(datetime_UTC, unit = "hour"),
    location_ID_number = as.character(location_ID_number)
  ) %>%

  # keep only relevant columns
  select(
    radiation_mean,
    mean_ratio,
    plotNumber,
    datetime_UTC,
    location_ID_number
  ) %>%

  # keep only requested plotnr
  filter(plotNumber == substr(plotnr, 2, nchar(plotnr))) %>%

  # calculate hourly means per sensor
  group_by(location_ID_number, plotNumber, datetime_UTC) %>%
  summarise(
    radiation_mean = mean(radiation_mean, na.rm = TRUE),
    mean_ratio = mean(mean_ratio, na.rm = TRUE),
    .groups = "drop"
  ) %>%

  # keep only requested timestep
  filter(datetime_UTC == time) %>%

  # add sensor locations
  left_join(
    sensor_locations %>% select(SN, x, y),
    by = c("location_ID_number" = "SN")
  )

# save as csv
write.csv(
  PAR_obs,
  paste0("Data/PAR_obs_", format(time, "%Y%m%d_%H"), "_", plotnr, ".csv"),
  row.names = FALSE
)
