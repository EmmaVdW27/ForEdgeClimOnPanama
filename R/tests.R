library(lubridate)
library(dplyr)

# input files
temp_input_file_rds <<- readRDS("Data/Data_Kasper/TMS_Emma.rds")
PAR_input_file_rds <<- readRDS("Data/Data_Kasper/PAR_data_Emma.rds")

datetime = as.POSIXct("2023-04-26 14:00:00", tz = "UTC")
plotnr <<- "C13" # or R11
sigma_SB <<- 5.67e-8

# mutate data set
PAR_input_file <- PAR_input_file_rds |>
  mutate(datetime = as.POSIXct(datetime, tz = "EST")) |>
  mutate(datetime_UTC = with_tz(datetime, "UTC")) |> # time zone conversion
  mutate(datetime_UTC_hour = floor_date(datetime_UTC, unit = "hour")) # floor the datetimes

# only keep one logger location ID
PAR_input_file <- PAR_input_file[PAR_input_file$location_ID == 1, ]

# only keep datetime of interest
PAR_input_file <- PAR_input_file[PAR_input_file$datetime_UTC_hour == datetime, ]

# mutate data set: create 2 new columns with direct and diffuse light
PAR_input_file <- PAR_input_file |>
  mutate(diffuse = ref_mean * mean_ratio) |>
  mutate(direct = ref_mean - diffuse)

# DIRECT SHORTWAVE RADIATION
# ref_mean in µmol/m²/s
F_sky_dir_init <<- mean(PAR_input_file$direct, na.rm = TRUE) * 0.45 # in W/m²; conversion factor is 0.45

# DIFFUSE SHORTWAVE RADIATION
F_sky_diff_init <<- mean(PAR_input_file$diffuse, na.rm = TRUE) * 0.45 # in W/m²; conversion factor is 0.45
