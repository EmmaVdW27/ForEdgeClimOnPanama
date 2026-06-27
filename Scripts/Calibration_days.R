start_analysis <- Sys.time()

library(dplyr)
library(lubridate)
library(readr)

# ---- PAR availability filter ----
par_data <- readRDS("Data/PAR_data_full_2023.rds") %>%
  filter(plot %in% c("8C", "7R")) %>%
  mutate(
    datetime_est = force_tz(datetime, tzone = "EST"),
    datetime_utc = with_tz(datetime_est, tzone = "UTC"),
    date = as.Date(datetime_utc),
    datetime_hour = floor_date(datetime_utc, "hour")
  )

valid_par_dates <- par_data %>%
  filter(!is.na(radiation_mean)) %>%
  group_by(date, plot, datetime_hour) %>%
  summarise(
    n_obs = n(),
    .groups = "drop"
  ) %>%
  filter(n_obs >= 10) %>%   # at least 10 observations = sensors per hour
  group_by(date, plot) %>%
  summarise(
    n_good_hours = n(),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(
    names_from = plot,
    values_from = n_good_hours,
    values_fill = 0
  ) %>%
  filter(`8C` >= 12, # at least 12h per day = daytime hours
         `7R` >= 12) %>%
  pull(date)

# ---- TEMPERATURE ----
temp_data <- readRDS("Data/TMS_macro_2023_C13_R11_R07_C08.rds") %>%
  mutate(datetime_hour = floor_date(datetime, "hour"))

temp_hourly <- temp_data %>%
  group_by(datetime_hour) %>%
  summarise(
    Tair = mean(temp_macro, na.rm = TRUE),
    .groups = "drop"
  )

# ---- RADIATION ----
hdr <- read_lines("Data/Table_CNR1_2023.dat", n_max = 2)[2]
col_names <- read_csv(I(hdr), col_names = FALSE, show_col_types = FALSE) |>
  unlist()

rad_data <- read_csv("Data/Table_CNR1_2023.dat", skip = 4,
                     col_names = col_names, show_col_types = FALSE) %>%
  mutate(
    TIMESTAMP = as.POSIXct(TIMESTAMP, tz = "EST"),
    datetime_hour = floor_date(with_tz(TIMESTAMP, "UTC"), "hour")
  )

rad_hourly <- rad_data %>%
  group_by(datetime_hour) %>%
  summarise(
    SW_total = mean(CM3Up_Avg, na.rm = TRUE),
    SW_diff  = mean(CM3Up_Avg * (PARdif_Avg / PARtot_Avg), na.rm = TRUE),
    SW_dir   = SW_total - SW_diff,
    .groups = "drop"
  )

# ---- RH ----
RH_hourly <- read_csv("Data/bci_cl_rh_elect.csv", show_col_types = FALSE) %>%
  mutate(
    datetime_est = dmy_hms(datetime, tz = "EST"),
    datetime_hour = floor_date(with_tz(datetime_est, "UTC"), "hour")
  ) %>%
  group_by(datetime_hour) %>%
  summarise(RH = mean(rh, na.rm = TRUE), .groups = "drop")

# ---- WIND ----
u_hourly <- read_csv("Data/bci_ava_ws_elect.csv", show_col_types = FALSE) %>%
  mutate(
    datetime_est = dmy_hms(datetime, tz = "EST"),
    datetime_hour = floor_date(with_tz(datetime_est, "UTC"), "hour")
  ) %>%
  group_by(datetime_hour) %>%
  summarise(u = mean(ws, na.rm = TRUE), .groups = "drop")


hourly_data <- temp_hourly %>%
  inner_join(rad_hourly, by = "datetime_hour") %>%
  inner_join(RH_hourly, by = "datetime_hour") %>%
  inner_join(u_hourly, by = "datetime_hour") %>%
  mutate(date = as.Date(datetime_hour))


df_daily <- hourly_data %>%
  group_by(date) %>%
  filter(sum(!is.na(Tair)) >= 20) %>%  # at least 20h of data
  summarise(
    sw_mean = mean(SW_total, na.rm = TRUE),
    sw_sd   = sd(SW_total, na.rm = TRUE),
    diffuse_frac = mean(SW_diff / SW_total, na.rm = TRUE),

    T_mean = mean(Tair, na.rm = TRUE),
    T_range = max(Tair, na.rm = TRUE) - min(Tair, na.rm = TRUE),

    RH_mean = mean(RH, na.rm = TRUE),
    u_mean  = mean(u, na.rm = TRUE),

    T_night = mean(Tair[hour(datetime_hour) %in% 0:5], na.rm = TRUE),
    T_day   = mean(Tair[hour(datetime_hour) %in% 12:15], na.rm = TRUE),
    T_daynight_diff = T_day - T_night,

    .groups = "drop"
  )

# make sure there is PAR data
df_daily <- df_daily %>%
  filter(date %in% valid_par_dates)

# ---- 1. Sunny and variable----
day_sunny_variable <- df_daily %>%
  slice_max(sw_mean * sw_sd, n = 1)

# ---- 2. Diffuse + humid ----
day_diffuse_humid <- df_daily %>%
  filter(date != day_sunny_variable$date) %>%
  mutate(score = diffuse_frac * RH_mean) %>%
  slice_max(score, n = 1)

# ---- 3. Windy ----
day_windy <- df_daily %>%
  filter(!date %in% c(day_sunny_variable$date,
                      day_diffuse_humid$date)) %>%
  slice_max(u_mean, n = 1)

# ---- 4. Day-night dynamics ----
day_night <- df_daily %>%
  filter(!date %in% c(day_sunny_variable$date,
                      day_diffuse_humid$date,
                      day_windy$date)) %>%
  slice_max(T_daynight_diff, n = 1)

# ---- Combine ----
selected_days <- c(
  day_sunny_variable$date,
  day_diffuse_humid$date,
  day_windy$date,
  day_night$date
)

print(selected_days)

end_analysis <- Sys.time()
cat(sprintf('\nTotal running time = %.2f min\n', as.numeric(end_analysis - start_analysis, units = 'mins')))
