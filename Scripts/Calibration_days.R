start_analysis <- Sys.time()

library(dplyr)
library(lubridate)
library(readr)

# ---- TEMPERATURE ----
temp_data <- readRDS("Data/TMS_macro_2023_C13_and_R11.rds") %>%
  mutate(datetime_hour = floor_date(datetime, "hour")) %>%
  filter(plot == plotnr)

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
  filter(sum(!is.na(Tair)) >= 20) %>%  # minstens 20 uren met data
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
