# --------------------------------------------------
# Investigation of LE data, given by Félicien
# --------------------------------------------------

data <- read.csv("Data/Data_Félicien/BCI_v5.1.csv", stringsAsFactors = FALSE)

data$date <- as.POSIXct(data$date, format = "%Y-%m-%d %H:%M")
data$hour <- as.numeric(format(data$date, "%H"))
night_data <- subset(data, hour >= 18 | hour < 6)

library(ggplot2)

p = ggplot(night_data, aes(x = date, y = Rnet)) +
  geom_line() +
  theme_minimal() +
  labs(title = "LE over time",
       x = "Date",
       y = "LE")

print(p)


# --------------------------------------------------
# Dynamic macro boundary, proxy nocturnal atmosphere-canopy decoupling
# --------------------------------------------------

# keep original macro_temp
if (iter == 1) {
  macro_temp_orig <- macro_temp
}

# Reset macro_temp every iteration towards original
macro_temp <- macro_temp_orig

# get top (5) layers per column
top5_per_col <- micro_grid |>
  dplyr::filter(den > 0) |>
  dplyr::group_by(x, y) |>
  dplyr::arrange(desc(z), .by_group = TRUE) |>
  dplyr::slice_head(n = 5) |>
  dplyr::ungroup()

mean_top_temp <- mean(top5_per_col$temperature, na.rm = TRUE)

# proxy for macro-canopy decoupling
if (!is.na(mean_top_temp) && mean_top_temp < macro_temp) {
  macro_temp <- mean_top_temp
}


# --------------------------------------
# PROXY WET SURFACE EVAPORATION
# --------------------------------------

P_sat <- 0.3
RH_threshold <- 90
k_wet <- 100

RH <- 95.9 # 29/08/2023 13h UTC
rain <- 0

# convert air temperature to Celsius
T_air_C <- T_air_vec - 273.15

# saturation vapour pressure (kPa)
es <- 0.611 * exp((17.27 * T_air_C)/(T_air_C + 237.3))

# VPD
VPD <- es * (1 - RH/100)
VPD <- pmax(VPD, 0)
VPD <- pmin(VPD, 0.5)  # stabilisation

# canopy wetness
f_wet <- pmin(1,
              rain / P_sat +
                pmax(0, (RH - RH_threshold)/(100 - RH_threshold)))

# wet evaporation
LE_wet <- pmin(k_wet * f_wet * VPD, 20)

latent_flux <- latent_flux + LE_wet


# -----------------------------
# MACROTEMPERATURE from clearing or Lutz
# -----------------------------

macro_temp_input <- read_csv(macro_input_file_csv, show_col_types = FALSE)
macro_temp_input$datetime_est <- dmy_hms(macro_temp_input$datetime, tz = "EST") # datetime is character string with time in EST = UTC-5 (no summer or winter time)

macro_temp_input$datetime_utc <- with_tz(macro_temp_input$datetime_est, tzone = "UTC") # convert to UTC

macro_hourly <- macro_temp_input |> # hourly averages
  mutate(datetime_hour = floor_date(datetime_utc, "hour")) |>
  group_by(datetime_hour) |>
  summarise(at_hour = mean(at, na.rm = TRUE),
            .groups = "drop")

macro_temp <<- macro_hourly$at_hour[macro_hourly$datetime_hour == time] + 273.15 # specific hour and °C to K

# -----------------------------
# MACROTEMPERATURE from AVA tower 48 m
# -----------------------------

T_48m <- mean(temp_input_file$temp_macro, na.rm = TRUE)  # in °C

hour_now <- as.numeric(format(time, "%H"))

# Define sunset and sunrise (local time assumption -> -5)
sunset  <- 23
sunrise <- 11

# Maximum nocturnal difference (°C)
deltaT_max <- 2
k <- 0.7  # stabilisation rate

if(hour_now >= sunset | hour_now < sunrise){

  # Compute hours since sunset
  if(hour_now >= sunset){
    t_since_sunset <- hour_now - sunset
  } else {
    t_since_sunset <- (24 - sunset) + hour_now
  }

  # Exponential growth of nocturnal stratification
  deltaT <- deltaT_max * (1 - exp(-k * t_since_sunset))

} else {

  deltaT <- 0

}

# Canopy-top macro temperature (35 m)
macro_temp <<- (T_48m - deltaT) + 273.15  # convert to Kelvin

