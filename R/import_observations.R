#' Function to import temperatures from observations
#'
#' @param time yyyy-mm-dd hh:mm:ss in UTC
#' @return macrotemperature, soil temperature, microtemperatures
#' @importFrom lubridate floor_date with_tz dmy_hms
#' @importFrom dplyr mutate group_by summarise select
#' @importFrom readr read_csv
#' @export
import_temp_observations <- function(time){

  # mutate data set
  temp_input_file <- temp_input_file_rds |>
    mutate(datetime_UTC_hour = floor_date(datetime, unit = "hour")) # datetime_UTC, floor the datetimes

  # only keep plot of interest
  temp_input_file <- temp_input_file[
    sub("^([A-Z])0+", "\\1", temp_input_file$plot) == plotnr,
  ]

  # only keep datetime of interest
  temp_input_file <- temp_input_file[temp_input_file$datetime_UTC_hour == time, ]

  # MACROTEMPERATURE
  macro_temp <<- mean(temp_input_file$temp_macro, na.rm = TRUE) + 273.15 # in Kelvin

  # SOILTEMPERATURE
  # only keep soil sensors
  soil_temp_input_file <- temp_input_file[temp_input_file$sensor_name == "TMS_T1", ]
  T_soil_deep <<- mean(soil_temp_input_file$value, na.rm = TRUE) + 273.15 # in Kelvin

  # MICROTEMPERATURE DATAFRAME
  # only keep air sensors
  air_temp_input_file <- temp_input_file[temp_input_file$sensor_name == "TMS_T3", ]
  # only keep columns of interest
  air_temp_obs <- air_temp_input_file[, c("SN", "value")]
  # get mean hourly value per sensor
  air_temp_obs <- air_temp_obs |>
    group_by(SN) |>
    summarise(
      value_mean = mean(value, na.rm = TRUE),
      .groups = "drop"
    )
  # add sensor locations from sensor_locations file
  air_temp_obs <- air_temp_obs|>
    left_join(sensor_locations |> select(SN, x, y), by = "SN")

  # save dataframe as CSV
  write.csv(air_temp_obs, paste0("Data/microtemp_air_obs_", format(time, "%Y%m%d_%H"), "_", plotnr, ".csv"), row.names = FALSE)

}

#' Function to import radiative observations (direct and diffuse shortwave and longwave radiation and measured PAR)
#'
#' @param time yyyy-mm-dd hh:mm:ss in UTC
#' @return Direct-beam and diffuse shortwave and longwave radiation and PAR measurements
#' @importFrom lubridate floor_date with_tz force_tz
#' @importFrom dplyr mutate rename select filter group_by left_join
#' @importFrom readr read_csv read_lines cols col_character
#' @export
import_rad_observations <- function(time){

  # TOP CANOPY RADIATION
  # read column names
  hdr <- read_lines(rad_input_file_csv, n_max = 2)[2]
  col_names <- read_csv(I(hdr), col_names = FALSE, show_col_types = FALSE) |>
    unlist(use.names = FALSE)

  rad_input_file <<- read_csv(rad_input_file_csv, skip = 4,
                              col_names = col_names, col_types = cols(TIMESTAMP = col_character()),
                              show_col_types = FALSE)

  # mutate data set
  rad_input_file <- rad_input_file |>
    mutate(TIMESTAMP = as.POSIXct(TIMESTAMP, tz = "EST")) |>
    mutate(TIMESTAMP_UTC = with_tz(TIMESTAMP, "UTC")) |> # time zone conversion
    mutate(TIMESTAMP_UTC_hour = floor_date(TIMESTAMP_UTC, unit = "hour")) # floor the TIMESTAMPs

  # only keep datetime of interest
  rad_input_file <- rad_input_file[rad_input_file$TIMESTAMP_UTC_hour == time, ]

  # mutate data set: calculate mean values
  # REMARK: Up and Dn (Down) refer to the direction in which the sensor is looking, so Up measure the downward radiation
  rad_input_file_mean <- rad_input_file |>
    summarise(
      CM3Up_Avg    = mean(CM3Up_Avg,    na.rm = TRUE),
      CG3UpCo_Avg  = mean(CG3UpCo_Avg,  na.rm = TRUE),
      CG3DnCo_Avg  = mean(CG3DnCo_Avg,  na.rm = TRUE),
      PARtot_Avg   = mean(PARtot_Avg,   na.rm = TRUE),
      PARdif_Avg   = mean(PARdif_Avg,   na.rm = TRUE)
    )

  # Diffuse shortwave fraction
  ratio = rad_input_file_mean$PARdif_Avg / rad_input_file_mean$PARtot_Avg

  # DIRECT-BEAM SHORTWAVE RADIATION
  F_sky_dir_init <<- rad_input_file_mean$CM3Up_Avg * (1 - ratio) # in W/m²

  # DIFFUSE SHORTWAVE RADIATION
  F_sky_diff_init <<- rad_input_file_mean$CM3Up_Avg * ratio # in W/m²

  # LONGWAVE RADIATION INCOMING
  F_sky_lw <<- rad_input_file_mean$CG3UpCo_Avg  # in W/m²

  # LONGWAVE RADIATION OUTGOING
  F_sky_lw_out <<- rad_input_file_mean$CG3DnCo_Avg  # in W/m²

  # BELOW CANOPY RADIATION
  sensor_locations <- sensor_locations |>
    mutate(
      plotnr_clean = gsub("plot0?", "", plot)
    ) |>
    filter(plotnr_clean == plotnr)

  # rename datetime column and create UTC datetime
  PAR_obs <- PAR_input_file_rds |>
    rename(datetime_EST = datetime) |>
    mutate(
      datetime_EST = force_tz(datetime_EST, tzone = "EST"),
      datetime_UTC = with_tz(datetime_EST, tzone = "UTC"),
      datetime_UTC = floor_date(datetime_UTC, unit = "hour"),
      location_ID_number = as.character(location_ID_number)
    ) |>

    # keep only relevant columns
    select(
      radiation_mean,
      mean_ratio,
      plotNumber,
      datetime_UTC,
      location_ID_number
    ) |>

    # keep only requested plotnr
    filter(plotNumber == substr(plotnr, 2, nchar(plotnr))) |>

    # calculate hourly means per sensor
    group_by(location_ID_number, plotNumber, datetime_UTC) |>
    summarise(
      radiation_mean = mean(radiation_mean, na.rm = TRUE),
      mean_ratio = mean(mean_ratio, na.rm = TRUE),
      .groups = "drop"
    ) |>

    # keep only requested timestep
    filter(datetime_UTC == time) |>

    # add sensor locations
    left_join(
      sensor_locations |> select(SN, x, y),
      by = c("location_ID_number" = "SN")
    )

  # save as csv
  write.csv(
    PAR_obs,
    paste0("Data/PAR_obs_", format(time, "%Y%m%d_%H"), "_", plotnr, ".csv"),
    row.names = FALSE
  )

}

#' Function to import relative humidity and wind speed observations from BCI stations: resp. forest clearing and AVA tower
#'
#' @param time yyyy-mm-dd hh:mm:ss in UTC
#' @return RH relative humidity
#' @return u_top wind speed above canopy
#' @importFrom lubridate floor_date with_tz dmy_hms
#' @importFrom dplyr mutate
#' @importFrom readr read_csv
#' @export
import_RH_u_observations <- function(time){

  # read files
  RH_input <- read_csv(RH_input_file_csv, show_col_types = FALSE)
  u_input <- read_csv(u_input_file_csv, show_col_types = FALSE)

  # mutate RH file
  RH_input$datetime_est <- dmy_hms(RH_input$datetime, tz = "EST") # datetime is character string with time in EST = UTC-5 (no summer or winter time)
  RH_input$datetime_utc <- with_tz(RH_input$datetime_est, tzone = "UTC") # convert to UTC
  RH_hourly <- RH_input |> # hourly averages
    mutate(datetime_hour = floor_date(datetime_utc, "hour")) |>
    group_by(datetime_hour) |>
    summarise(rh_hour = mean(rh, na.rm = TRUE),
              .groups = "drop")

  # mutate u file
  u_input$datetime_est <- dmy_hms(u_input$datetime, tz = "EST") # datetime is character string with time in EST = UTC-5 (no summer or winter time)
  u_input$datetime_utc <- with_tz(u_input$datetime_est, tzone = "UTC") # convert to UTC
  u_hourly <- u_input |> # hourly averages
    mutate(datetime_hour = floor_date(datetime_utc, "hour")) |>
    group_by(datetime_hour) |>
    summarise(u_hour = mean(ws, na.rm = TRUE),
              .groups = "drop")

  # return values
  RH <<- RH_hourly$rh_hour[RH_hourly$datetime_hour == time]
  u_top <<- u_hourly$u_hour[u_hourly$datetime_hour == time]

}
