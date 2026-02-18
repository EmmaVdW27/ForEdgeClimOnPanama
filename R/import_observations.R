#' Function to import temperatures from observations
#'
#' @param datetime yyyy-mm-dd hh:mm:ss in UTC
#' @return macrotemperature, soil temperature, microtemperatures
#' @importFrom lubridate floor_date with_tz
#' @importFrom dplyr mutate group_by summarise select
#' @export
import_temp_observations <- function(time){

  # mutate data set
  temp_input_file <- temp_input_file_rds |>
    #mutate(datetime_UTC = with_tz(datetime_EST, "UTC")) |> # time zone conversion (datetime_EST is already a POSIXct object)
    mutate(datetime_UTC_hour = floor_date(datetime, unit = "hour")) # datetime_UTC, floor the datetimes

  # only keep plot of interest
  temp_input_file <- temp_input_file[temp_input_file$plot == plotnr, ]

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

#' Function to import radiative observations (direct and diffuse shortwave and longwave radiation)
#'
#' @param datetime yyyy-mm-dd hh:mm:ss in UTC
#' @return Direct-beam and diffuse shortwave and longwave radiation
#' @importFrom lubridate floor_date with_tz
#' @importFrom dplyr mutate
#' @importFrom readr read_csv read_lines cols col_character
#' @export
import_rad_observations <- function(datetime){

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
  rad_input_file <- rad_input_file[rad_input_file$TIMESTAMP_UTC_hour == datetime, ]

  # mutate data set: calculate mean values
  # REMARK: Up and Dn (Down) refer to the direction in which the sensor is looking, so Up measure the downward radiation
  rad_input_file_mean <- rad_input_file |>
    summarise(
      CM3Up_Avg    = mean(CM3Up_Avg,    na.rm = TRUE),
      CG3UpCo_Avg  = mean(CG3UpCo_Avg,  na.rm = TRUE),
      PARtot_Avg   = mean(PARtot_Avg,   na.rm = TRUE),
      PARdif_Avg   = mean(PARdif_Avg,   na.rm = TRUE)
    )

  # Diffuse shortwave fraction
  ratio = rad_input_file_mean$PARdif_Avg / rad_input_file_mean$PARtot_Avg

  # DIRECT-BEAM SHORTWAVE RADIATION
  F_sky_dir_init <<- rad_input_file_mean$CM3Up_Avg * (1 - ratio) # in W/m²

  # DIFFUSE SHORTWAVE RADIATION
  F_sky_diff_init <<- rad_input_file_mean$CM3Up_Avg * ratio # in W/m²

  # DIFFUSE SHORTWAVE RADIATION
  F_sky_lw <<- rad_input_file_mean$CG3UpCo_Avg # in W/m²

}


