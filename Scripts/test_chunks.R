library(data.table)

############################################
# 1. Dummy lw_two_stream functie
############################################

lw_two_stream_dummy <- function(soil_reflect, density, F_sky_lw,
                                temperature, T_soil, T_macro,
                                Kd, omega, beta) {

  n <- length(density)

  list(
    F_d_down = density * 1,
    F_d_up   = density * 2,
    net_lw   = density * 3
  )
}

############################################
# 2. Functie 1 (data.table versie)
############################################

longwave_two_stream_RTM_dt <- function(voxel_grid, micro_grid, lw_two_stream,
                                       F_sky_lw, omega_g_lw_v, macro_temp,
                                       Kd_lw_v, omega_lw, beta_lw) {

  voxel_grid <- as.data.table(voxel_grid)
  micro_grid <- as.data.table(micro_grid)

  setnames(micro_grid, old = c("x", "y", "z"), new = c("X", "Y", "Z"))

  micro_grid <- micro_grid[, .(X, Y, Z, temperature,
                               T_soil = ifelse(Z == 1, T_soil, NA))]

  voxel_grid <- voxel_grid[micro_grid, on = .(X, Y, Z)]

  final_results_lw_vertical <- voxel_grid[, {

    T_soil_value <- na.omit(T_soil)[1]
    n_layers <- .N

    fluxes <- lw_two_stream(
      soil_reflect = omega_g_lw_v,
      density = density,
      F_sky_lw = F_sky_lw,
      temperature = temperature,
      T_soil = T_soil_value,
      T_macro = macro_temp,
      Kd = Kd_lw_v,
      omega = omega_lw,
      beta = beta_lw
    )

    .(Z = seq_len(n_layers),
      density, temperature,
      F_d_down = fluxes$F_d_down[1:n_layers],
      F_d_up   = fluxes$F_d_up[1:n_layers],
      net_lw   = fluxes$net_lw[1:n_layers])

  }, by = .(X, Y)]

  setorder(final_results_lw_vertical, X, Y, Z)

  return(final_results_lw_vertical)
}

shortwave_two_stream_RTM_dt <- function(datetime, lat, lon, voxel_grid,
                                     F_sky_diff_init, F_sky_dir_init,
                                     omega_g_v,
                                     Kd_v, Kb_v,
                                     omega, betad, beta0) {

  voxel_grid <- as.data.table(voxel_grid)

  # Calculate solar angles
  calculate_sun_angle <- function(datetime, lat, lon) {
    sun_angle <- getSunlightPosition(date = datetime, lat = lat, lon = lon)
    sun_altitude <- sun_angle$altitude
    sun_azimuth <- sun_angle$azimuth
    in_eastern_half <- ifelse(sun_azimuth < 0, 1, 0)
    # Is sun above the horizon?
    # This is 1 if altitude is between 0 and pi/2 radians
    above_horizon <- ifelse(sun_altitude >= 0, 1, 0)
    return(list(above_horizon = above_horizon,
                sun_altitude = sun_altitude, sun_azimuth = sun_azimuth))
  }

  sun_angles <- calculate_sun_angle(datetime, lat, lon)

  ##########################
  # Vertical shortwave RTM #
  ##########################

  # Diffuse radiation from canopy top
  F_sky_diff_v <- F_sky_diff_init
  # Direct radiation from canopy top
  F_sky_dir_v <<- F_sky_dir_init * sun_angles$above_horizon

  # Calculate fluxes per unique XY combination
  final_results_vertical <- voxel_grid[, {
    # Number of layers
    n_layers <- .N

    fluxes <- sw_two_stream(
      soil_reflect = omega_g_v, density = density,
      F_sky_diff = F_sky_diff_v, F_sky_dir = F_sky_dir_v,
      Kd = Kd_v, Kb = Kb_v,
      omega = omega, beta = betad, beta0 = beta0
    )

    .(X, Y, Z = seq_len(n_layers), density,
      F_d_down = fluxes$F_d_down[1:n_layers],
      F_d_up = fluxes$F_d_up[1:n_layers],
      F_b_down = fluxes$F_b_down[1:n_layers],
      net_sw = fluxes$net_sw[1:n_layers])
  }, by = .(X, Y)]


  setorder(final_results_vertical, X, Y, Z)

  # Select relevant columns
  final_results <- final_results_vertical[, .(X, Y, Z, density, F_d_down, F_d_up, F_b_down, net_sw)]

  return(final_results)
}

############################################
# 3. Functie 2 (snelle loop versie)
############################################

longwave_two_stream_RTM_fast <- function(voxel_grid, micro_grid, lw_two_stream,
                                         F_sky_lw, omega_g_lw_v, macro_temp,
                                         Kd_lw_v, omega_lw, beta_lw) {

  X <- voxel_grid$X
  Y <- voxel_grid$Y
  Z <- voxel_grid$Z
  density <- voxel_grid$density

  temperature <- micro_grid$temperature
  T_soil <- micro_grid$T_soil

  n <- length(X)

  F_d_down <- numeric(n)
  F_d_up   <- numeric(n)
  net_lw   <- numeric(n)

  # precompute kolommen (GEEN data copies)
  xy_id <- X + Y * max(X)
  split_idx <- split(seq_len(n), xy_id)

  for (idx in split_idx) {

    density_col <- density[idx]
    temp_col    <- temperature[idx]

    T_soil_value <- na.omit(T_soil[idx])[1]

    fluxes <- lw_two_stream(
      soil_reflect = omega_g_lw_v,
      density = density_col,
      F_sky_lw = F_sky_lw,
      temperature = temp_col,
      T_soil = T_soil_value,
      T_macro = macro_temp,
      Kd = Kd_lw_v,
      omega = omega_lw,
      beta = beta_lw
    )

    n_layers <- length(idx)

    F_d_down[idx] <- fluxes$F_d_down[1:n_layers]
    F_d_up[idx]   <- fluxes$F_d_up[1:n_layers]
    net_lw[idx]   <- fluxes$net_lw[1:n_layers]
  }

  out <- data.frame(
    X = X,
    Y = Y,
    Z = Z,
    density = density,
    temperature = temperature,
    F_d_down = F_d_down,
    F_d_up   = F_d_up,
    net_lw   = net_lw
  )

  out <- out[order(out$X, out$Y, out$Z), ]

  return(out)
}

shortwave_two_stream_RTM_fast <- function(datetime, lat, lon, voxel_grid,
                                     F_sky_diff_init, F_sky_dir_init,
                                     omega_g_v,
                                     Kd_v, Kb_v,
                                     omega, betad, beta0) {

  # extract vectors
  X <- voxel_grid$X
  Y <- voxel_grid$Y
  Z <- voxel_grid$Z
  density <- voxel_grid$density

  n <- length(X)

  # solar angles
  sun_angle <- getSunlightPosition(date = datetime, lat = lat, lon = lon)
  sun_altitude <- sun_angle$altitude
  above_horizon <- ifelse(sun_altitude >= 0, 1, 0) # This is 1 if altitude is between 0 and pi/2 radians

  F_sky_diff_v <- F_sky_diff_init
  F_sky_dir_v  <- F_sky_dir_init * above_horizon

  # pre-allocation
  F_d_down <- numeric(n)
  F_d_up   <- numeric(n)
  F_b_down <- numeric(n)
  net_sw   <- numeric(n)

  # precompute columns
  xy_id <- interaction(X, Y, drop = TRUE)
  split_idx <- split(seq_len(n), xy_id)

  for (idx in split_idx) {

    density_col <- density[idx]

    ##########################
    # Vertical shortwave RTM #
    ##########################

    fluxes <- sw_two_stream(
      soil_reflect = omega_g_v,
      density = density_col,
      F_sky_diff = F_sky_diff_v,
      F_sky_dir = F_sky_dir_v,
      Kd = Kd_v,
      Kb = Kb_v,
      omega = omega,
      beta = betad,
      beta0 = beta0
    )

    n_layers <- length(idx)

    F_d_down[idx] <- fluxes$F_d_down[1:n_layers]
    F_d_up[idx]   <- fluxes$F_d_up[1:n_layers]
    F_b_down[idx] <- fluxes$F_b_down[1:n_layers]
    net_sw[idx]   <- fluxes$net_sw[1:n_layers]
  }

  return(data.frame(
    X = X,
    Y = Y,
    Z = Z,
    density = density,
    F_d_down = F_d_down,
    F_d_up = F_d_up,
    F_b_down = F_b_down,
    net_sw = net_sw
  ))

}


############################################
# 4. Testdata (toy example)
############################################

voxel_grid <- data.frame(
  X = c(1,1,1, 2,2,2),
  Y = c(1,1,1, 1,1,1),
  Z = c(1,2,3, 1,2,3),
  density = c(0, 0.5, 1, 0, 0.5, 1)
)

micro_grid <- data.frame(
  x = voxel_grid$X,
  y = voxel_grid$Y,
  z = voxel_grid$Z,
  temperature = c(10, 15, 20, 10, 15, 20),
  T_soil = c(12, NA, NA, 12, NA, NA)
)

############################################
# 5. Run beide functies
############################################

out1_lw <- longwave_two_stream_RTM_dt(
  voxel_grid, micro_grid, lw_two_stream_dummy,
  F_sky_lw = 1, omega_g_lw_v = 0.2, macro_temp = 15,
  Kd_lw_v = 0.5, omega_lw = 0.3, beta_lw = 0.1
)

out2_lw <- longwave_two_stream_RTM_fast(
  voxel_grid, micro_grid, lw_two_stream_dummy,
  F_sky_lw = 1, omega_g_lw_v = 0.2, macro_temp = 15,
  Kd_lw_v = 0.5, omega_lw = 0.3, beta_lw = 0.1
)

out1_sw <- shortwave_two_stream_RTM_dt(
  Sys.time(), 50, 3,
  voxel_grid,
  1, 2, 0.2,
  0.5, 0.6,
  0.3, 0.1, 0.05
)

out2_sw <- shortwave_two_stream_RTM_fast(
  Sys.time(), 50, 3,
  voxel_grid,
  1, 2, 0.2,
  0.5, 0.6,
  0.3, 0.1, 0.05
)


############################################
# 6. Vergelijking
############################################

out1 = out1_lw
out2 = out2_lw

cat("STRUCTUUR:\n")
str(out1)
str(out2)

cat("\nALL.EQUAL:\n")
print(all.equal(
  as.data.frame(out1[order(out1$X, out1$Y, out1$Z), ]),
  out2[order(out2$X, out2$Y, out2$Z), ],
  tolerance = 1e-12
))

cat("\nMAX verschil net rad:\n")
print(max(abs(out1$net_lw - out2$net_lw)))

############################################
# 7. Extra: shuffle test (detecteert bugs)
############################################

set.seed(42)

voxel_shuffled <- voxel_grid[sample(nrow(voxel_grid)), ]
micro_shuffled <- micro_grid[sample(nrow(micro_grid)), ]

out1_s <- longwave_two_stream_RTM_dt(
  voxel_shuffled, micro_shuffled, lw_two_stream_dummy,
  1, 0.2, 15, 0.5, 0.3, 0.1
)

out2_s <- longwave_two_stream_RTM_fast(
  voxel_shuffled, micro_shuffled, lw_two_stream_dummy,
  1, 0.2, 15, 0.5, 0.3, 0.1
)

cat("\nSHUFFLE TEST:\n")
print(all.equal(
  out1_s[order(out1_s$X, out1_s$Y, out1_s$Z), ],
  out2_s[order(out2_s$X, out2_s$Y, out2_s$Z), ]
))
