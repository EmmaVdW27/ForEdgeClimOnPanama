#' Two-stream RTM for shortwave radiation and for a single column existing of layers/voxels with density values (%).
#' This model is inspired by the shortwave RTM of the Ecosystem Demography model (ED 2.2).
#' @param soil_reflect Soil albedo (%)
#' @param density Density of each layer (%)
#' @param F_sky_diff Diffuse radiation (W/m2)
#' @param F_sky_dir Direct beam radiation (W/m2)
#' @param Kd Diffuse extinction coefficient, per unit density
#' @param Kb Direct beam extinction coefficient, per unit density
#' @param omega Shortwave scattering coefficient
#' @param beta Fraction of scattering in the backward direction for diffuse radiation (%)
#' @param beta0 Fraction of scattering in the backward direction for direct beam radiation (%)
#' @return Dataframe containing the combined results of the 2D shortwave RTM
#' @export

sw_two_stream <- function(soil_reflect, density, F_sky_diff, F_sky_dir, Kd, Kb, omega, beta, beta0) {

  omega_g = soil_reflect
  # Number of (cohort) layers
  n <- length(density)
  # Add atmosphere layer
  density <- c(density, 0)
  # Inverse optical depths
  # mu atmosphere (density = 0, no forest structure) is defined to be 1
  mu0 <- ifelse(density == 0, 1, 1 / (Kb * density)) # Inverse optical depth direct beam radiation
  mu <- ifelse(density == 0, 1, 1 / (Kd * density))  # Inverse optical depth diffuse radiation

  # Initialization matrices
  S <- matrix(0, nrow = 2 * n + 2, ncol = 2 * n + 2)
  y <- numeric(2 * n + 2)

  # Calculate direct radiation profile (exp decay from top to bottom)
  direct_radiation <- F_sky_dir * rev(cumprod(rev(exp(-density / mu0))))

  # Auxiliary variables
  # Gamma
  gamma_plus <- rep(0.5 * (1 + sqrt((1 - omega) / (1 - (1 - 2 * beta) * omega))) , n)
  gamma_plus[n+1] = 1
  gamma_minus <- rep(0.5 * (1 - sqrt((1 - omega) / (1 - (1 - 2 * beta) * omega))) , n)
  gamma_minus[n+1] = 0
  # Lambda
  lambda <- sqrt((1 - (1 - 2 * beta) * omega) * (1 - omega) / (mu)^2)
  lambda[n+1] = 0
  # Kappa
  kappa_plus <- -((1 - (1 - 2 * beta) * omega)/ mu[1:n] + (1 - 2 * beta0) / mu0[1:n]) * omega * direct_radiation[2:(n+1)] / mu0[1:n]
  kappa_plus[n+1] =  -F_sky_dir
  kappa_minus = - ( (1-2*beta0)*(1-omega)/mu[1:n] + 1/mu0[1:n]) * omega*direct_radiation[2:(n+1)]/mu0[1:n]
  kappa_minus[n+1] = -F_sky_dir
  # Delta
  delta_plus <- (kappa_plus + kappa_minus) * (mu0)^2 / (2 * (1 - lambda^2 * (mu0)^2))
  delta_plus[n+1] = 0.5*(kappa_plus[n+1]+kappa_minus[n+1])
  delta_minus <- (kappa_plus - kappa_minus) * (mu0)^2 / (2 * (1 - lambda^2 * (mu0)^2))
  delta_minus[n+1] = 0.5*(kappa_plus[n+1]-kappa_minus[n+1])

  # Boundary conditions
  y[2*n+2] <- F_sky_diff - delta_plus[n+1] # Incoming diffuse radiation at top canopy
  y[1] <- omega_g*direct_radiation[1] - (delta_minus[1] - omega_g*delta_plus[1])*exp(-density[1]/mu0[1]) # Outgoing diffuse radiation from ground surface

  S[2 * n + 2, 2 * n + 1] <- gamma_plus[n+1]
  S[2 * n + 2, 2 * n + 2] <- gamma_minus[n+1]
  S[1, 1] <- (gamma_minus[1] - omega_g*gamma_plus[1])*exp(-lambda[1]*density[1])
  S[1, 2] <- (gamma_plus[1] - omega_g*gamma_minus[1])*exp(lambda[1]*density[1])

  # **Vectorized matrix filling**
  indices <- 1:n
  idx_even <- 2 * indices
  idx_odd <- idx_even - 1
  idx_next <- indices + 1

  # Calculate exponential terms
  exp_lambda_density_neg <- exp(-lambda[idx_next] * density[idx_next])
  exp_lambda_density_pos <- exp(lambda[idx_next] * density[idx_next])

  exp_density_mu0_neg <- exp(-density[idx_next] / mu0[idx_next])

  # Populate the matrix S
  S[cbind(idx_even, idx_odd)] <- gamma_plus[indices]
  S[cbind(idx_even, idx_even)] <- gamma_minus[indices]
  S[cbind(idx_even, idx_even + 1)] <- -gamma_plus[idx_next] * exp_lambda_density_neg
  S[cbind(idx_even, idx_even + 2)] <- -gamma_minus[idx_next] * exp_lambda_density_pos

  S[cbind(idx_even + 1, idx_odd)] <- gamma_minus[indices]
  S[cbind(idx_even + 1, idx_even)] <- gamma_plus[indices]
  S[cbind(idx_even + 1, idx_even + 1)] <- -gamma_minus[idx_next] * exp_lambda_density_neg
  S[cbind(idx_even + 1, idx_even + 2)] <- -gamma_plus[idx_next] * exp_lambda_density_pos

  # Populate the vector y
  y[idx_even] <- delta_plus[idx_next] * exp_density_mu0_neg - delta_plus[indices]
  y[idx_even + 1] <- delta_minus[idx_next] * exp_density_mu0_neg - delta_minus[indices]

  # Solving matrix equation using arma::solve from C++
  x <- solve_rtm(S, y)

  # Calculate diffuse radiative fluxes
  # Radiative fluxes refer to the rate of energy transfer by electromagnetic radiation through a surface or area,
  # typically measured in W/m², encompassing all or specific wavelengths (e.g., shortwave or longwave radiation).
  # These values are mainly of interest when modelling e.g. energy balance, heat fluxes...
  # Every down[k] and up[k] respectively refer to the downward and upward flux below layer k.
  indices <- seq_len(n + 1)
  k2 <- 2 * indices
  k2m1 <- k2 - 1

  down <- x[k2m1] * gamma_plus[indices] * exp(-lambda[indices]*density[indices]) +
      x[k2] * gamma_minus[indices] * exp(lambda[indices]*density[indices]) +
      delta_plus[indices] * exp(-density[indices]/mu0[indices])

  up <- x[k2m1] * gamma_minus[indices] * exp(-lambda[indices]*density[indices]) +
      x[k2] * gamma_plus[indices] * exp(lambda[indices]*density[indices]) +
      delta_minus[indices] * exp(-density[indices]/mu0[indices])

  # Calculate light levels and netto absorbed or emitted radiation
  # Light levels refer to the intensity of photosynthetically active radiation (PAR, 400-700 nm) available for plants.
  # These values are mainly of interest when modelling e.g. photosynthesis, plant growth, shadow structure...
  indices2 <- seq_len(n)
  kp1 <- indices2 + 1

  light_level <- 0.5 * (down[indices2] + down[kp1] + direct_radiation[indices2] + direct_radiation[kp1])
  light_beam_level <- 0.5 * (direct_radiation[indices2] + direct_radiation[kp1])
  light_diff_level <- 0.5 * (down[indices2] + down[kp1])
  net_sw <- (down[kp1] - down[indices2]) + (up[indices2] - up[kp1]) + (direct_radiation[kp1] - direct_radiation[indices2] )

  # There are no light levels for atmosphere
  light_level = c(light_level, NA)
  light_beam_level = c(light_beam_level, NA)
  light_diff_level = c(light_diff_level, NA)

  return(list(
    layer = c(paste("layer", 1:n), "atmosphere"),
    density = density,
    F_d_down = round(down, 0),
    F_d_up = round(up, 0),
    F_b_down = round(direct_radiation, 0),
    net_sw = c(round(net_sw, 0), 0),  # Net radiation in atmosphere is set to zero (no structure here)
    albedo = c(rep(' ', n), round( up[length(up)]/(F_sky_diff + F_sky_dir), 3 ) )
  ))


}


#' 1D (vertical) shortwave two-stream RTM
#'
#' This function performs a shortwave
#' two-stream RTM to estimate radiation fluxes in a forest voxel grid.
#'
#' @param datetime Datetime object representing the current simulation time
#' @param lat Latitude of the location
#' @param lon Longitude of the location
#' @param voxel_grid Dataframe containing the voxel grid with density values
#' @param F_sky_diff_init Initial diffuse shortwave radiation (W/m2)
#' @param F_sky_dir_init Initial direct shortwave radiation (W/m2)
#' @param omega_g_v Ground reflectance for vertical RTM
#' @param Kd_v Diffuse extinction coefficient for vertical RTM
#' @param Kb_v Direct extinction coefficient for vertical RTM
#' @param betad Fraction of scattered diffuse radiation in the backward direction
#' @param beta0 Fraction of scattered direct radiation in the backward direction
#' @return Dataframe containing the combined results of the shortwave RTM
#' A lot of these parameters are global parameters and actually do not need to be stated as arguments.
#' @importFrom suncalc getSunlightPosition
#' @importFrom data.table as.data.table setorder
#' @export
shortwave_two_stream_RTM <- function(datetime, lat, lon, voxel_grid,
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
  above_horizon <- as.integer(sun_altitude >= 0) # This is 1 if altitude is between 0 and pi/2 radians

  F_sky_diff_v <- F_sky_diff_init
  F_sky_dir_v  <- F_sky_dir_init * above_horizon

  # pre-allocation
  F_d_down <- numeric(n)
  F_d_up   <- numeric(n)
  F_b_down <- numeric(n)
  net_sw   <- numeric(n)

  # precompute columns
  xy_id <- X + Y * max(X)
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



