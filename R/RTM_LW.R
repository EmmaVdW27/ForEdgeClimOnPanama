#' Two-stream RTM for longwave radiation and for a single column existing of layers/voxels with density values (%).
#' This model is inspired by the longwave RTM of the Ecosystem Demography model (ED 2.2).
#' @param soil_reflect Soil albedo (%)
#' @param density Density of each layer (%)
#' @param F_sky_lw Longwave radiation at top canopy (W/m2)
#' @param temperature Forest surface temperature of each layer (K)
#' @param T_soil Soil temperature at 8 cm depth (K)
#' @param T_macro Macrotemperature outside the forest (K)
#' @param Kd Longwave extinction coefficient, per unit density
#' @param omega Longwave scattering coefficient
#' @param beta Fraction of scattering in the backward direction for longwave radiation (%)
#' @return Dataframe containing the combined results of the 2D longwave RTM
#' @export

lw_two_stream <- function(soil_reflect, density, F_sky_lw, temperature, T_soil, T_macro, Kd, omega, beta) {
  omega_g = soil_reflect
  F_sky_diff = F_sky_lw

  # Calculate black body emission at the temperature of the layer.
  # Black body emission of atmosphere is defined to be 0 because F_sky_lw has already been defined as above canopy input.
  black = c(e_forest*sigma_SB * temperature^4, 0)
  # Black body emission of the soil
  # Emissivity of soil = absorption of soil for lw radiation => emissivity = 1 - omega_g
  black_g = (1-omega_g)*sigma_SB*T_soil^4

  # Number of (cohort) layers
  n <- length(density)
  # Add atmosphere layer
  density <- c(density, 0)
  # Inverse optical depths -> moved to run_foredgeclim()
  # mu atmosphere (density = 0, no forest structure) is defined to be 1
  # mu <- ifelse(density == 0, 1, 1 / (Kd * density))  # Inverse optical depth diffuse radiation

  # Initialization matrices
  S <- matrix(0, nrow = 2 * n + 2, ncol = 2 * n + 2)
  y <- numeric(2 * n + 2)

  # # Auxiliary variables -> moved to run_foredgeclim()
  # # Gamma
  # gamma_plus <<- rep(0.5 * (1 + sqrt((1 - omega) / (1 - (1 - 2 * beta) * omega))) , n)
  # gamma_plus[n+1] <<- 1
  # gamma_minus <<- rep(0.5 * (1 - sqrt((1 - omega) / (1 - (1 - 2 * beta) * omega))) , n)
  # gamma_minus[n+1] <<- 0
  # # Lambda
  # lambda <<- sqrt((1 - (1 - 2 * beta) * omega) * (1 - omega) / (mu)^2)
  # lambda[n+1] <<- 0
  # no Kappa and Delta as in shortwave RTM

  # Boundary conditions
  y[2*n+2] <- F_sky_diff #- black[n+1] # Downward (incoming) longwave radiation above canopy.
  y[1] <- black_g - (1 - omega_g)*black[1] # Upward (outgoing) longwave radiation from the ground surface.

  S[2 * n + 2, 2 * n + 1] <- gamma_plus_lw[n+1]
  S[2 * n + 2, 2 * n + 2] <- gamma_minus_lw[n+1]
  S[1, 1] <- (gamma_minus_lw[1] - omega_g*gamma_plus_lw[1])*exp(-lambda_lw[1]*density[1])
  S[1, 2] <- (gamma_plus_lw[1] - omega_g*gamma_minus_lw[1])*exp(lambda_lw[1]*density[1])

  # Vectorized matrix filling
  indices <- 1:n
  idx_even <- 2 * indices
  idx_odd <- idx_even - 1
  idx_next <- indices + 1

  # Calculate exponential terms
  exp_lambda_density_neg <- exp(-lambda_lw[idx_next] * density[idx_next])
  exp_lambda_density_pos <- exp(lambda_lw[idx_next] * density[idx_next])

  # Populate the S matrix
  S[cbind(idx_even, idx_odd)] <- gamma_plus_lw[indices]
  S[cbind(idx_even, idx_even)] <- gamma_minus_lw[indices]
  S[cbind(idx_even, idx_even + 1)] <- -gamma_plus_lw[idx_next] * exp_lambda_density_neg
  S[cbind(idx_even, idx_even + 2)] <- -gamma_minus_lw[idx_next] * exp_lambda_density_pos

  S[cbind(idx_even + 1, idx_odd)] <- gamma_minus_lw[indices]
  S[cbind(idx_even + 1, idx_even)] <- gamma_plus_lw[indices]
  S[cbind(idx_even + 1, idx_odd + 2)] <- -gamma_minus_lw[idx_next] * exp_lambda_density_neg
  S[cbind(idx_even + 1, idx_even + 2)] <- -gamma_plus_lw[idx_next] * exp_lambda_density_pos

  # Populate the vector y
  y[idx_even] <- black[idx_next] - black[indices]
  y[idx_even + 1] <- black[idx_next] - black[indices]


  # Solving matrix equation using arma::solve from C++
  x <- solve_rtm(S, y)

  # Calculate longwave radiative fluxes
  # Radiative fluxes refer to the rate of energy transfer by electromagnetic radiation through a surface or area,
  # typically measured in W/m², encompassing all or specific wavelengths (e.g., shortwave or longwave radiation).
  # These values are mainly of interest when modelling e.g. energy balance, heat fluxes...
  # Every down[k] and up[k] respectively refers to the downward and upward flux below layer k.
  indices <- seq_len(n + 1)
  k2 <- 2 * indices
  k2m1 <- k2 - 1

  down <- x[k2m1] * gamma_plus_lw[indices] * exp(-lambda_lw[indices] * density[indices]) +
    x[k2] * gamma_minus_lw[indices] * exp(lambda_lw[indices] * density[indices]) +
    black[indices]

  up <- x[k2m1] * gamma_minus_lw[indices] * exp(-lambda_lw[indices] * density[indices]) +
    x[k2] * gamma_plus_lw[indices] * exp(lambda_lw[indices] * density[indices]) +
    black[indices]

  # Calculate light levels and netto absorbed or emitted radiation
  # Light levels refer to the intensity of photosynthetically active radiation (PAR, 400-700 nm) available for plants.
  # These values are mainly of interest when modelling e.g. photosynthesis, plant growth, shadow structure...
  indices2 <- seq_len(n)
  kp1 <- indices2 + 1

  #light_diff_level <- 0.5 * (down[indices2] + down[kp1])
  net_lw <- (down[kp1] - down[indices2]) + (up[indices2] - up[kp1])

  # There are no light levels for atmosphere
  #light_diff_level <- c(light_diff_level, NA)

  #result <- data.frame(
  return(list(
    #layer = c(paste("layer", 1:n), "atmosphere"),
    density = density,
    temperature = c(temperature, NA),
    F_d_down = round(down, 0),
    F_d_up = round(up, 0),
    #BB_emission = black,
    net_lw = c(round(net_lw, 0), 0 )  # Net radiation in atmosphere is set to zero (no structure here)
  ))


}


#' 1D (vertical) longwave two-stream RTM
#'
#' This function performs a longwave
#' two-stream RTM to estimate radiation fluxes in a forest voxel grid.
#'
#' @param voxel_grid Dataframe containing voxel density data (columns: X, Y, Z, density).
#' @param micro_grid Dataframe containing microclimate temperature data (columns: x, y, z, temperature, T_soil).
#' @param lw_two_stream Function that calculates longwave two-stream radiation for a single column.
#' @param F_sky_lw Longwave sky flux.
#' @param omega_g_lw_v Ground reflectance for vertical LW RTM.
#' @param macro_temp Macroclimate temperature outside the forest.
#' @param Kd_lw_v Diffuse extinction coefficient for vertical LW RTM.
#' @param omega_lw Single-scattering albedo for longwave RTM.
#' @param beta_lw Phase function parameter for longwave RTM.
#' @return A dataframe with computed longwave fluxes and net longwave radiation.
#' A lot of these parameters are global parameters and actually do not need to be stated as arguments.
#' @export

longwave_two_stream_RTM <- function(voxel_grid, micro_grid, lw_two_stream,
                                         F_sky_lw, omega_g_lw_v, macro_temp,
                                         Kd_lw_v, omega_lw, beta_lw) {

  # create vectors, no data tables
  X <- voxel_grid$X
  Y <- voxel_grid$Y
  Z <- voxel_grid$Z
  density <- voxel_grid$density

  temperature <- micro_grid$temperature
  T_soil <- micro_grid$T_soil

  n <- length(X)

  # pre-allocation
  F_d_down <- numeric(n)
  F_d_up   <- numeric(n)
  net_lw   <- numeric(n)

  # precompute columns
  xy_id <- X + Y * max(X)
  split_idx <- split(seq_len(n), xy_id)

  for (idx in split_idx) {

    # extract column
    density_col <- density[idx]
    temp_col    <- temperature[idx]
    # single value for soil temp
    T_soil_value <- NA_real_
    for (j in idx) {
      if (!is.na(T_soil[j])) {
        T_soil_value <- T_soil[j]
        break
      }
    }

    ##########################
    # Vertical longwave RTM  #
    ##########################

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

  return(data.frame(
    X = X,
    Y = Y,
    Z = Z,
    density = density,
    temperature = temperature,
    F_d_down = F_d_down,
    F_d_up = F_d_up,
    net_lw = net_lw
  ))

}


