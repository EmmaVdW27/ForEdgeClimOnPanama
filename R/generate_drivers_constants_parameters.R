#' Generate input driver variables
#'
#' @return input driver variables
#' @importFrom readxl read_excel
#' @export
create_input_drivers <- function() {

  # input files
  temp_input_file_rds <<- readRDS("Data/TMS_macro_2023_4C_4R.rds")
  rad_input_file_csv <<- "Data/Table_CNR1_2023.dat"
  PAR_input_file_rds <<- readRDS("Data/PAR_data_full_2023.rds")
  sensor_locations <<- read_excel("Data/sensor_locations_all_plots.xlsx")
  RH_input_file_csv <<- "Data/bci_cl_rh_elect.csv"
  u_input_file_csv <<- "Data/bci_ava_ws_elect.csv"

  # spacetime specifics
  lat <<- 9.108
  lon <<- -79.851

  # parameters for energy balance convergence
  energy_balance_tolerance <<- 2  # Maximum energy balance closure error between successive iteration steps to reach convergence (W/m2)

}

#' Generate (physical) constants
#'
#' @return (physical) constants
#' @export
create_physical_constants <- function() {

  # constants for longwave RTM
  sigma_SB <<- 5.67e-8  # Stefan-Boltzmann constant (W/m²/K⁴)

  # constants for air to air heat convection (C)
  voxel_length <<- 1    # Voxel edge length (m)
  Cp <<- 1000           # Specific heat air (J/kg/K)
  rho <<- 1.225         # Density air (kg/m3)

  # constants for ground heat flux (G) and soil temperature
  stable_soil_depth <<- 0.08 # Depth of measured soil temperature

  # constants (empirical) for latent heat flux (LE)
  alpha_PT <<- 1.26         # Replaces aerodynamic terms from Penman-Monteith equation (unitless)
                            # 1.26 is the value for non-water-stressed conditions (arid regions lower, wetlands higher)
  gamma_psy <<- 0.066       # Psychrometric constant (kPa/K = kPa/°C)

}


#' Generate ForEdgeClim model parameters
#'
#' @return model parameters
#' @export
#
# These model parameters are the uncalibrated values as found in literature.
# I.e., these are the averages of the uniform distributions as found in literature.
create_model_parameters <- function() {

  # parameters for shortwave RTM
  betad <<- 0.325         # Fraction of scattered diffuse radiation in backward direction
  beta0 <<- 0.325         # Fraction of scattered direct beam radiation in backward direction
  omega <<- 0.52          # Shortwave scattering coefficient
  # -> Vertical RTM
  Kd_v <<- 0.775          # Diffuse extinction coefficient for vertical radiation, per unit density
  Kb_v <<- 1.25           # Direct beam extinction coefficient for vertical radiation, per unit density
  omega_g_v <<- 0.13      # Ground scattering for shortwave radiation

  # parameters for longwave RTM
  e_forest <<- 0.965      # Emissivity of the forest (leaves and wood)
  beta_lw <<- 0.325       # Fraction of scattered longwave radiation in backward direction
  omega_lw <<- 0.035      # Longwave scattering coefficient
  # -> Vertical RTM
  Kd_lw_v <<- 0.3         # Longwave extinction coefficient for vertical radiation, per unit density
  omega_g_lw_v <<- 0.055  # Ground scattering for longwave radiation

  # parameters for air to air heat convection (C)
  h <<- 10                # Convection coefficient of air (W/m2/K), higher value means more wind/more turbulence

  # parameters for update linearization air temperature
  g_macro_ref <<- 25      # Convective heat transfer coefficient between (voxel) air and (macro) air (W/m2/K)
  infl_macro_ref <<- 32.5      # Distance over which the influence of macro temp on air temp is reduced by 50% (m)
  infl_soil_ref <<- 5          # Distance over which the influence of soil temp on air temp is reduced by 50% (m)
  infl_forest_ref <<- 5       # Distance over which the influence of forest temp on air temp is reduced by 50% (m)

  # parameters for sensible heat flux (H)
  g_forest_ref <<- 12.5       # Combined conductive & convective heat transfer coefficient between (voxel) air and structure (leaf) (W/m2/K)

  # parameters for ground heat flux (G) and soil temperature
  p_ground <<- 0.225      # Fraction of net ground radiation to define ground flux
  g_soil_ref <<- 10           # Convective heat transfer coefficient between (voxel) air of air layer just above the ground and ground surface (W/m2/K)
  k_soil <<- 1.225        # Thermal conductance soil (W/m/K)


  # PM LE parameters
  gs_max <<- 0.004  # maximal stomatal conductance (mol/m2/s)
  stomata_closure <<- 0.75 # stomatal closure parameter that controls the sensitivity of stomatal conductance to increasing VPD (dimensionless)

  # g(u) parameters
  k_u <<- 0.7 # wind attenuation coefficient (dimensionless)
  u_ref <<- 1.6 # typical wind speed around forests (above and below canopy) (m/s)
}

