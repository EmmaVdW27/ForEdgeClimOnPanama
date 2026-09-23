#' Function to calculate the derivative of the energy balance closure error to surface temperature with Priestley-Taylor approximation.
#'
#' @param temp Forest surface temperature
#' @param net_rad Net radiation
#' @return Derivative of energy balance wrt forest surface temperature
#' @export
#'
de_dT <- function(temp, net_rad){
  # derivative of net radiation to temperature:
  dR_dT = -4*e_forest*sigma_SB*temp^3
  # derivative of sensible heat transfer to temperature:
  dH_dT = den*g_forest
  # derivative of latent heat flux to temperature:
  T_c = temp - 273.15 # temperature in degrees Celsius
  s = 4098 * saturated_vapor_pressure(T_c) / (T_c + 237.3)^2 # slope in Priestly-Taylor equation
  des_dT = 0.6108 * exp(17.27*T_c/(T_c+237.3)) * ((T_c+237.3)*17.27 - 17.27*T_c)/(T_c+237.3)^2 # derivative of saturated vapor pressure to temperature
  ds_dT = (4098*des_dT*(T_c+237.3)^2 - 2*(T_c+237.3)*4098*saturated_vapor_pressure(T_c))/(T_c+237.3)^4 # derivative of slope to temperature
  df_dT = (ds_dT*(s+gamma_psy)-s*ds_dT)/(s+gamma_psy)^2 # derivative of s/(s+gamma_psy) to temperature
  dLE_dT = den*alpha_PT * (dR_dT)* s / (s + gamma_psy) + den*alpha_PT * (net_rad) * df_dT

  return(dR_dT - dH_dT - dLE_dT)

}

#' Derivative of the surface energy-balance residual with respect to
#' surface temperature using Penman-Monteith latent heat flux
#'
#' @param temp Surface temperature (K)
#' @param SW_down Downward shortwave radiation (W/m2)
#' @param net_rad Net radiation (W/m2)
#' @param G Ground heat flux (W/m2)
#' @return Derivative of the energy-balance residual with respect to
#'   surface temperature
#'
de_dT_PM <- function(temp, SW_down, net_rad, G) {

  # Surface temperature in degrees Celsius
  T_c <- temp - 273.15

  # ------------------------------------------------------------------
  # Radiation and sensible heat
  # ------------------------------------------------------------------

  # dRn / dTsurf
  dR_dT <- -4 * e_forest * sigma_SB * temp^3

  # dH / dTsurf, with Tair fixed during the Newton step
  dH_dT <- den * g_forest


  # ------------------------------------------------------------------
  # Saturation vapour pressure
  # ------------------------------------------------------------------

  es <- saturated_vapor_pressure(T_c)

  # des/dT = slope of saturation vapour pressure curve
  s <- 4098 * es / (T_c + 237.3)^2

  # ds/dT = d²es/dT²
  ds_dT <- s * (
    17.27 * 237.3 / (T_c + 237.3)^2 -
      2 / (T_c + 237.3)
  )


  # ------------------------------------------------------------------
  # Vapour pressure deficit
  # ------------------------------------------------------------------

  ea <- saturated_vapor_pressure(macro_temp - 273.15) * (RH / 100)

  VPD <- pmax(0, es - ea)

  dVPD_dT <- ifelse(VPD > 0, s, 0)


  # ------------------------------------------------------------------
  # Stomatal conductance
  # ------------------------------------------------------------------

  PAR <- 0.45 * 4.57 * SW_down

  f_light <- PAR / (PAR + 100)
  f_vpd   <- exp(-stomata_closure * VPD)

  gs_raw <- gs_max * f_light * f_vpd
  gs_eff <- pmax(0.0001, gs_raw)

  # Derivative is zero when gs is limited by its imposed minimum
  dgs_dT <- ifelse(
    VPD > 0 & gs_raw > 0.0001,
    -stomata_closure * gs_eff * s,
    0
  )


  # ------------------------------------------------------------------
  # Effective psychrometric constant
  # ------------------------------------------------------------------

  alpha_a = 1 # in m.K/(W.s)
  ga <- g_forest * alpha_a

  gamma_eff <- gamma_psy * (1 + ga / gs_eff)

  dgamma_eff_dT <-
    -gamma_psy * ga / gs_eff^2 * dgs_dT


  # ------------------------------------------------------------------
  # Penman-Monteith
  # ------------------------------------------------------------------

  numerator <-
    s * (net_rad - G) +
    rho * Cp * ga * VPD

  denominator <- s + gamma_eff

  # d numerator / dT
  #
  # net_rad changes with surface temperature;
  # G is treated as fixed within one Newton step.
  dnumerator_dT <-
    ds_dT * (net_rad - G) +
    s * dR_dT +
    rho * Cp * ga * dVPD_dT

  # d denominator / dT
  ddenominator_dT <-
    ds_dT + dgamma_eff_dT

  # Quotient rule
  dLE_dT <- den * (
    dnumerator_dT * denominator -
      numerator * ddenominator_dT
  ) / denominator^2


  # ------------------------------------------------------------------
  # Energy-balance residual
  #
  # E = Rn - H - LE
  # ------------------------------------------------------------------

  dE_dT <- dR_dT - dH_dT - dLE_dT

  return(dE_dT)
}


