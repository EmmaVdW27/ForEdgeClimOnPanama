#' Function to calculate the derivative of the energy balance closure error to surface temperature
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


#' Function to calculate the derivative of the energy balance closure error to surface temperature with PM LE
#'
#' @param temp Forest surface temperature (K)
#' @param T_air_K Air temperature in Kelvin (K)
#' @param net_rad Net radiation (W/m2)
#' @param G Ground heat flux (W/m2)
#' @param gamma_eff effective psychrometric constant (kPa/K)
#' @return Derivative of energy balance wrt surface temperature
de_dT_PM <- function(temp, T_air_K, net_rad, G, gamma_eff) {
  # 1. conversion to Celsius for s-calculation
  T_c = temp - 273.15
  T_air_c = T_air_K - 273.15

  # 2. derivative of radiation and sensible heat
  dR_dT = -4 * e_forest * sigma_SB * temp^3
  dH_dT = den * g_forest

  # 3. saturation vapour pressure and derivative
  # es is the saturation vapour pressure at surface temperature
  es = 0.6108 * exp(17.27 * T_c / (T_c + 237.3))
  s = 4098 * es / (T_c + 237.3)^2

  # 4. actual vapour pressure of air
  # ea is constant with respect to the surface temperature
  ea = 0.6108 * exp(17.27 * T_air_c / (T_air_c + 237.3)) * (RH / 100)

  # 5. second derivative of vapour pressure
  ds_dT = s * (17.27 / (T_c + 237.3) - (17.27 * T_c) / (T_c + 237.3)^2 - 2 / (T_c + 237.3))

  # 6. Penman-Monteith quotient rule: (f'g - fg') / g^2
  Noemer = s + gamma_eff
  dNoemer_dT = ds_dT

  VPD_term = rho * Cp * g_forest * (es - ea)
  Teller = s * (net_rad - G) + VPD_term
  dTeller_dT = (ds_dT * (net_rad - G) + s * dR_dT) + (rho * Cp * g_forest * s)

  dLE_dT = den * ((dTeller_dT * Noemer) - (Teller * dNoemer_dT)) / (Noemer^2)

  # 7. total derivative of energy balance error
  return(dR_dT - dH_dT - dLE_dT)
}


