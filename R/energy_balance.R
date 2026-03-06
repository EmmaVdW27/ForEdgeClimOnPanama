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


#' Aangepaste afgeleide van de energiebalans voor Penman-Monteith
#'
#' @param temp Forest surface temperature (K)
#' @param net_rad Net radiation (W/m2)
#' @param den Voxel density (0-1)
#' @param gs_base Stomatal conductance (m/s)
#' @return Derivative of energy balance wrt surface temperature
de_dT_PM <- function(temp, net_rad, den, gs_base = 0.005) {
  T_c = temp - 273.15
  dR_dT = -4 * e_forest * sigma_SB * temp^3
  dH_dT = den * g_forest

  # Vereenvoudigde maar stabiele dLE/dT
  # We focussen op de dominante slope-verandering
  s = 4098 * saturated_vapor_pressure(T_c) / (T_c + 237.3)^2
  gamma_eff = 0.066 * (1 + g_forest / gs_base)

  # De verandering in LE wordt 'gedempt' door de weerstand
  dLE_dT = den * (s / (s + gamma_eff)) * dR_dT

  # Voeg een kleine epsilon toe om deling door exact nul te voorkomen
  return(dR_dT - dH_dT - dLE_dT - 1e-6)
}
