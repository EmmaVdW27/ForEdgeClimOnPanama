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
#' @param net_rad Net radiation (W/m2)
#' @param den Voxel density (0-1)
#' @return Derivative of energy balance wrt surface temperature
de_dT_PM <- function(temp, T_air_K, net_rad, G, gamma_eff) {
  # 1. Omzetting naar Celsius voor s-berekening
  T_c = temp - 273.15
  T_air_c = T_air_K - 273.15

  # 2. Straling en Voelbare warmte afgeleiden
  dR_dT = -4 * e_forest * sigma_SB * temp^3
  dH_dT = den * g_forest

  # 3. Verzadigingsdampdruk en de afgeleide (s)
  # es is de verzadigingsdruk bij de oppervlakte-temperatuur (temp)
  es = 0.6108 * exp(17.27 * T_c / (T_c + 237.3))
  s = 4098 * es / (T_c + 237.3)^2

  # 4. Actuele dampdruk van de lucht (ea) - DIT ONTBAK
  # ea is constant t.o.v. de oppervlaktetemperatuur in deze stap
  ea = 0.6108 * exp(17.27 * T_air_c / (T_air_c + 237.3)) * (RH / 100)

  # 5. Tweede afgeleide van de verzadigingscurve (ds_dT)
  # Nodig voor de quotiëntregel van de PM-vergelijking
  ds_dT = s * (17.27 / (T_c + 237.3) - (17.27 * T_c) / (T_c + 237.3)^2 - 2 / (T_c + 237.3))

  # 6. Penman-Monteith Quotiëntregel: dLE/dT
  # LE = den * [ s*(Rn - G) + rho*Cp*ga*(es - ea) ] / [ s + gamma_eff ]

  Noemer = s + gamma_eff
  dNoemer_dT = ds_dT

  # Teller en de afgeleide van de teller
  VPD_term = rho * Cp * g_forest * (es - ea)
  Teller = s * (net_rad - G) + VPD_term

  # dTeller/dT: d/dT van [s*(Rn-G)] + d/dT van [VPD_term]
  # d/dT(es - ea) = s (omdat ea constant is t.o.v. T_surf)
  dTeller_dT = (ds_dT * (net_rad - G) + s * dR_dT) + (rho * Cp * g_forest * s)

  # De quotiëntregel: (f'g - fg') / g^2
  dLE_dT = den * ((dTeller_dT * Noemer) - (Teller * dNoemer_dT)) / (Noemer^2)

  # 7. Totale afgeleide van de energiebalans-fout
  # f(T) = Rn - H - LE - G (G is hier vaak constant verondersteld)
  return(dR_dT - dH_dT - dLE_dT)
}


