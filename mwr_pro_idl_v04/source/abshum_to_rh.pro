;+
;*****************
PRO ABSHUM_TO_RH,$
;*****************
;INPUT
T,$                    ;temperature in K 
rho_w,$                ;absolute humidity in kg/m^3
;OUTPUT
rh,$                   ;relative humidity in %
;KEYWORDS:
verbose=verbose
; Abstract:
;* calculate relative humidity from absolute humidity and temperature
;  using analytical approximation of Clausius-Clapeyron equation
; Author:
; U. Loehnert
; Date:
; 2011-11-09
; Dependencies:
; -   
; Changes:
;-

;**specific gas constant for water vapor (J/kg K)
Rw = 462.

;**vapor pressure e0 (Pa) at T0 (K)
T0 = 273.15
e0 = 611.

;**specific heat for evaporation (J/kg)
L = (2500.-2.42*(T-273.15))*1000.
es = e0*EXP((L/(Rw*T0))*((T-T0)/T)) ; saturation pressure in Pa 
e = rho_w*Rw*T ; water vapor pressure
rh = e/es 
rh = rh*100. ; relative humidity
END
