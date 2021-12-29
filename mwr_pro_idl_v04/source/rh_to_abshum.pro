;+
;*****************
PRO RH_TO_ABSHUM,$
;*****************
;INPUT
T,$                    ;temperature in K 
rh,$                   ;relative humidity in %
;OUTPUT
rho_w,$                ;absolute humidity in kg/m^3
;KEYWORDS:
verbose=verbose
; Abstract:
;* calculate absolute humidity from relative humidity and temperature
;  using analytical approximation of Clausius-Clapeyron equation
; Author:
; U. Loehnert
; Date:
; 2014-10-31
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
rh = rh/100.
e = es*rh
rho_w = e/(Rw*T); absolute humidity
END
