PRO IWV_CALC, q, z, iwv

;calculate IWV
;input: 
;q - abs. hum. (SI)
;z - height above ground (SI)
;output:
;iwv in kgm^-2

iwv = 0.
n = N_ELEMENTS(z)
FOR i_iwv = 0, n-2 DO BEGIN
 iwv = iwv + 0.5*(q(i_iwv)+q(i_iwv+1))*(z(i_iwv+1)-z(i_iwv))
ENDFOR

RETURN
END