;***********
PRO READ_SCANCLDRAD_NC,$
;***********
;INPUT:
filename,$	
;OUPUT:
time,$		;unix-time(s since 1970-01-01 00UTC)
range,$		;range (from antenna) [m]
Ze,$		;equivalent reflectivity 
ldr,$
vel,$ 
rms,$
elv,$
elvv,$
azi,$
aziv


infid   = NCDF_OPEN(filename, /NOWRITE)

;get variables
NCDF_VARGET, infid, NCDF_VARID(infid, 'time'), time
NCDF_VARGET, infid, NCDF_VARID(infid, 'range'), range
NCDF_VARGET, infid, NCDF_VARID(infid, 'Ze'), Ze
NCDF_VARGET, infid, NCDF_VARID(infid, 'LDR'), ldr
NCDF_VARGET, infid, NCDF_VARID(infid, 'VEL'), vel
NCDF_VARGET, infid, NCDF_VARID(infid, 'RMS'), rms
NCDF_VARGET, infid, NCDF_VARID(infid, 'elv'), elv
NCDF_VARGET, infid, NCDF_VARID(infid, 'elvv'), elvv
NCDF_VARGET, infid, NCDF_VARID(infid, 'azi'), azi
NCDF_VARGET, infid, NCDF_VARID(infid, 'aziv'), aziv

NCDF_CLOSE, infid

END
