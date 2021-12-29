;+
;********************
PRO READ_LEVEL2B_NC,$
;********************
;INPUT:
filename,$                      ;name of netcdf file to be read
;OUTPUT:
algo,$                          ;retrieval filename
comment,$                       ;any specific comments
time,$                          ;array of decimal hours
z,$                             ;array of heights
tprof,$                         ;2D array of temperature profile (time x height)
qprof,$                         ;2D array of humidity profile (time x height)
dqprof,$                        ;2D array of humidity profile offset correction (org-oc) based on TB (time x height)
flag,$                          ;data flagging
temp,$                          ;surface temperature in K
pres,$                          ;surface pressure in hPa
relh                            ;surface relative humidity in %
; $Id: read_hatpro_level2b_nc.pro,v 1.2 2010/12/13 15:17:32 hatpro Exp $
; Abstract:
; * read HATPRO level2b netcdf files
; Authors:
; U. Loehnert
; Date:
; 2008-09-25
; Dependencies:
; -
; Changes:
; XXXX-XX-XX:
; changed program ...
;-

infid = NCDF_OPEN(filename(0), /NOWRITE)
tid = NCDF_VARID(infid,'time')
latid  = NCDF_VARID(infid, 'lat')
longid = NCDF_VARID(infid, 'lon')
altid = NCDF_VARID(infid, 'zsl')
zid = NCDF_VARID(infid, 'height')
tprofid = NCDF_VARID(infid, 'ta')
dtprofid = NCDF_VARID(infid, 'ta_offset')
qprofid = NCDF_VARID(infid, 'hua')
dqprofid = NCDF_VARID(infid, 'hua_offset')
rid = NCDF_VARID(infid, 'flag')
  
NCDF_VARGET, infid, latid, latitude
NCDF_VARGET, infid, longid, longitude
NCDF_VARGET, infid, altid, altitude
NCDF_VARGET, infid, tid, time
NCDF_VARGET, infid, zid, z
NCDF_VARGET, infid, tprofid, tprof
NCDF_VARGET, infid, dtprofid, dtprof
NCDF_VARGET, infid, qprofid, qprof
NCDF_VARGET, infid, dqprofid, dqprof
NCDF_VARGET, infid, rid, flag

; get 1st attribute name = NCDF_ATTNAME(infid,/GLOBAL,1)
NCDF_ATTGET,infid,/GLOBAL,'comments',com
comment = STRING(BYTE(com))

END
