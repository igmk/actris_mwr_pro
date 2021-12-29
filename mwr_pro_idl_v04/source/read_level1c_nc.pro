;+
;***************************
PRO READ_LEVEL1C_NC,$
;***************************
;INPUT:
filename,$                      ;name of netcdf file to be read
;OUTPUT:
comment,$                       ;any specific comments
latitude,longitude,altitude,$   ;three FLOAT
angles,$                        ;array of elevation angles
freq,$                          ;array of HATPRO frequencies
time,$                          ;array of decimal hours
tb,$                            ;3D TB array (freq x angles x time)
flag,$                          ;data flagging
temp,$                          ;surface temperature in K
pres,$                          ;surface pressure in hPa
relh                            ;surface relative humidity in %
; $Id: read_hatpro_level0c_nc.pro,v 1.2 2010/12/27 17:35:32 loehnert Exp $
; Abstract:
; * read HATPRO level1c netcdf files
; Authors:
; U. Loehnert
; Date:
; 2008-09-19
; Dependencies:
; -
; Changes:
; XXXX-XX-XX:
; changed program ...
;-

infid  = NCDF_OPEN(filename, /NOWRITE)
tid    = NCDF_VARID(infid,'time')
fid    = NCDF_VARID(infid,'frequencies')
tbid   = NCDF_VARID(infid,'brightness_temperature')
latid  = NCDF_VARID(infid, 'latitude')
longid = NCDF_VARID(infid, 'longitude')
altid  = NCDF_VARID(infid, 'altitude')
anglid = NCDF_VARID(infid, 'elevation_angle')
tid    = NCDF_VARID(infid, 'time')
fid    = NCDF_VARID(infid, 'frequencies')
rid    = NCDF_VARID(infid, 'flag')
teid   = NCDF_VARID(infid, 'air_temperature')
pid    = NCDF_VARID(infid, 'air_pressure')
rhid   = NCDF_VARID(infid, 'relative_humidity')
  
NCDF_VARGET, infid, latid,  latitude
NCDF_VARGET, infid, longid, longitude
NCDF_VARGET, infid, altid,  altitude
NCDF_VARGET, infid, anglid, angles
NCDF_VARGET, infid, tid,    time
NCDF_VARGET, infid, fid,    freq
NCDF_VARGET, infid, tbid,   tb
NCDF_VARGET, infid, rid,    flag
NCDF_VARGET, infid, teid,   temp
NCDF_VARGET, infid, pid,    pres
NCDF_VARGET, infid, rhid,   relh

; get 1st attribute name = NCDF_ATTNAME(infid,/GLOBAL,1)
NCDF_ATTGET,infid,/GLOBAL,'comments',com
comment = STRING(BYTE(com))

NCDF_CLOSE, infid

n = N_ELEMENTS(time)

END
