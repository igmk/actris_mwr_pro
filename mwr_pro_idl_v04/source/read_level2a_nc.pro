;+
;********************
PRO READ_LEVEL2A_NC,$
;********************
;INPUT:
filename,$                      ;name of netcdf file to be read
;OUTPUT:
algo,$                          ;retrieval filename
comment,$                       ;any specific comments
time,$                          ;array of decimal hours
ele,$                           ;array of elevation angles 
azi,$                           ;array of azimuth angles
iwv,$                           ;array integrated water vapor (time)
lwp,$                           ;array of liquid water path (time)
flag,$                          ;data flagging
temp,$                          ;surface temperature in K
pres,$                          ;surface pressure in hPa
relh                            ;surface relative humidity in %
; $Id: $
; Abstract:
; * read HATPRO level2a netcdf files
; Authors:
; U. Loehnert
; Date:
; 2012-09-07
; Dependencies:
; -
; Changes:
; XXXX-XX-XX:
; changed program ...
;-

infid  = NCDF_OPEN(filename(0), /NOWRITE)
tid    = NCDF_VARID(infid,'time')
latid  = NCDF_VARID(infid, 'latitude')
longid = NCDF_VARID(infid, 'longitude')
altid  = NCDF_VARID(infid, 'altitude')
elid   = NCDF_VARID(infid, 'elevation_angle')
azid   = NCDF_VARID(infid, 'azimuth_angle')
iwvid  = NCDF_VARID(infid, 'atmosphere_water_vapor_content')
lwpid  = NCDF_VARID(infid, 'atmosphere_liquid_water_content')
rid    = NCDF_VARID(infid, 'flag')
teid   = NCDF_VARID(infid, 'air_temperature')
pid    = NCDF_VARID(infid, 'air_pressure')
rhid   = NCDF_VARID(infid, 'relative_humidity')
  
NCDF_VARGET, infid, latid,  latitude
NCDF_VARGET, infid, longid, longitude
NCDF_VARGET, infid, altid,  altitude
NCDF_VARGET, infid, tid,    time
NCDF_VARGET, infid, elid,   ele
NCDF_VARGET, infid, azid,   azi
NCDF_VARGET, infid, iwvid,  iwv
NCDF_VARGET, infid, lwpid,  lwp
NCDF_VARGET, infid, rid,    flag
NCDF_VARGET, infid, teid,   temp
NCDF_VARGET, infid, pid,    pres
NCDF_VARGET, infid, rhid,   relh

; get 1st attribute name = NCDF_ATTNAME(infid,/GLOBAL,1)
NCDF_ATTGET,infid,/GLOBAL,'comments',com
comment = STRING(BYTE(com))

END
