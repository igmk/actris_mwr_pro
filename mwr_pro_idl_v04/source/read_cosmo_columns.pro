PRO READ_COSMO_COLUMNS, station, date, time, z, T, p, q

if station eq 'jue' or station eq 'sel' then station_str = 'Selhausen_7x7'

year = '20'  + strmid(date, 0, 2)

cosmo_path = '/data/gop/final/gop9/cols_lmk/'  + year + '/'

filename  = '20' + date + '_gop9_lmk_' + station_str + '.nc'

tarname = cosmo_path + '20' + date + '_gop9_lmk_cols.tar'
spawn, 'ls ' + cosmo_path + '20' + date + '_gop9_lmk_cols.tar', list
if list ne tarname then begin
    print, 'no cosmo profile available'
    z = -999.
    goto, ABORT
endif

spawn, 'tar -xvf '  + cosmo_path + '20' + date + '_gop9_lmk_cols.tar '  + filename + '.gz'
spawn, 'gunzip ' + filename + '.gz'

infid  = NCDF_OPEN(filename, /NOWRITE)

time      = NCDF_VARID(infid, 'time')  
longitude = NCDF_VARID(infid, 'longitude')  
latitude  = NCDF_VARID(infid, 'latitude')  
z 	  = NCDF_VARID(infid, 'hfl')  
T 	  = NCDF_VARID(infid, 'temperature')  
p         = NCDF_VARID(infid, 'p')  
q         = NCDF_VARID(infid, 'qv')  
 
NCDF_VARGET, infid, time, time
NCDF_VARGET, infid, longitude, longitude
NCDF_VARGET, infid, latitude, latitude
NCDF_VARGET, infid, z, z
NCDF_VARGET, infid, T, T
NCDF_VARGET, infid, p, p
NCDF_VARGET, infid, q, q

spawn, 'rm ' + filename + '*'

ABORT:

END