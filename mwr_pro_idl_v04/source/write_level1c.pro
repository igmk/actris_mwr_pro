;+
;******************
PRO WRITE_LEVEL1C,$
;******************[B
;INPUT:
out_path,$                      ;path where .nc file will be written
date,$                          ;YYMMDD (STRING)
measurement,$                   ;three letter code identifying the measurement (STRING)
lat, lon, zsl,$                 ;three FLOAT
freq,$                          ;array of MWR + IRT frequencies
time,$                          ;array of decimal hours
time_bnds,$                     ;time array bounds
tb,$                            ;3D TB array dimension = N(freq)xN(angles)xN(time)
offset_tb,$                     ;3D TB array of TB offset if specified
offset_tb_type,$                ;CHARACTER defining source for offset correction
freq_shift,$                    ;frequency shift applied to original TBs by RPG software
el,$                            ;array of elevation angles (angles)
az,$                            ;array of azimuth angles (angles)
flag,$                          ;data flagging, see attributes for details
temp,$                          ;surface temperature in K (if available)
pres,$                          ;surface pressure in hPa (if available)
relh,$                          ;surface relative humidity in % (if available)
threshold,$                     ;TB thresholds
fn,$                            ;file naming parameters hdcp2 compliant
ga,$                            ;global attributes 
tb_bias,$                       ;systematic TB calibration uncertainty, 1 stddev
tb_cov                          ;TB error covariance matrix of MWR channels (random uncertainty)
; $Id: $
; Abstract:
; * write HATPRO level1c TB files to netcdf format
; Authors:
; U. Loehnert
; Date:
; 2011-09-011
; Dependencies:
; Changes:
; 2013-09-09 (UL)
; make written level1c file HDCP(2) compliant  
; 2015-01-08 (UL)
; pressure now written to file as Pa (formerly hPa)
; -

date_create = SYSTIME(0, /UTC)

;****create date_string YYYYMMDDhhmmss 
year = '20'+STRMID(date,0,2)
month= STRMID(date,2,2)
day  = STRMID(date,4,2)
IF time(0) LT 10. THEN BEGIN
 hh_c = '0'+STRING(time(0), format = '(i1)')
ENDIF ELSE BEGIN
 hh_c = STRING(time(0), format = '(i2)')
ENDELSE
mm = time(0) - FIX(time(0))
mm = mm*60.
IF mm LT 10. THEN BEGIN
 mm_c = '0'+STRING(mm, format = '(i1)')
ENDIF ELSE BEGIN
 mm_c = STRING(mm, format = '(i2)')
ENDELSE
ss = mm - FIX(mm)
ss = ss*60.
IF ss LT 10. THEN BEGIN
 ss_c = '0'+STRING(ss, format = '(i1)')
ENDIF ELSE BEGIN
 ss_c = STRING(ss, format = '(i2)')
ENDELSE

date_string = year+month+day+hh_c+mm_c+ss_c 
date_long = LONG(date_string)

date_string = year+month+day+hh_c+mm_c+ss_c
date_long = LONG(date_string)

time_conv = time
time_bnds_conv = time_bnds

;****transform time to seconds since 19700101
SECONDS_SINCE_19700101, year+month+day, time_conv, time

;****transform time_bnds to seconds since 19700101
SECONDS_SINCE_19700101, year+month+day, time_bnds_conv, time_bnds
time_bnds = TRANSPOSE(time_bnds)

;****initiate dimensions
n_time = N_ELEMENTS(time)
n_freq = N_ELEMENTS(freq)
n_band = N_ELEMENTS(threshold)
n_angle = N_ELEMENTS(el)

;****create netcdf output files
outfile = out_path + fn.kkk + '_' + fn.sss + '_' + fn.inst + fn.instnn + '_'+ fn.lll1 + '_' + fn.var + '_' + fn.vnn + '_' + date_string + '.nc'
outfid = NCDF_CREATE(outfile, /clobber)

;****define dimensions
timdid = NCDF_DIMDEF(outfid, 'time', n_time)
nvdid = NCDF_DIMDEF(outfid, 'nv', 2)
fredid = NCDF_DIMDEF(outfid, 'n_freq', n_freq)
fredid2 = NCDF_DIMDEF(outfid, 'n_freq2', n_freq)
bandid = NCDF_DIMDEF(outfid, 'n_band', n_band)
angdid = NCDF_DIMDEF(outfid, 'n_angle', n_angle)

;****define variables & attributes
timid = NCDF_VARDEF(outfid, 'time',[timdid], /double)
NCDF_ATTPUT, outfid, timid, 'units', 'seconds since 1970-01-01 00:00:00 UTC'
NCDF_ATTPUT, outfid, timid, 'standard_name', 'time'
NCDF_ATTPUT, outfid, timid, 'bounds', 'time_bnds'

;tbnid = NCDF_VARDEF(outfid, 'time_bnds',[timdid, nvdid], /double)
tbnid = NCDF_VARDEF(outfid, 'time_bnds',[nvdid, timdid], /double)

latid = NCDF_VARDEF(outfid, 'lat', /float)
NCDF_ATTPUT, outfid, latid, 'units','degrees_north'
NCDF_ATTPUT, outfid, latid, 'long_name','latitude'

lonid = NCDF_VARDEF(outfid, 'lon', /float)
NCDF_ATTPUT, outfid, lonid, 'units','degrees_east'
NCDF_ATTPUT, outfid, lonid, 'long_name','longitude'

altid = NCDF_VARDEF(outfid, 'zsl', /float)
NCDF_ATTPUT, outfid, altid, 'units', 'm'
NCDF_ATTPUT, outfid, altid, 'standard_name','altitude'
NCDF_ATTPUT, outfid, altid, 'long_name','altitude above mean sea level'

freid = NCDF_VARDEF(outfid, 'freq_sb',[fredid], /float)
NCDF_ATTPUT, outfid, freid, 'units', 'GHz'
NCDF_ATTPUT, outfid, freid, 'standard_name','sensor_band_central_radiation_frequency'
NCDF_ATTPUT, outfid, freid, 'long_name','frequency of microwave channels'

aziid = NCDF_VARDEF(outfid, 'azi',[angdid], /float)
NCDF_ATTPUT, outfid, aziid, 'units', 'degree'
NCDF_ATTPUT, outfid, aziid, 'standard_name','sensor_azimuth_angle'
NCDF_ATTPUT, outfid, aziid, 'long_name','sensor azimuth angle'
NCDF_ATTPUT, outfid, aziid, '_FillValue',-999.
NCDF_ATTPUT, outfid, aziid, 'comment', '0=North, 90=East, 180=South, 270=West'

eleid = NCDF_VARDEF(outfid, 'ele',[angdid], /float)
NCDF_ATTPUT, outfid, eleid, 'units', 'degree'
NCDF_ATTPUT, outfid, eleid, 'long_name','sensor elevation angle'
NCDF_ATTPUT, outfid, eleid, '_FillValue',-999.

tbrid = NCDF_VARDEF(outfid, 'tb',[fredid,angdid,timdid], /float)
NCDF_ATTPUT, outfid, tbrid, 'units', 'K'
NCDF_ATTPUT, outfid, tbrid, 'standard_name', 'brightness_temperature'
NCDF_ATTPUT, outfid, tbrid, 'long_name', 'brightness temperatures'
NCDF_ATTPUT, outfid, tbrid, '_FillValue',-999.
NCDF_ATTPUT, outfid, tbrid, 'valid_min', threshold(0)
NCDF_ATTPUT, outfid, tbrid, 'valid_max', threshold(1)

tbocrid = NCDF_VARDEF(outfid, 'offset_tb',[fredid,angdid,timdid], /float)
NCDF_ATTPUT, outfid, tbocrid, 'units', 'K'
NCDF_ATTPUT, outfid, tbocrid, 'long_name', 'brightness temperatures offset correction subtracted'
NCDF_ATTPUT, outfid, tbocrid, '_FillValue',-999.
NCDF_ATTPUT, outfid, tbocrid, 'comment', 'Some types of MWR require a systemmatic adjustement of the measured TBs. This variable gives the offset which was subtracted from each measurement. The offset was deteremined from ' + offset_tb_type + '.'

fsrid = NCDF_VARDEF(outfid, 'freq_shift',[fredid], /float)
NCDF_ATTPUT, outfid, fsrid, 'units', '1'
NCDF_ATTPUT, outfid, fsrid, 'long_name', 'frequency shift applied to correct measured brightness temperature for frequency offset of microwave radiometer channel'
NCDF_ATTPUT, outfid, fsrid, '_FillValue',-999.

tbiid = NCDF_VARDEF(outfid, 'tb_bias',[fredid], /float)
NCDF_ATTPUT, outfid, tbiid, 'units', 'K'
NCDF_ATTPUT, outfid, tbiid, 'long_name', 'systematic calibration uncertainty of brightness temperature, one standard deviation'
NCDF_ATTPUT, outfid, tbiid, '_FillValue',-999.
NCDF_ATTPUT, outfid, tbiid, 'comment', 'This variable is an estimate of the one-standard-deviation calibration error to be expected from an absolute system calibration, i.e. the likely systematic error of brightness temperature. As a reference see Maschwitz et al. 2012, AMT (Tab. 5). However, these numbers differ from instrument to instrument and should be adapted accordingly. Values only valid for elevation angles larger than 20deg.'

tbcid = NCDF_VARDEF(outfid, 'tb_cov',[fredid, fredid], /float)
NCDF_ATTPUT, outfid, tbcid, 'units', 'K*K'
NCDF_ATTPUT, outfid, tbcid, 'long_name', 'error covariance matrix of brightness temperature channels'
NCDF_ATTPUT, outfid, tbcid, '_FillValue',-999.
NCDF_ATTPUT, outfid, tbcid, 'comment', 'This variable is calculated from brightness temperature observations of an internal black body whose physical temperature is known. The square root of the matrix diagonal gives the brightness temperature random error of each frequency channel. Values only valid for elevation angles larger than 20deg.' 

temid = NCDF_VARDEF(outfid, 'ta',[timdid], /float)
NCDF_ATTPUT, outfid, temid, 'units', 'K'
NCDF_ATTPUT, outfid, temid, 'standard_name','air_temperature'
NCDF_ATTPUT, outfid, temid, 'long_name','air temperature'
NCDF_ATTPUT, outfid, temid, '_FillValue',-999.
NCDF_ATTPUT, outfid, temid, 'comments', 'ambient air temperature measured by sensor on microwave radiometer'
NCDF_ATTPUT, outfid, temid, 'valid_min', 200.
NCDF_ATTPUT, outfid, temid, 'valid_max', 330.

preid = NCDF_VARDEF(outfid, 'pa',[timdid], /float)
NCDF_ATTPUT, outfid, preid, 'units', 'Pa'
NCDF_ATTPUT, outfid, preid, 'standard_name','air_pressure'
NCDF_ATTPUT, outfid, preid, 'long_name','air pressure'
NCDF_ATTPUT, outfid, preid, '_FillValue',-999.
NCDF_ATTPUT, outfid, preid, 'comments', 'ambient air pressure measured by sensor on microwave radiometer'
NCDF_ATTPUT, outfid, preid, 'valid_min', 90000.
NCDF_ATTPUT, outfid, preid, 'valid_max', 104000.

rhuid = NCDF_VARDEF(outfid, 'hur',[timdid], /float)
NCDF_ATTPUT, outfid, rhuid, 'units', '1'
NCDF_ATTPUT, outfid, rhuid, 'standard_name','relative_humidity'
NCDF_ATTPUT, outfid, rhuid, 'long_name','relative humidity'
NCDF_ATTPUT, outfid, rhuid, '_FillValue',-999.
NCDF_ATTPUT, outfid, rhuid, 'comments', 'ambient relative humidity measured by sensor on microwave radiometer'
NCDF_ATTPUT, outfid, rhuid, 'valid_min', 0.
NCDF_ATTPUT, outfid, rhuid, 'valid_max', 1.1

flaid = NCDF_VARDEF(outfid, 'flag',[timdid], /short)
NCDF_ATTPUT, outfid, flaid, 'long_name','quality control flags'
NCDF_ATTPUT, outfid, flaid, '_FillValue', 0
NCDF_ATTPUT, outfid, flaid, 'flag_masks', [1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
NCDF_ATTPUT, outfid, flaid, 'flag_meanings', "manual_filter_band_1 manual_filter_band2 manual_filter_band3 rain_flag sanity_receiver_band1 sanity_receiver_band2 sun_in_beam tb_threshold_band1 tb_threshold_band2 tb_threshold_band3"

c1 = 'Flags indicate data that the user should only use with care. In cases of doubt, please refer to the contact person. '
c2 = 'A Fillvalue of 0 means that data has not been flagged. ' 
c3 = 'Bands refer to the measurement ranges (if applicable) of the microwave radiometer; i.e band 1: 20-30 GHz, band 2: 50-60 GHz, band 3: 90 GHz; '
c4 = 'tb valid range: ['+ STRING(threshold(0), format = '(f6.2)') + ', ' + STRING(threshold(1), format = '(f6.2)')+'] in K; ' 
 
NCDF_ATTPUT, outfid, flaid, 'comment', c1+c2+c3+c4

;****define attributes for variables

;........ global attributes
NCDF_ATTPUT, outfid, /GLOBAL, 'Title', 'Microwave radiometer brightness temperatures boundary layer scans' 
NCDF_ATTPUT, outfid, /GLOBAL, 'Institution', ga.i
NCDF_ATTPUT, outfid, /GLOBAL, 'Contact_person', ga.cp
NCDF_ATTPUT, outfid, /GLOBAL, 'Source', ga.s
NCDF_ATTPUT, outfid, /GLOBAL, 'History', ga.h
NCDF_ATTPUT, outfid, /GLOBAL, 'Conventions', ga.cv
NCDF_ATTPUT, outfid, /GLOBAL, 'Processing_date', date_create 
NCDF_ATTPUT, outfid, /GLOBAL, 'Author', ga.au
NCDF_ATTPUT, outfid, /GLOBAL, 'Comments', ga.co
NCDF_ATTPUT, outfid, /GLOBAL, 'License', ga.l
NCDF_ATTPUT, outfid, /GLOBAL, 'Measurement_site', ga.ms

;......... fill all variables with dummy values
NCDF_CONTROL, outfid, /fill

;..........close define mode
NCDF_CONTROL, outfid, /verbose
NCDF_CONTROL, outfid, /endef

;..........assign data to variables
NCDF_VARPUT, outfid, timid, time
NCDF_VARPUT, outfid, tbnid, time_bnds
NCDF_VARPUT, outfid, latid, lat
NCDF_VARPUT, outfid, lonid, lon
NCDF_VARPUT, outfid, altid, zsl
NCDF_VARPUT, outfid, freid, freq
NCDF_VARPUT, outfid, aziid, az
NCDF_VARPUT, outfid, eleid, el
NCDF_VARPUT, outfid, tbrid, tb
NCDF_VARPUT, outfid, tbocrid, offset_tb
NCDF_VARPUT, outfid, fsrid, freq_shift
NCDF_VARPUT, outfid, tbiid, tb_bias
NCDF_VARPUT, outfid, tbcid, tb_cov
NCDF_VARPUT, outfid, temid, temp
NCDF_VARPUT, outfid, preid, pres*100.
NCDF_VARPUT, outfid, rhuid, relh
NCDF_VARPUT, outfid, flaid, flag

;..........close netcdf file
NCDF_CLOSE, outfid
END
