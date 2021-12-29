;+
;******************
PRO WRITE_LEVEL2C,$
;******************
;INPUT:
out_path,$                      ;path where .nc file will be written
date,$                          ;YYMMDD (STRING)
measurement,$                   ;three letter code identifying the measurement (STRING)
lat, lon, zsl,$                 ;three FLOAT
time,$                          ;array of decimal hours
time_bnds,$                     ;time array bounds
flag,$                          ;data flagging, see attributes for details
temp,$                          ;surface temperature in K
pres,$                          ;surface pressure in hPa
relh,$                          ;surface relative humidity in %
az,$                            ;array of azimuth angles (if available)
el,$                            ;array of elevation angles
threshold,$                     ;TB & level2a/b/c thresholds
fn_l1c,$                         ;file naming parameters l1c, hdcp2 compliant
fn_tel,$                        ;file naming parameters l2c, hdcp2 compliant
ga,$                            ;global attributes
rets,$                          ;retrieval algorithm used
z_final,$                       ;height grid in m
;KEYWORDS:
T2c=T2c,$                       ;array of temperature 2c profile
dT2c=dT2c,$                     ;array of temp. offsets (T_org-T_oc)
err_T2c=err_T2c,$               ;error of temperature 2c
dep_l2c_t=dep_l2c_t             ;dependency on level1 data version

; $Id: $
; Abstract:
; * write HATPRO level2a netcdf files
; Authors:
; 
; Date:
; 2011-05-25
; Dependencies:
; -
; Changes:
; 2013-11-11, 2014-02-07, 2014-11-06:
; UL: make written level2a file HDCP(2) compliant

;-

;-->the option for writing taa as a function of frequency still needs to be included!


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

time_conv = time
time_bnds_conv = time_bnds

;****transform time to seconds since 19700101
SECONDS_SINCE_19700101, year+month+day, time_conv, time

;****transform time_bnds to seconds since 19700101
SECONDS_SINCE_19700101, year+month+day, time_bnds_conv, time_bnds
time_bnds = TRANSPOSE(time_bnds)

;****initiate dimensions
n_time = N_ELEMENTS(time)
n_freq_taa = N_ELEMENTS(f_taa)
n_ang = N_ELEMENTS(el)
n_z = N_ELEMENTS(z_final)

;****create netcdf output files

;****write tel file if algorithm specified
IF rets(6) NE '-' THEN BEGIN

 outfile_tel = out_path + fn_tel.kkk + '_' + fn_tel.sss + '_' + fn_tel.inst + fn_tel.instnn + '_'+ fn_tel.lll2 + '_' + fn_tel.var + '_' + fn_tel.vnn + '_' + date_string + '.nc'
 dep_l1 = fn_l1c.kkk + '_' + fn_l1c.sss + '_' + fn_l1c.inst + fn_l1c.instnn + '_'+ fn_l1c.lll1 + '_' + fn_l1c.var + '_' + dep_l2c_t
 outfid_tel = NCDF_CREATE(outfile_tel, /clobber)

;*****define dimensions
 timdid = NCDF_DIMDEF(outfid_tel,'time', n_time)
 hgdid = NCDF_DIMDEF(outfid_tel,'height', n_z)
 nvdid = NCDF_DIMDEF(outfid_tel, 'nv', 2)
 angdid = NCDF_DIMDEF(outfid_tel, 'n_ang', n_ang) 

;*****define variables & attributes
 timid = NCDF_VARDEF(outfid_tel, 'time',[timdid], /double)
 NCDF_ATTPUT, outfid_tel, timid, 'units', 'seconds since 1970-01-01 00:00:00 UTC'
 NCDF_ATTPUT, outfid_tel, timid, 'standard_name', 'time'
 NCDF_ATTPUT, outfid_tel, timid, 'bounds', 'time_bnds'

 tbnid = NCDF_VARDEF(outfid_tel, 'time_bnds',[nvdid, timdid], /double)

 latid = NCDF_VARDEF(outfid_tel, 'lat', /float)
 NCDF_ATTPUT, outfid_tel, latid, 'units','degree_north'
 NCDF_ATTPUT, outfid_tel, latid, 'standard_name','latitude'

 lonid = NCDF_VARDEF(outfid_tel, 'lon', /float)
 NCDF_ATTPUT, outfid_tel, lonid, 'units','degree_east'
 NCDF_ATTPUT, outfid_tel, lonid, 'standard_name','longitude'

 zslid = NCDF_VARDEF(outfid_tel, 'zsl', /float)
 NCDF_ATTPUT, outfid_tel, zslid, 'units', 'm'
 NCDF_ATTPUT, outfid_tel, zslid, 'standard_name','altitude'

 aziid = NCDF_VARDEF(outfid_tel, 'azi', /float)
 NCDF_ATTPUT, outfid_tel, aziid, 'units', 'degree'
 NCDF_ATTPUT, outfid_tel, aziid, 'standard_name','sensor_azimuth_angle'
 NCDF_ATTPUT, outfid_tel, aziid, '_FillValue',-999.
 NCDF_ATTPUT, outfid_tel, aziid, 'comment', '0=North, 90=East, 180=South, 270=West'
 
 eleid = NCDF_VARDEF(outfid_tel, 'ele', [angdid], /float)
 NCDF_ATTPUT, outfid_tel, eleid, 'units', 'degree'
 NCDF_ATTPUT, outfid_tel, eleid, 'long_name','sensor elevation angle'
 NCDF_ATTPUT, outfid_tel, eleid, '_FillValue',-999.

 hgid = NCDF_VARDEF(outfid_tel, 'height', [hgdid], /float)
 NCDF_ATTPUT, outfid_tel, hgid, 'units', 'm'
 NCDF_ATTPUT, outfid_tel, hgid, 'standard_name','height'
 NCDF_ATTPUT, outfid_tel, hgid, '_FillValue',-999.

 telid = NCDF_VARDEF(outfid_tel, 'ta',[hgdid, timdid], /float)
 NCDF_ATTPUT, outfid_tel, telid, 'units', 'K'
 NCDF_ATTPUT, outfid_tel, telid, 'standard_name','air_temperature'
 NCDF_ATTPUT, outfid_tel, telid, '_FillValue',-999.
 NCDF_ATTPUT, outfid_tel, telid, 'comment', 'ta profiles are vertical profiles over the measurement site. Accurcay is especially high in the lowest 1km.' 

 teloffid = NCDF_VARDEF(outfid_tel, 'ta_offset',[hgdid, timdid], /float)
 NCDF_ATTPUT, outfid_tel, teloffid, 'units', 'K'
 NCDF_ATTPUT, outfid_tel, teloffid, 'long_name','air_temperature offset correction based on brightness temperature offset'
 NCDF_ATTPUT, outfid_tel, teloffid, '_FillValue', 0.
 NCDF_ATTPUT, outfid_tel, teloffid, 'comment', 'In order to obtain the un-corrected ta profile, add this offset to ta. This variable is intended for expert use only.'

 teleid  = NCDF_VARDEF(outfid_tel, 'ta_err', [hgdid], /float)
 NCDF_ATTPUT, outfid_tel, teleid, 'units', 'K'
 NCDF_ATTPUT, outfid_tel, teleid, 'long_name','standard error of air_temperature'
 NCDF_ATTPUT, outfid_tel, teleid, '_FillValue',-999.

 flaid = NCDF_VARDEF(outfid_tel, 'flag',[timdid], /short)
 NCDF_ATTPUT, outfid_tel, flaid, 'long_name','quality control flags'
 NCDF_ATTPUT, outfid_tel, flaid, 'flag_masks', [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]
 NCDF_ATTPUT, outfid_tel, flaid, 'flag_meanings', "manual_filter_band_1 manual_filter_band2 manual_filter_band3 rain_flag sanity_receiver_band1 sanity_receiver_band2 sun_in_beam tb_threshold_band1 tb_threshold_band2 tb_threshold_band3 iwv_lwp_threshold temperature_threshold"
 NCDF_ATTPUT, outfid_tel, flaid, '_FillValue', 0

 c1 = 'Flags indicate data that the user should only use with care. In cases of doubt, please refer to the contact person. '
 c2 = 'A Fillvalue of 0 means that data has not been flagged. ' 
 c3 = 'Bands refer to the measurement ranges (if applicable) of the microwave radiometer; i.e band 1: 20-30 GHz, band 2: 50-60 GHz, band 3: 90 GHz; '
 c4 = 'tb valid range: ['+ STRING(threshold(0), format = '(f6.2)') + ', ' + STRING(threshold(1), format = '(f6.2)')+'] in K; ' 
 c5 = 'prw valid range: ['+ STRING(threshold(2), format = '(f5.0)') + ', ' + STRING(threshold(3), format = '(f5.0)')+'] in kgm-2; '
 c6 = 'clwvi (zeroing not applied) valid range: ['+ STRING(threshold(4), format = '(f4.1)') + ', ' + STRING(threshold(5), format = '(f4.1)')+'] in kgm-2; '
 c7 = 'ta valid range: ['+ STRING(threshold(6), format = '(f6.2)') + ', ' + STRING(threshold(7), format = '(f6.2)')+'] in K'     
 
 NCDF_ATTPUT, outfid_tel, flaid, 'comment', c1+c2+c3+c4+c5+c6+c7

;****global attributes

 NCDF_ATTPUT, outfid_tel, /GLOBAL, 'Title', 'Microwave radiometer retrieved temperature profile from boundary layers scans'
 NCDF_ATTPUT, outfid_tel, /GLOBAL, 'Institution', ga.i
 NCDF_ATTPUT, outfid_tel, /GLOBAL, 'Contact_person', ga.cp
 NCDF_ATTPUT, outfid_tel, /GLOBAL, 'Source', ga.s
 NCDF_ATTPUT, outfid_tel, /GLOBAL, 'History', ga.h + 'Retrieval version: ' + rets(6)
 NCDF_ATTPUT, outfid_tel, /GLOBAL, 'Dependencies', dep_l1
 NCDF_ATTPUT, outfid_tel, /GLOBAL, 'Conventions', ga.cv
 NCDF_ATTPUT, outfid_tel, /GLOBAL, 'License', ga.l
 NCDF_ATTPUT, outfid_tel, /GLOBAL, 'Processing_date', date_create
 NCDF_ATTPUT, outfid_tel, /GLOBAL, 'Measurement_site', ga.ms
 NCDF_ATTPUT, outfid_tel, /GLOBAL, 'Author', ga.au
 NCDF_ATTPUT, outfid_tel, /GLOBAL, 'Comment', ga.co

;****fill all variables with dummy values
 NCDF_CONTROL, outfid_tel, /fill

;****close define mode
 NCDF_CONTROL, outfid_tel, /verbose
 NCDF_CONTROL, outfid_tel, /endef

;****assign data to variables and write to file
 NCDF_VARPUT, outfid_tel, timid, time 
 NCDF_VARPUT, outfid_tel, tbnid, time_bnds
 NCDF_VARPUT, outfid_tel, latid, lat
 NCDF_VARPUT, outfid_tel, lonid, lon
 NCDF_VARPUT, outfid_tel, zslid, zsl
 NCDF_VARPUT, outfid_tel, aziid, az
 NCDF_VARPUT, outfid_tel, eleid, el
 NCDF_VARPUT, outfid_tel, hgid, z_final
 NCDF_VARPUT, outfid_tel, telid, T2c
 NCDF_VARPUT, outfid_tel, teloffid, dT2c
 NCDF_VARPUT, outfid_tel, teleid, err_T2c
 NCDF_VARPUT, outfid_tel, flaid, flag

;****close netcdf file
 NCDF_CLOSE, outfid_tel

ENDIF ; tel file

END

