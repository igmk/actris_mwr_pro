;+
;******************
PRO WRITE_LEVEL2B,$
;******************
;INPUT:
out_path,$                      ;path where .nc file will be written
date,$                          ;YYMMDD (STRING)
measurement,$                   ;three letter code identifying the measurement (STRING)
lat, lon, zsl,$                 ;three FLOAT
time,$                          ;array of decimal hours
time_bnds,$                     ;time array bounds
flag_T2b,$                      ;data flagging T, see attributes for details
flag_q2b,$                      ;data flagging q, see attributes for details 
temp,$                          ;surface temperature in K
pres,$                          ;surface pressure in hPa
relh,$                          ;surface relative humidity in %
az,$                            ;array of azimuth angles (if available)
el,$                            ;array of elevation angles
threshold,$                     ;TB & level2a/b thresholds
fn_l1b,$                        ;file naming parameters l1b, hdcp2 compliant
fn_tze,$                        ;file naming parameters tze, hdcp2 compliant
fn_hze,$                        ;file naming parameters hze, hdcp2 compliant
ga,$                            ;global attributes
rets,$                          ;retrieval algorithms used
z_final,$                       ;height grid in m
;KEYWORDS:
T2b=T2b,$                       ;array of temperature 2b profile
dT2b=dT2b,$                     ;temp. offset after TB offset correction T_org-T_oc
err_T2b=err_T2b,$               ;error of temperature 2b
ang_ret_T2b=ang_ret_T2b,$       ;retrieval angles T
q2b=q2b,$                       ;array of abs. humidity 2b profile
dq2b=dq2b,$                     ;hum. offset after TB offset correction q_org-q_oc
err_q2b=err_q2b,$               ;error of abs. humidity 2b
ang_ret_q2b=ang_ret_q2b,$       ;retrieval angles q
dep_l2b_t=dep_l2b_t,$           ;dependency on level1 data version
dep_l2b_q=dep_l2b_q             ;dependency on level1 data version

; $Id: $
; Abstract:
; * write HATPRO level2b netcdf files
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
n_z = N_ELEMENTS(z_final)

;****create netcdf output files

;****write hze file if algorithm specified
IF rets(5) NE '-' THEN BEGIN

 n_ret = N_ELEMENTS(ang_ret_q2b)
 ang_ret = ang_ret_q2b
 outfile_hze = out_path + fn_hze.kkk + '_' + fn_hze.sss + '_' + fn_hze.inst + fn_hze.instnn + '_'+ fn_hze.lll2 + '_' + fn_hze.var + '_' + fn_hze.vnn + '_' + date_string + '.nc'
 dep_l1 = fn_l1b.kkk + '_' + fn_l1b.sss + '_' + fn_l1b.inst + fn_l1b.instnn + '_'+ fn_l1b.lll1 + '_' + fn_l1b.var + '_' + dep_l2b_q
 outfid_hze = NCDF_CREATE(outfile_hze, /clobber)

;*****define dimensions
 timdid = NCDF_DIMDEF(outfid_hze,'time', n_time)
 hgdid = NCDF_DIMDEF(outfid_hze,'height', n_z)
 nvdid = NCDF_DIMDEF(outfid_hze, 'nv', 2)
 rfdid = NCDF_DIMDEF(outfid_hze, 'n_ret', n_ret) 

;*****define variables & attributes
 timid = NCDF_VARDEF(outfid_hze, 'time',[timdid], /double)
 NCDF_ATTPUT, outfid_hze, timid, 'units', 'seconds since 1970-01-01 00:00:00 UTC'
 NCDF_ATTPUT, outfid_hze, timid, 'standard_name', 'time'
 NCDF_ATTPUT, outfid_hze, timid, 'bounds', 'time_bnds'

 tbnid = NCDF_VARDEF(outfid_hze, 'time_bnds',[nvdid, timdid], /double)

 latid = NCDF_VARDEF(outfid_hze, 'lat', /float)
 NCDF_ATTPUT, outfid_hze, latid, 'units','degree_north'
 NCDF_ATTPUT, outfid_hze, latid, 'standard_name','latitude'

 lonid = NCDF_VARDEF(outfid_hze, 'lon', /float)
 NCDF_ATTPUT, outfid_hze, lonid, 'units','degree_east'
 NCDF_ATTPUT, outfid_hze, lonid, 'standard_name','longitude'

 zslid = NCDF_VARDEF(outfid_hze, 'zsl', /float)
 NCDF_ATTPUT, outfid_hze, zslid, 'units', 'm'
 NCDF_ATTPUT, outfid_hze, zslid, 'standard_name','altitude'
 NCDF_ATTPUT, outfid_hze, zslid, 'long_name','altitude above mean sea level'

 aziid = NCDF_VARDEF(outfid_hze, 'azi',[timdid], /float)
 NCDF_ATTPUT, outfid_hze, aziid, 'units', 'degree'
 NCDF_ATTPUT, outfid_hze, aziid, 'standard_name','sensor_azimuth_angle'
 NCDF_ATTPUT, outfid_hze, aziid, '_FillValue',-999.
 NCDF_ATTPUT, outfid_hze, aziid, 'comment', '0=North, 90=East, 180=South, 270=West'
 
 eleid = NCDF_VARDEF(outfid_hze, 'ele',[timdid], /float)
 NCDF_ATTPUT, outfid_hze, eleid, 'units', 'degree'
 NCDF_ATTPUT, outfid_hze, eleid, 'long_name','sensor elevation angle'
 NCDF_ATTPUT, outfid_hze, eleid, '_FillValue',-999.

 releid = NCDF_VARDEF(outfid_hze, 'ele_ret', [rfdid], /float)
 NCDF_ATTPUT, outfid_hze, releid, 'units', 'degree'
 NCDF_ATTPUT, outfid_hze, releid, 'long_name','retrieval elevation angle'
 NCDF_ATTPUT, outfid_hze, releid, '_FillValue',-999.
 NCDF_ATTPUT, outfid_hze, releid, 'comment', 'This variable specifies the elevation angle at which retrievals have been derived.'

 hgid = NCDF_VARDEF(outfid_hze, 'height', [hgdid], /float)
 NCDF_ATTPUT, outfid_hze, hgid, 'units', 'm'
 NCDF_ATTPUT, outfid_hze, hgid, 'standard_name','height'
 NCDF_ATTPUT, outfid_hze, hgid, '_FillValue',-999.

 hzeid = NCDF_VARDEF(outfid_hze, 'hua',[hgdid, timdid], /float)
 NCDF_ATTPUT, outfid_hze, hzeid, 'units', 'kg m-3'
 NCDF_ATTPUT, outfid_hze, hzeid, 'long_name','absolute humidity'
 NCDF_ATTPUT, outfid_hze, hzeid, '_FillValue',-999.
 NCDF_ATTPUT, outfid_hze, hzeid, 'comment', 'hua profiles are given for arbitrary viewing directions in elevation and azimuth.'

 hzeoffid = NCDF_VARDEF(outfid_hze, 'hua_offset',[hgdid, timdid], /float)
 NCDF_ATTPUT, outfid_hze, hzeoffid, 'units', 'kg m-3'
 NCDF_ATTPUT, outfid_hze, hzeoffid, 'long_name','absolute humidity offset correction based on brightness temperature offset'
 NCDF_ATTPUT, outfid_hze, hzeoffid, '_FillValue', 0.
 NCDF_ATTPUT, outfid_hze, hzeoffid, 'comment', 'In order to obtain the un-corrected hua profile, add this offset to hua. This variable is intended for expert use only.'

 hzeeid  = NCDF_VARDEF(outfid_hze, 'hua_err', [rfdid, hgdid], /float)
 NCDF_ATTPUT, outfid_hze, hzeeid, 'units', 'kg m-3'
 NCDF_ATTPUT, outfid_hze, hzeeid, 'long_name','standard error of absolute humidity'
 NCDF_ATTPUT, outfid_hze, hzeeid, '_FillValue',-999.
 NCDF_ATTPUT, outfid_hze, hzeeid, 'comment', 'This variable specifies the uncertainty of hua as a function of height above ground and elevation angle'

 flaid = NCDF_VARDEF(outfid_hze, 'flag',[timdid], /short)
 NCDF_ATTPUT, outfid_hze, flaid, 'long_name','quality control flags'
 NCDF_ATTPUT, outfid_hze, flaid, 'flag_masks', [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]
 NCDF_ATTPUT, outfid_hze, flaid, 'flag_meanings', 'manual_filter_band_1 manual_filter_band2 manual_filter_band3 rain_flag sanity_receiver_band1 sanity_receiver_band2 sun_in_beam tb_threshold_band1 tb_threshold_band2 tb_threshold_band3 iwv_lwp_threshold humidity_threshold'
 NCDF_ATTPUT, outfid_hze, flaid, '_FillValue', 0

 c1 = 'Flags indicate data that the user should only use with care. In cases of doubt, please refer to the contact person. '
 c2 = 'A Fillvalue of 0 means that data has not been flagged. ' 
 c3 = 'Bands refer to the measurement ranges (if applicable) of the microwave radiometer; i.e band 1: 20-30 GHz, band 2: 50-60 GHz, band 3: 90 GHz; '
 c4 = 'tb valid range: ['+ STRING(threshold(0), format = '(f6.2)') + ', ' + STRING(threshold(1), format = '(f6.2)')+'] in K; ' 
 c5 = 'prw valid range: ['+ STRING(threshold(2), format = '(f5.0)') + ', ' + STRING(threshold(3), format = '(f5.0)')+'] in kgm-2; '
 c6 = 'clwvi (zeroing not applied) valid range: ['+ STRING(threshold(4), format = '(f4.1)') + ', ' + STRING(threshold(5), format = '(f4.1)')+'] in kgm-2; '
 c7 = 'hua valid range: ['+ STRING(threshold(8), format = '(f7.4)') + ', ' + STRING(threshold(9), format = '(f7.4)')+'] in kgm-3'
 
 NCDF_ATTPUT, outfid_hze, flaid, 'comment', c1+c2+c3+c4+c5+c6+c7

;****global attributes

 NCDF_ATTPUT, outfid_hze, /GLOBAL, 'Title', 'Microwave radiometer retrieved humidity profile'
 NCDF_ATTPUT, outfid_hze, /GLOBAL, 'Institution', ga.i
 NCDF_ATTPUT, outfid_hze, /GLOBAL, 'Contact_person', ga.cp
 NCDF_ATTPUT, outfid_hze, /GLOBAL, 'Source', ga.s
 NCDF_ATTPUT, outfid_hze, /GLOBAL, 'History', ga.h + 'Retrieval version: ' + rets(5)
 NCDF_ATTPUT, outfid_hze, /GLOBAL, 'Dependencies', dep_l1
 NCDF_ATTPUT, outfid_hze, /GLOBAL, 'Conventions', ga.cv
 NCDF_ATTPUT, outfid_hze, /GLOBAL, 'License', ga.l
 NCDF_ATTPUT, outfid_hze, /GLOBAL, 'Processing_date', date_create
 NCDF_ATTPUT, outfid_hze, /GLOBAL, 'Measurement_site', ga.ms
 NCDF_ATTPUT, outfid_hze, /GLOBAL, 'Author', ga.au

;****fill all variables with dummy values
 NCDF_CONTROL, outfid_hze, /fill

;****close define mode
 NCDF_CONTROL, outfid_hze, /verbose
 NCDF_CONTROL, outfid_hze, /endef

;****assign data to variables and write to file
 NCDF_VARPUT, outfid_hze, timid, time 
 NCDF_VARPUT, outfid_hze, tbnid, time_bnds
 NCDF_VARPUT, outfid_hze, latid, lat
 NCDF_VARPUT, outfid_hze, lonid, lon
 NCDF_VARPUT, outfid_hze, zslid, zsl
 NCDF_VARPUT, outfid_hze, aziid, az
 NCDF_VARPUT, outfid_hze, eleid, el
 NCDF_VARPUT, outfid_hze, releid, ang_ret
 NCDF_VARPUT, outfid_hze, hgid, z_final
 NCDF_VARPUT, outfid_hze, hzeid, q2b
 NCDF_VARPUT, outfid_hze, hzeoffid, dq2b
 NCDF_VARPUT, outfid_hze, hzeeid, err_q2b
 NCDF_VARPUT, outfid_hze, flaid, flag_q2b

;****close netcdf file
 NCDF_CLOSE, outfid_hze

ENDIF ; hze file

;****write tze file if algorithm specified
IF rets(4) NE '-' THEN BEGIN

 n_ret = N_ELEMENTS(ang_ret_T2b)
 ang_ret = ang_ret_T2b

 outfile_tze = out_path + fn_tze.kkk + '_' + fn_tze.sss + '_' + fn_tze.inst + fn_tze.instnn + '_'+ fn_tze.lll2 + '_' + fn_tze.var + '_' + fn_tze.vnn + '_' + date_string + '.nc'
 dep_l1 = fn_l1b.kkk + '_' + fn_l1b.sss + '_' + fn_l1b.inst + fn_l1b.instnn + '_'+ fn_l1b.lll1 + '_' + fn_l1b.var + '_' + dep_l2b_t
 outfid_tze = NCDF_CREATE(outfile_tze, /clobber)

;*****define dimensions
 timdid = NCDF_DIMDEF(outfid_tze,'time', n_time)
 hgdid = NCDF_DIMDEF(outfid_tze, 'height', n_z)
 nvdid = NCDF_DIMDEF(outfid_tze, 'nv', 2)
 rfdid = NCDF_DIMDEF(outfid_tze, 'n_ret', n_ret) 

;*****define variables & attributes
 timid = NCDF_VARDEF(outfid_tze, 'time',[timdid], /double)
 NCDF_ATTPUT, outfid_tze, timid, 'units', 'seconds since 1970-01-01 00:00:00 UTC'
 NCDF_ATTPUT, outfid_tze, timid, 'standard_name', 'time'
 NCDF_ATTPUT, outfid_tze, timid, 'bounds', 'time_bnds'

 tbnid = NCDF_VARDEF(outfid_tze, 'time_bnds',[nvdid, timdid], /double)

 latid = NCDF_VARDEF(outfid_tze, 'lat', /float)
 NCDF_ATTPUT, outfid_tze, latid, 'units','degree_north'
 NCDF_ATTPUT, outfid_tze, latid, 'standard_name','latitude'

 lonid = NCDF_VARDEF(outfid_tze, 'lon', /float)
 NCDF_ATTPUT, outfid_tze, lonid, 'units','degree_east'
 NCDF_ATTPUT, outfid_tze, lonid, 'standard_name','longitude'

 zslid = NCDF_VARDEF(outfid_tze, 'zsl', /float)
 NCDF_ATTPUT, outfid_tze, zslid, 'units', 'm'
 NCDF_ATTPUT, outfid_tze, zslid, 'standard_name','altitude'
 NCDF_ATTPUT, outfid_hze, zslid, 'long_name','altitude above mean sea level'

 aziid = NCDF_VARDEF(outfid_tze, 'azi',[timdid], /float)
 NCDF_ATTPUT, outfid_tze, aziid, 'units', 'degree'
 NCDF_ATTPUT, outfid_tze, aziid, 'standard_name','sensor_azimuth_angle'
 NCDF_ATTPUT, outfid_tze, aziid, '_FillValue',-999.
 NCDF_ATTPUT, outfid_tze, aziid, 'comment', '0=North, 90=East, 180=South, 270=West'
 
 eleid = NCDF_VARDEF(outfid_tze, 'ele',[timdid], /float)
 NCDF_ATTPUT, outfid_tze, eleid, 'units', 'degree'
 NCDF_ATTPUT, outfid_tze, eleid, 'long_name','sensor elevation angle'
 NCDF_ATTPUT, outfid_tze, eleid, '_FillValue',-999.

 releid = NCDF_VARDEF(outfid_tze, 'ele_ret', [rfdid], /float)
 NCDF_ATTPUT, outfid_tze, releid, 'units', 'degree'
 NCDF_ATTPUT, outfid_tze, releid, 'long_name','retrieval elevation angle'
 NCDF_ATTPUT, outfid_tze, releid, '_FillValue',-999.
 NCDF_ATTPUT, outfid_tze, releid, 'comment', 'This variable specifies the elevation angle at which retrievals have been derived.'

 hgid = NCDF_VARDEF(outfid_hze, 'height', [hgdid], /float)
 NCDF_ATTPUT, outfid_hze, hgid, 'units', 'm'
 NCDF_ATTPUT, outfid_hze, hgid, 'standard_name','height'
 NCDF_ATTPUT, outfid_hze, hgid, '_FillValue',-999.

 tzeid = NCDF_VARDEF(outfid_tze, 'ta',[hgdid, timdid], /float)
 NCDF_ATTPUT, outfid_tze, tzeid, 'units', 'K'
 NCDF_ATTPUT, outfid_tze, tzeid, 'standard_name','air_temperature'
 NCDF_ATTPUT, outfid_tze, tzeid, '_FillValue',-999.
 NCDF_ATTPUT, outfid_tze, tzeid, 'comment', 'ta profiles are given for arbitrary viewing directions in elevation and azimuth.'

 tzeoffid = NCDF_VARDEF(outfid_tze, 'ta_offset',[hgdid, timdid], /float)
 NCDF_ATTPUT, outfid_tze, tzeoffid, 'units', 'K'
 NCDF_ATTPUT, outfid_tze, tzeoffid, 'long_name','air_temperature offset correction based on brightness temperature offset'
 NCDF_ATTPUT, outfid_tze, tzeoffid, '_FillValue', 0.
 NCDF_ATTPUT, outfid_tze, tzeoffid, 'comment', 'In order to obtain the un-corrected ta profile, add this offset to ta. This variable is intended for expert use only.'

 tzeeid  = NCDF_VARDEF(outfid_tze, 'ta_err', [rfdid, hgdid], /float)
 NCDF_ATTPUT, outfid_tze, tzeeid, 'units', 'K'
 NCDF_ATTPUT, outfid_tze, tzeeid, 'long_name','standard error of temperature'
 NCDF_ATTPUT, outfid_tze, tzeeid, '_FillValue',-999.
 NCDF_ATTPUT, outfid_tze, tzeeid, 'comment', 'This variable specifies the uncertainty of ta as a function of height above ground and elevation angle'
 
 flaid = NCDF_VARDEF(outfid_tze, 'flag',[timdid], /short)
 NCDF_ATTPUT, outfid_tze, flaid, 'long_name','quality control flags'
 NCDF_ATTPUT, outfid_tze, flaid, 'flag_masks', [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048]
 NCDF_ATTPUT, outfid_tze, flaid, 'flag_meanings', 'manual_filter_band_1 manual_filter_band2 manual_filter_band3 rain_flag sanity_receiver_band1 sanity_receiver_band2 sun_in_beam tb_threshold_band1 tb_threshold_band2 tb_threshold_band3 iwv_lwp_threshold temperature_threshold'
 NCDF_ATTPUT, outfid_tze, flaid, '_FillValue', 0

 c1 = 'Flags indicate data that the user should only use with care. In cases of doubt, please refer to the contact person. '
 c2 = 'A Fillvalue of 0 means that data has not been flagged. ' 
 c3 = 'Bands refer to the measurement ranges (if applicable) of the microwave radiometer; i.e band 1: 20-30 GHz, band 2: 50-60 GHz, band 3: 90 GHz; '
 c4 = 'tb valid range: ['+ STRING(threshold(0), format = '(f6.2)') + ', ' + STRING(threshold(1), format = '(f6.2)')+'] in K; ' 
 c5 = 'prw valid range: ['+ STRING(threshold(2), format = '(f5.0)') + ', ' + STRING(threshold(3), format = '(f5.0)')+'] in kgm-2; '
 c6 = 'clwvi (zeroing not applied) valid range: ['+ STRING(threshold(4), format = '(f4.1)') + ', ' + STRING(threshold(5), format = '(f4.1)')+'] in kgm-2; '
 c7 = 'ta valid range: ['+ STRING(threshold(6), format = '(f6.2)') + ', ' + STRING(threshold(7), format = '(f6.2)')+'] in K'
 
 NCDF_ATTPUT, outfid_tze, flaid, 'comment', c1+c2+c3+c4+c5+c6+c7

;****global attributes

 NCDF_ATTPUT, outfid_tze, /GLOBAL, 'Title', 'Microwave radiometer retrieved temeprature profile'
 NCDF_ATTPUT, outfid_tze, /GLOBAL, 'Institution', ga.i
 NCDF_ATTPUT, outfid_tze, /GLOBAL, 'Contact_person', ga.cp
 NCDF_ATTPUT, outfid_tze, /GLOBAL, 'Source', ga.s
 NCDF_ATTPUT, outfid_tze, /GLOBAL, 'History', ga.h + 'Retrieval version: ' + rets(4)
 NCDF_ATTPUT, outfid_tze, /GLOBAL, 'Dependencies', dep_l1
 NCDF_ATTPUT, outfid_tze, /GLOBAL, 'Conventions', ga.cv
 NCDF_ATTPUT, outfid_tze, /GLOBAL, 'License', ga.l
 NCDF_ATTPUT, outfid_tze, /GLOBAL, 'Processing_date', date_create
 NCDF_ATTPUT, outfid_tze, /GLOBAL, 'Measurement_site', ga.ms
 NCDF_ATTPUT, outfid_tze, /GLOBAL, 'Author', ga.au

;****fill all variables with dummy values
 NCDF_CONTROL, outfid_tze, /fill

;****close define mode
 NCDF_CONTROL, outfid_tze, /verbose
 NCDF_CONTROL, outfid_tze, /endef

;****assign data to variables and write to file
 NCDF_VARPUT, outfid_tze, timid, time 
 NCDF_VARPUT, outfid_tze, tbnid, time_bnds
 NCDF_VARPUT, outfid_tze, latid, lat
 NCDF_VARPUT, outfid_tze, lonid, lon
 NCDF_VARPUT, outfid_tze, zslid, zsl
 NCDF_VARPUT, outfid_tze, aziid, az
 NCDF_VARPUT, outfid_tze, eleid, el
 NCDF_VARPUT, outfid_tze, releid, ang_ret
 NCDF_VARPUT, outfid_tze, hgid, z_final 
 NCDF_VARPUT, outfid_tze, tzeid, T2b
 NCDF_VARPUT, outfid_tze, tzeoffid, dT2b
 NCDF_VARPUT, outfid_tze, tzeeid, err_T2b
 NCDF_VARPUT, outfid_tze, flaid, flag_T2b

;****close netcdf file
 NCDF_CLOSE, outfid_tze

ENDIF ; tze file

END

