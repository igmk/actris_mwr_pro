;+
;******************
PRO WRITE_LEVEL2A,$
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
threshold,$                     ;TB & level2a thresholds
fn_l1b,$                        ;file naming parameters for level1b, hdcp2 compliant
fn_iwv,$                        ;file naming parameters for iwv, hdcp2 compliant
fn_lwp,$                        ;file naming parameters for lwp, hdcp2 compliant
fn_wdl,$                        ;file naming parameters for wdl, hdcp2 compliant
fn_taa,$                        ;file naming parameters for taa, hdcp2 compliant
ga,$                            ;global attributes
rets,$                          ;retrieval algorithms used
ang_ret,$                       ;retrieval elevation angles
;KEYWORDS:
iwv=iwv,$                       ;array of Integrated Water Vapor
diwv=diwv,$                     ;array of iwv offset; iwv_org-iwv_oc
err_iwv=err_iwv,$               ;error of iwv 
lwpo=lwpo,$                     ;array of Liquid Water Path
lwp_cor=lwp_cor,$               ;array of offset corrected Liquid Water Path (zeroing)
dlwp=dlwp,$                     ;array of offset lwp offset (TB corrected); lwp_org-lwp_oc
err_lwp=err_lwp,$               ;error of lwp 
wdl=wdl,$                       ;array of Zenith Wet Dely
dwdl=dwdl,$                     ;array of wdl offset; wdl_org-wdl_oc
err_wdl=err_wdl,$               ;error of wdl
taa_a=taa_a,$                   ;array of Atmospheric Attenuation
dtaa_a=dtaa_a,$                 ;array of taa offset; taa_org-taa_oc
err_taa=err_taa,$               ;error of taa_a
f_taa=f_taa,$                   ;frequecies of total atmospheric attenuation
dep_l2a_iwv=dep_l2a_iwv,$       ;dependency on level1 data version
dep_l2a_lwp=dep_l2a_lwp         ;dependency on level1 data version

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
n_ret = N_ELEMENTS(ang_ret)

;****differentiate between zenith and off-zenith variables
iwv_zen = REPLICATE(-999., n_time)
iwv_off_zen = REPLICATE(-999., n_time)
diwv_zen = REPLICATE(0., n_time)
diwv_off_zen = REPLICATE(0., n_time) 
lwp_zen = REPLICATE(-999., n_time)
lwp_off_zen = REPLICATE(-999., n_time)
lwp_offset = REPLICATE(-999., n_time)
dlwp_zen = REPLICATE(0., n_time)
dlwp_off_zen = REPLICATE(0., n_time) 

i_off_zen = WHERE(el LE 89.4 OR el GE 90.6)
i_zen = WHERE(el GT 89.4 AND el LT 90.6)

IF i_zen(0) NE -1 THEN BEGIN 
 iwv_zen(i_zen) = iwv(i_zen)
 lwp_zen(i_zen) = lwp_cor(i_zen)
 diwv_zen(i_zen) = diwv(i_zen)
 dlwp_zen(i_zen) = dlwp(i_zen)

;**lwp_offset
 i_zen_off = WHERE(el GT 89.4 AND el LT 90.6 AND lwpo NE -999. AND lwp_cor NE -999.)
 IF i_zen_off(0) NE -1 THEN lwp_offset(i_zen_off) = lwpo(i_zen_off) - lwp_cor(i_zen_off)  
ENDIF 

IF i_off_zen(0) NE -1 THEN BEGIN 
 iwv_off_zen(i_off_zen) = iwv(i_off_zen)
 lwp_off_zen(i_off_zen) = lwpo(i_off_zen)
 diwv_off_zen(i_off_zen) = diwv(i_off_zen)
 dlwp_off_zen(i_off_zen) = dlwp(i_off_zen)
ENDIF 

;****create netcdf output files

;****write IWV file if algorithm specified
IF rets(0) NE '-' THEN BEGIN

 outfile_iwv = out_path + fn_iwv.kkk + '_' + fn_iwv.sss + '_' + fn_iwv.inst + fn_iwv.instnn + '_'+ fn_iwv.lll2 + '_' + fn_iwv.var + '_' + fn_iwv.vnn + '_' + date_string + '.nc'
 dep_l1 = fn_l1b.kkk + '_' + fn_l1b.sss + '_' + fn_l1b.inst + fn_l1b.instnn + '_'+ fn_l1b.lll1 + '_' + fn_l1b.var + '_' + dep_l2a_iwv
 outfid_iwv = NCDF_CREATE(outfile_iwv, /clobber)

;*****define dimensions
 timdid = NCDF_DIMDEF(outfid_iwv,'time', n_time)
 nvdid = NCDF_DIMDEF(outfid_iwv, 'nv', 2)
 rfdid = NCDF_DIMDEF(outfid_iwv, 'n_ret', n_ret)

;*****define variables & attributes
 timid = NCDF_VARDEF(outfid_iwv, 'time',[timdid], /double)
 NCDF_ATTPUT, outfid_iwv, timid, 'units', 'seconds since 1970-01-01 00:00:00 UTC'
 NCDF_ATTPUT, outfid_iwv, timid, 'standard_name', 'time'
 NCDF_ATTPUT, outfid_iwv, timid, 'bounds', 'time_bnds'

 tbnid = NCDF_VARDEF(outfid_iwv, 'time_bnds',[nvdid, timdid], /double)

 latid = NCDF_VARDEF(outfid_iwv, 'lat', /float)
 NCDF_ATTPUT, outfid_iwv, latid, 'units','degree_north'
 NCDF_ATTPUT, outfid_iwv, latid, 'standard_name','latitude'

 lonid = NCDF_VARDEF(outfid_iwv, 'lon', /float)
 NCDF_ATTPUT, outfid_iwv, lonid, 'units','degree_east'
 NCDF_ATTPUT, outfid_iwv, lonid, 'standard_name','longitude'

 zslid = NCDF_VARDEF(outfid_iwv, 'zsl', /float)
 NCDF_ATTPUT, outfid_iwv, zslid, 'units', 'm'
 NCDF_ATTPUT, outfid_iwv, zslid, 'standard_name','altitude'
 NCDF_ATTPUT, outfid_iwv, zslid, 'long_name','altitude above mean sea level'

 aziid = NCDF_VARDEF(outfid_iwv, 'azi',[timdid], /float)
 NCDF_ATTPUT, outfid_iwv, aziid, 'units', 'degree'
 NCDF_ATTPUT, outfid_iwv, aziid, 'standard_name','sensor_azimuth_angle'
 NCDF_ATTPUT, outfid_iwv, aziid, '_FillValue',-999.
 NCDF_ATTPUT, outfid_iwv, aziid, 'comment', '0=North, 90=East, 180=South, 270=West'
 
 eleid = NCDF_VARDEF(outfid_iwv, 'ele',[timdid], /float)
 NCDF_ATTPUT, outfid_iwv, eleid, 'units', 'degree'
 NCDF_ATTPUT, outfid_iwv, eleid, 'long_name','sensor elevation angle'
 NCDF_ATTPUT, outfid_iwv, eleid, '_FillValue',-999.

 releid = NCDF_VARDEF(outfid_iwv, 'ele_ret',[rfdid], /float)
 NCDF_ATTPUT, outfid_iwv, releid, 'units', 'degree'
 NCDF_ATTPUT, outfid_iwv, releid, 'long_name','retrieval elevation angle'
 NCDF_ATTPUT, outfid_iwv, releid, '_FillValue',-999.

 iwvid = NCDF_VARDEF(outfid_iwv, 'prw',[timdid], /float)
 NCDF_ATTPUT, outfid_iwv, iwvid, 'units', 'kg m-2'
 NCDF_ATTPUT, outfid_iwv, iwvid, 'standard_name','atmosphere_mass_content_of_water_vapor'
 NCDF_ATTPUT, outfid_iwv, iwvid, '_FillValue',-999.
 NCDF_ATTPUT, outfid_iwv, iwvid, 'comment', 'These values denote the vertically integrated amount of water vapor from the surface to TOA.'

 iwvoffid  = NCDF_VARDEF(outfid_iwv, 'prw_offset',[timdid], /float)
 NCDF_ATTPUT, outfid_iwv, iwvoffid, 'units', 'kg m-2'
 NCDF_ATTPUT, outfid_iwv, iwvoffid, 'long_name','atmosphere_mass_content_of_water_vapor offset correction based on brightness temperature offset'
 NCDF_ATTPUT, outfid_iwv, iwvoffid, '_FillValue', -999.
 NCDF_ATTPUT, outfid_iwv, iwvoffid, 'comment', 'This value has been subtracted from the original prw value to account for instrument calibration drifts. The information is designated for expert user use.'

 iwvozid = NCDF_VARDEF(outfid_iwv, 'prw_off_zenith',[timdid], /float)
 NCDF_ATTPUT, outfid_iwv, iwvozid, 'units', 'kg m-2'
 NCDF_ATTPUT, outfid_iwv, iwvozid, 'long_name','off zenith path integrated water vapor'
 NCDF_ATTPUT, outfid_iwv, iwvozid, '_FillValue', -999.
 NCDF_ATTPUT, outfid_iwv, iwvozid, 'comment', 'These values denote the path integrated amount of water vapor at arbitrary viewing directions in elevation and azimuth.'

 iwvozoffid  = NCDF_VARDEF(outfid_iwv, 'prw_off_zenith_offset',[timdid], /float)
 NCDF_ATTPUT, outfid_iwv, iwvozoffid, 'units', 'kg m-2'
 NCDF_ATTPUT, outfid_iwv, iwvozoffid, 'long_name','off zenith path integrated water vapor offset correction based on brightness temperature offset'
 NCDF_ATTPUT, outfid_iwv, iwvozoffid, '_FillValue', -999.
 NCDF_ATTPUT, outfid_iwv, iwvozoffid, 'comment', 'This value has been subtracted from the original iwv value to account for instrument calibration drifts. The information is designated for expert user use.'

 iwveid  = NCDF_VARDEF(outfid_iwv, 'prw_err',[rfdid], /float)
 NCDF_ATTPUT, outfid_iwv, iwveid, 'units', 'kg m-2'
 NCDF_ATTPUT, outfid_iwv, iwveid, 'long_name','These values denote the standard error of atmosphere mass content of water vapor or off zenith path integrated water vapor at the corresponding retrieval elevation angle.'
 NCDF_ATTPUT, outfid_iwv, iwveid, '_FillValue',-999.

 flaid = NCDF_VARDEF(outfid_iwv, 'flag',[timdid], /short)
 NCDF_ATTPUT, outfid_iwv, flaid, 'long_name','quality control flags'
 NCDF_ATTPUT, outfid_iwv, flaid, 'flag_masks', [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
 NCDF_ATTPUT, outfid_iwv, flaid, 'flag_meanings', "manual_filter_band_1 manual_filter_band2 manual_filter_band3 rain_flag sanity_receiver_band1 sanity_receiver_band2 sun_in_beam tb_threshold_band1 tb_threshold_band2 tb_threshold_band3 iwv_lwp_threshold"
 NCDF_ATTPUT, outfid_iwv, flaid, '_FillValue', 0

 c1 = 'Flags indicate data that the user should only use with care. In cases of doubt, please refer to the contact person. '
 c2 = 'A Fillvalue of 0 means that data has not been flagged. ' 
 c3 = 'Bands refer to the measurement ranges (if applicable) of the microwave radiometer; i.e band 1: 20-30 GHz, band 2: 50-60 GHz, band 3: 90 GHz; '
 c4 = 'tb valid range: ['+ STRING(threshold(0), format = '(f6.2)') + ', ' + STRING(threshold(1), format = '(f6.2)')+'] in K; ' 
 c5 = 'prw valid range: ['+ STRING(threshold(2), format = '(f5.0)') + ', ' + STRING(threshold(3), format = '(f5.0)')+'] in kgm-2; '
 c6 = 'clwvi (zeroing not applied) valid range: ['+ STRING(threshold(4), format = '(f4.1)') + ', ' + STRING(threshold(5), format = '(f4.1)')+'] in kgm-2; '
 
 NCDF_ATTPUT, outfid_iwv, flaid, 'comment', c1+c2+c3+c4+c5+c6

;****global attributes

 NCDF_ATTPUT, outfid_iwv, /GLOBAL, 'Title', 'Microwave radiometer retrieved prw'
 NCDF_ATTPUT, outfid_iwv, /GLOBAL, 'Institution', ga.i
 NCDF_ATTPUT, outfid_iwv, /GLOBAL, 'Contact_person', ga.cp
 NCDF_ATTPUT, outfid_iwv, /GLOBAL, 'Source', ga.s
 NCDF_ATTPUT, outfid_iwv, /GLOBAL, 'History', ga.h + 'Retrieval version: ' + rets(0)
 NCDF_ATTPUT, outfid_iwv, /GLOBAL, 'Dependencies', dep_l1
 NCDF_ATTPUT, outfid_iwv, /GLOBAL, 'Conventions', ga.cv
 NCDF_ATTPUT, outfid_iwv, /GLOBAL, 'License', ga.l
 NCDF_ATTPUT, outfid_iwv, /GLOBAL, 'Processing_date', date_create
 NCDF_ATTPUT, outfid_iwv, /GLOBAL, 'Measurement_site', ga.ms
 NCDF_ATTPUT, outfid_iwv, /GLOBAL, 'Author', ga.au
 NCDF_ATTPUT, outfid_iwv, /GLOBAL, 'Comment', 'prw errors are given as a function of retrieval elevation angle. Note, observations may occur at different angles (interpolation).'

;****fill all variables with dummy values
 NCDF_CONTROL, outfid_iwv, /fill

;****close define mode
 NCDF_CONTROL, outfid_iwv, /verbose
 NCDF_CONTROL, outfid_iwv, /endef

;****assign data to variables and write to file
 NCDF_VARPUT, outfid_iwv, timid, time 
 NCDF_VARPUT, outfid_iwv, tbnid, time_bnds
 NCDF_VARPUT, outfid_iwv, latid, lat
 NCDF_VARPUT, outfid_iwv, lonid, lon
 NCDF_VARPUT, outfid_iwv, zslid, zsl
 NCDF_VARPUT, outfid_iwv, aziid, az
 NCDF_VARPUT, outfid_iwv, eleid, el
 NCDF_VARPUT, outfid_iwv, releid, ang_ret
 NCDF_VARPUT, outfid_iwv, iwvid, iwv_zen
 NCDF_VARPUT, outfid_iwv, iwvoffid, diwv_zen
 NCDF_VARPUT, outfid_iwv, iwvozid, iwv_off_zen
 NCDF_VARPUT, outfid_iwv, iwvozoffid, diwv_off_zen
 NCDF_VARPUT, outfid_iwv, iwveid, err_iwv
 NCDF_VARPUT, outfid_iwv, flaid, flag

;****close netcdf file
 NCDF_CLOSE, outfid_iwv

ENDIF ; iwv file

;****write LWP file if algorithm specified
IF rets(1) NE '-' THEN BEGIN

 outfile_lwp = out_path + fn_lwp.kkk + '_' + fn_lwp.sss + '_' + fn_lwp.inst + fn_lwp.instnn + '_'+ fn_lwp.lll2 + '_' + fn_lwp.var + '_' + fn_lwp.vnn + '_' + date_string + '.nc'
 dep_l1 = fn_l1b.kkk + '_' + fn_l1b.sss + '_' + fn_l1b.inst + fn_l1b.instnn + '_'+ fn_l1b.lll1 + '_' + fn_l1b.var + '_' + dep_l2a_lwp
 outfid_lwp = NCDF_CREATE(outfile_lwp, /clobber)

;*****define dimensions
 timdid = NCDF_DIMDEF(outfid_lwp,'time', n_time)
 nvdid = NCDF_DIMDEF(outfid_lwp, 'nv', 2)
 rfdid = NCDF_DIMDEF(outfid_iwv, 'n_ret', n_ret)

;*****define variables & attributes
 timid = NCDF_VARDEF(outfid_lwp, 'time',[timdid], /double)
 NCDF_ATTPUT, outfid_lwp, timid, 'units', 'seconds since 1970-01-01 00:00:00 UTC'
 NCDF_ATTPUT, outfid_lwp, timid, 'standard_name', 'time'
 NCDF_ATTPUT, outfid_lwp, timid, 'bounds', 'time_bnds'

 tbnid = NCDF_VARDEF(outfid_lwp, 'time_bnds',[nvdid, timdid], /double)

 latid = NCDF_VARDEF(outfid_lwp, 'lat', /float)
 NCDF_ATTPUT, outfid_lwp, latid, 'units','degree_north'
 NCDF_ATTPUT, outfid_lwp, latid, 'standard_name','latitude'

 lonid = NCDF_VARDEF(outfid_lwp, 'lon', /float)
 NCDF_ATTPUT, outfid_lwp, lonid, 'units','degree_east'
 NCDF_ATTPUT, outfid_lwp, lonid, 'standard_name','longitude'

 zslid = NCDF_VARDEF(outfid_lwp, 'zsl', /float)
 NCDF_ATTPUT, outfid_lwp, zslid, 'units', 'm'
 NCDF_ATTPUT, outfid_lwp, zslid, 'standard_name','altitude'
 NCDF_ATTPUT, outfid_lwp, zslid, 'long_name','altitude above mean sea level'

 aziid = NCDF_VARDEF(outfid_lwp, 'azi',[timdid], /float)
 NCDF_ATTPUT, outfid_lwp, aziid, 'units', 'degree'
 NCDF_ATTPUT, outfid_lwp, aziid, 'standard_name','sensor_azimuth_angle'
 NCDF_ATTPUT, outfid_lwp, aziid, '_FillValue',-999.
 NCDF_ATTPUT, outfid_lwp, aziid, 'comment', '0=North, 90=East, 180=South, 270=West'
 
 eleid = NCDF_VARDEF(outfid_lwp, 'ele',[timdid], /float)
 NCDF_ATTPUT, outfid_lwp, eleid, 'units', 'degree'
 NCDF_ATTPUT, outfid_lwp, eleid, 'long_name','sensor elevation angle'
 NCDF_ATTPUT, outfid_lwp, eleid, '_FillValue',-999.

 releid = NCDF_VARDEF(outfid_lwp, 'ele_ret',[rfdid], /float)
 NCDF_ATTPUT, outfid_lwp, releid, 'units', 'degree'
 NCDF_ATTPUT, outfid_lwp, releid, 'long_name','retrieval elevation angle'
 NCDF_ATTPUT, outfid_lwp, releid, '_FillValue',-999.

 lwpid = NCDF_VARDEF(outfid_lwp, 'clwvi',[timdid], /float)
 NCDF_ATTPUT, outfid_lwp, lwpid, 'units', 'kg m-2'
 NCDF_ATTPUT, outfid_lwp, lwpid, 'standard_name','atmosphere_cloud_liquid_water_content'
 NCDF_ATTPUT, outfid_lwp, lwpid, '_FillValue',-999.
 NCDF_ATTPUT, outfid_lwp, lwpid, 'comment', 'These values denote the vertically integrated amount of condensed water from the surface to TOA.'

 lwpozid = NCDF_VARDEF(outfid_iwv, 'clwvi_off_zenith',[timdid], /float)
 NCDF_ATTPUT, outfid_lwp, lwpozid, 'units', 'kg m-2'
 NCDF_ATTPUT, outfid_lwp, lwpozid, 'long_name','off zenith path integrated liquid water'
 NCDF_ATTPUT, outfid_lwp, lwpozid, '_FillValue', -999.
 NCDF_ATTPUT, outfid_lwp, lwpozid, 'comment', 'These values denote the path integrated amount of condensed water at arbitrary viewing directions in elevation and azimuth. Values not available yet in this data version.'

 lwpeid = NCDF_VARDEF(outfid_lwp, 'clwvi_err',[rfdid], /float)
 NCDF_ATTPUT, outfid_lwp, lwpeid, 'units', 'kg m-2'
 NCDF_ATTPUT, outfid_lwp, lwpeid, 'long_name','standard error of atmosphere cloud liquid water content at retrieval elevation angle'
 NCDF_ATTPUT, outfid_lwp, lwpeid, '_FillValue',-999.

 lwpoffzid  = NCDF_VARDEF(outfid_lwp, 'clwvi_offset_zeroing',[timdid], /float)
 NCDF_ATTPUT, outfid_lwp, lwpoffzid, 'units', 'kg m-2'
 NCDF_ATTPUT, outfid_lwp, lwpoffzid, 'long_name','atmosphere_cloud_liquid_water_content clear sky zeroing offset correction'
 NCDF_ATTPUT, outfid_lwp, lwpoffzid, '_FillValue', -999.
 NCDF_ATTPUT, outfid_lwp, lwpoffzid, 'comment', 'This value has been subtracted from the original clwvi (only zenith!) value to account for instrument calibration drifts. The information is designated for expert user use.'

 lwpoffid  = NCDF_VARDEF(outfid_lwp, 'clwvi_offset',[timdid], /float)
 NCDF_ATTPUT, outfid_lwp, lwpoffid, 'units', 'kg m-2'
 NCDF_ATTPUT, outfid_lwp, lwpoffid, 'long_name','atmosphere_cloud_liquid_water_content offset correction based on brightness temperature offset'
 NCDF_ATTPUT, outfid_lwp, lwpoffid, '_FillValue', -999.
 NCDF_ATTPUT, outfid_lwp, lwpoffid, 'comment', 'This value has been subtracted from the original clwvi value to account for instrument calibration drifts. The information is designated for expert user use.'

 lwpozoffid  = NCDF_VARDEF(outfid_lwp, 'clwvi_off_zenith_offset',[timdid], /float)
 NCDF_ATTPUT, outfid_lwp, lwpozoffid, 'units', 'kg m-2'
 NCDF_ATTPUT, outfid_lwp, lwpozoffid, 'long_name','off zenith path integrated liquid water offset correction based on brightness temperature offset'
 NCDF_ATTPUT, outfid_lwp, lwpozoffid, '_FillValue', -999.
 NCDF_ATTPUT, outfid_lwp, lwpozoffid, 'comment', 'This value has been subtracted from the original path integrated amount of condensed water to account for instrument calibration drifts. The information is designated for expert user use.'

 flaid = NCDF_VARDEF(outfid_lwp, 'flag',[timdid], /short)
 NCDF_ATTPUT, outfid_lwp, flaid, 'long_name','quality control flags'
 NCDF_ATTPUT, outfid_lwp, flaid, 'flag_masks', [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
 NCDF_ATTPUT, outfid_lwp, flaid, 'flag_meanings', 'manual_filter_band_1 manual_filter_band2 manual_filter_band3 rain_flag sanity_receiver_band1 sanity_receiver_band2 sun_in_beam tb_threshold_band1 tb_threshold_band2 tb_threshold_band3 lwp_lwp_threshold'
 NCDF_ATTPUT, outfid_lwp, flaid, '_FillValue', 0
  
 c1 = 'Flags indicate data that the user should only use with care. In cases of doubt, please refer to the contact person. '
 c2 = 'A Fillvalue of 0 means that data has not been flagged. ' 
 c3 = 'Bands refer to the measurement ranges (if applicable) of the microwave radiometer; i.e band 1: 20-30 GHz, band 2: 50-60 GHz, band 3: 90 GHz; '
 c4 = 'tb valid range: ['+ STRING(threshold(0), format = '(f6.2)') + ', ' + STRING(threshold(1), format = '(f6.2)')+'] in K; ' 
 c5 = 'prw valid range: ['+ STRING(threshold(2), format = '(f5.0)') + ', ' + STRING(threshold(3), format = '(f5.0)')+'] in kgm-2; '
 c6 = 'clwvi (zeroing not applied) valid range: ['+ STRING(threshold(4), format = '(f4.1)') + ', ' + STRING(threshold(5), format = '(f4.1)')+'] in kgm-2; '
 
 NCDF_ATTPUT, outfid_lwp, flaid, 'comment', c1+c2+c3+c4+c5+c6

;****global attributes

 NCDF_ATTPUT, outfid_lwp, /GLOBAL, 'Title', 'Microwave radiometer retrieved clwvi'
 NCDF_ATTPUT, outfid_lwp, /GLOBAL, 'Institution', ga.i
 NCDF_ATTPUT, outfid_lwp, /GLOBAL, 'Contact_person', ga.cp
 NCDF_ATTPUT, outfid_lwp, /GLOBAL, 'Source', ga.s
 NCDF_ATTPUT, outfid_lwp, /GLOBAL, 'History', ga.h + 'Retrieval version: ' + rets(1)
 NCDF_ATTPUT, outfid_lwp, /GLOBAL, 'Dependencies', dep_l1
 NCDF_ATTPUT, outfid_lwp, /GLOBAL, 'Conventions', ga.cv
 NCDF_ATTPUT, outfid_lwp, /GLOBAL, 'License', ga.l
 NCDF_ATTPUT, outfid_lwp, /GLOBAL, 'Processing_date', date_create
 NCDF_ATTPUT, outfid_lwp, /GLOBAL, 'Measurement_site', ga.ms
 NCDF_ATTPUT, outfid_lwp, /GLOBAL, 'Author', ga.au
 NCDF_ATTPUT, outfid_lwp, /GLOBAL, 'Comment', 'clwvi errors are given as a function of retrieval elevation angle. Note, observations may occur at different angles (interpolation).'

;****fill all variables with dummy values
 NCDF_CONTROL, outfid_lwp, /fill

;****close define mode
 NCDF_CONTROL, outfid_lwp, /verbose
 NCDF_CONTROL, outfid_lwp, /endef

;****assign data to variables and write to file
 NCDF_VARPUT, outfid_lwp, timid, time 
 NCDF_VARPUT, outfid_lwp, tbnid, time_bnds
 NCDF_VARPUT, outfid_lwp, latid, lat
 NCDF_VARPUT, outfid_lwp, lonid, lon
 NCDF_VARPUT, outfid_lwp, zslid, zsl
 NCDF_VARPUT, outfid_lwp, aziid, az
 NCDF_VARPUT, outfid_lwp, eleid, el
 NCDF_VARPUT, outfid_lwp, releid, ang_ret
 NCDF_VARPUT, outfid_lwp, lwpid, lwp_zen
 NCDF_VARPUT, outfid_lwp, lwpoffid, dlwp
 NCDF_VARPUT, outfid_lwp, lwpozid, lwp_off_zen
 NCDF_VARPUT, outfid_lwp, lwpozoffid, dlwp_off_zen
 NCDF_VARPUT, outfid_lwp, lwpeid, err_lwp
 NCDF_VARPUT, outfid_lwp, lwpoffzid, lwp_offset
 NCDF_VARPUT, outfid_lwp, flaid, flag

;****close netcdf file
 NCDF_CLOSE, outfid_lwp

ENDIF ; lwp file

END

