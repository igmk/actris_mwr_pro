;+
;**********
PRO PL_MK,$
;**********
;KEYWORDS
date=date,$                         ; YYMMDD to be processed 
mwr_meas=mwr_meas,$                 ; MWR measurement (user-defined three letter code)
mwr_dir=mwr_dir,$                   ; directory containing MWR data
home_path=home_path                 ; home directory for MWR_PRO       
; Abstract:
;* read raw radiometer data 
;* make quicklooks of TBs and retrievals (optional)
;* write data to netcdf files (optional)
; Author:
; U. Loehnert
; Date:
; 2011-01-18
; Dependencies:
; -
; Changes:
; 2013-09-10 (UL): 
; 1.) removed faulty lwp offset correction at off-zenith angles
; 2.) Consequently AZTP now plotted with lwp instead of lwp_cor
; 2014-10-24
; 1.) corrected azimuth correction formulation - parentheses were set incorrectly
; 2.) updated get_brt.pro and get_blb.pro for a correct reading of azimuth and elevation angles
; 3.) enhanced search range (in elevation) for sun-influenced TBs by +10deg
;-

!EXCEPT=2

IF N_ELEMENTS(date) EQ 0 THEN BEGIN
 print, 'ABORT: no date specified!'
 STOP
ENDIF

IF N_ELEMENTS(mwr_dir) EQ 0 THEN BEGIN
 print, 'ABORT: no mwr data directory specified!'
 STOP
ENDIF

;***restore user defined parameters
@par_mwr_pro

if verbose GT 0 THEN print, 'entering PL_MK ', date

thresholds = [TB_min, TB_max, IWV_min, IWV_max, LWP_min, LWP_max, T_min, T_max, q_min, q_max] 

;***split date and construct pathes
cal_path = mwr_dir+'/data/calibration/'

month_c = STRMID(date, 2, 2)
date_y = STRMID(date, 0, 2)
date_m = STRMID(date, 0, 4)
date_d = FIX(STRMID(date, 4, 2)) 
date_md = STRMID(date, 2, 4)

;***calculate date of previous day
GET_PREV_DATE, '20'+date, yyyy_prev, mm_prev, dd_prev
date_m_prev = STRMID(yyyy_prev, 2, 2) + mm_prev
date_prev = STRMID(yyyy_prev, 2, 2) + mm_prev + dd_prev

;***create paths
;**check if raw data exists in daily directories
raw_path = mwr_dir+'/data/raw/'+date_m+'/'+date+'/'
raw_path_prev =  mwr_dir+'/data/raw/'+date_m_prev+'/'+date_prev+'/'
filetest = FILE_SEARCH(mwr_dir+'/data/raw/'+date_m+'/'+date+'/*.*')
IF filetest(0) EQ '' THEN BEGIN
;**if no raw data in daily directories found - check if raw data exists in monthly directories
 filetest = FILE_SEARCH(mwr_dir+'/data/raw/'+date_m+'/*.*')
 IF filetest(0) EQ '' THEN BEGIN
  print, 'NO raw data found for '+date
 ENDIF ELSE BEGIN
  raw_path = mwr_dir+'/data/raw/'+date_m+'/'
  raw_path_prev =  mwr_dir+'/data/raw/'+date_m_prev+'/'
 ENDELSE
ENDIF

plot_path1 = mwr_dir+'/plots/level1/'+date_m+'/'
plot_path2 = mwr_dir+'/plots/level2/'+date_m+'/'

data_path0 = mwr_dir+'/data/level0/'+date_m+'/'
data_path1 = mwr_dir+'/data/level1/'+date_m+'/'
data_path2 = mwr_dir+'/data/level2/'+date_m+'/'

;***handle tb offfset correction

index0 = -1
index1 = -1

IF N_ELEMENTS(tb_offset_file) NE 0 THEN BEGIN
;specify offset correction file
 RESTORE, home_path + '/mwr_pro/scripts/'+tb_offset_file
 offset = {date:date_offset, freq:freq_offset, ang:ang_offset, tb:tb_offset}
;**correct only where offset_index is set to 1 
 index0 = WHERE(offset_index EQ 0)
 index1 = WHERE(offset_index EQ 1)
 IF index0(0) NE -1 THEN BEGIN
  offset.tb(*, index0, *) = 0d
 ENDIF
ENDIF ELSE BEGIN
 offset = {date:-999., freq:-999., ang:-999., tb:0d}
ENDELSE

;***EXTRACT daily MWR data

;**assign dummy values
tb = {time:-999., tb:[[-999.], [-999.]], el:-999., az:-999., f:-999., r:-999.}
plot = {time:-999., tb:[[-999.], [-999.]], r:-999.}
l1c = {time:-999., tb:[[-999.], [-999.], [-999.]], el:-999., az:-999., f:-999., r:-999.}
met = {time:-999., p:-999., t:-999., q:-999., r:-999.}
hkd = {time:-999., rec_san:-999.}
irt = {time:-999., tb:-999., el:-999., wavel:-999.}
z_final = -999.

;**RPG data
IF rd_file EQ 'wvl' OR rd_file EQ 'brt' THEN BEGIN
 GET_RPG, raw_path, raw_path_prev, ang_low, ang_high, rd_file, date, tb, plot, l1c, met, hkd, irt, verbose=verbose

;***azimuth correction for HATPRO level1b (transform to "geographical" coordinates)
 az = tb.az
 IF az_cor NE -999. THEN BEGIN
  i_0_180 = WHERE(az GE 0. AND az LE 180.)
  i_180_360 = WHERE(az GT 180. AND az LE 360.)
  IF i_0_180(0) NE -1 THEN az(i_0_180) = az_cor - az(i_0_180)
  IF i_180_360(0) NE -1 THEN az(i_180_360) = 360. + az_cor - az(i_180_360)

  i_neg = WHERE(az LT 0.)
  IF i_neg(0) NE -1 THEN az(i_neg) = az(i_neg) + 360.
 ENDIF

;***azimuth correction for HATPRO level1c (transform to "geographical" coordinates)
 az_l1c = l1c.az
 IF az_cor NE -999. THEN BEGIN
  i_0_180 = WHERE(az_l1c GE 0. AND az_l1c LE 180.)
  i_180_360 = WHERE(az_l1c GT 180. AND az_l1c LE 360.)
  IF i_0_180(0) NE -1 THEN az_l1c(i_0_180) = az_cor - az_l1c(i_0_180)
  IF i_180_360(0) NE -1 THEN az_l1c(i_180_360) = 360. + az_cor - az_l1c(i_180_360)

  i_neg = WHERE(az_l1c LT 0.)
  IF i_neg(0) NE -1 THEN az_l1c(i_neg) = az_l1c(i_neg) + 360.
 ENDIF

;***define elevation angle for RPG data
 el = tb.el

;***perform tb offset correction to level1b data (if specified in par_mwr_pro.pro).

 n_meas = N_ELEMENTS(el)
 n_f = N_ELEMENTS(tb.f)

;**define offset-corrected tb
 tb_l1b_oc = tb.tb
 tb_l1b = tb.tb
 nx_l1b = N_ELEMENTS(tb.tb(*, 0))   
 ny_l1b = N_ELEMENTS(tb.tb(0, *))
 l1b_oc = REPLICATE(0., nx_l1b, ny_l1b)

 n_offset_date = N_ELEMENTS(offset.date)
 n_offset_freq = N_ELEMENTS(offset.freq)
 n_offset_ang = N_ELEMENTS(offset.ang)

;**loop over measurements
 FOR i = 0l, LONG(n_meas-1) DO BEGIN

;*initialize offset indices
  i_offset_t = -99
  i_offset_f = -99

;*find offset correction time index (i_offset_t)
  DATE_TIME_TO_JULDAT, '20'+date, tb.time(i), jd_meas
  IF N_ELEMENTS(offset.date) EQ 1 AND offset.date(0) NE -999. THEN BEGIN
   IF jd_meas GE offset.date(0) THEN i_offset_t = 0
  ENDIF ELSE IF N_ELEMENTS(offset.date) GT 1 THEN BEGIN
   FOR j = 0, n_offset_date-2 DO BEGIN
    ii = WHERE(jd_meas GE offset.date(j) AND jd_meas LT offset.date(j+1))
    IF ii(0) NE -1 THEN i_offset_t = j
   ENDFOR
   ii = WHERE(jd_meas GE offset.date(n_offset_date-1))
   IF ii(0) NE -1 THEN i_offset_t = n_offset_date-1
  ENDIF

;*offset correct tb_l1b_oc
  IF i_offset_t GE 0 THEN BEGIN
   FOR j = 0, n_f-1 DO BEGIN
    ii = WHERE(offset.freq LT tb.f(j)+0.01 AND offset.freq GT tb.f(j)-0.01)
    IF ii(0) NE -1 THEN i_offset_f = ii(0)

    IF i_offset_f GE 0 THEN BEGIN
     IF n_offset_ang EQ 1 THEN BEGIN
      IF el(i) LT offset.ang(0)+0.6 AND el(i) GT offset.ang(0)-0.6 THEN BEGIN
       tb_l1b_oc(j, i) = tb.tb(j, i) - offset.tb(i_offset_t, i_offset_f, 0)
       l1b_oc(j, i) = offset.tb(i_offset_t, i_offset_f, 0)
      ENDIF
     ENDIF ELSE BEGIN
      FOR k = 0, n_offset_ang-2 DO BEGIN
       x = [offset.ang(k), offset.ang(k+1)]
       y = x(SORT(x))
       IF el(i) GE y(0) AND el(i) LE y(1) THEN BEGIN
;**linearly interpolate tb offset correction to elevation angle of measurement
        steig = (offset.tb(i_offset_t, i_offset_f, k+1)-offset.tb(i_offset_t, i_offset_f, k))/(offset.ang(k+1)-offset.ang(k))
        tbx = offset.tb(i_offset_t, i_offset_f, k) + steig*(el(i)-offset.ang(k))
        tb_l1b_oc(j, i) = tb.tb(j, i) - tbx
        l1b_oc(j, i) = tbx 
       ENDIF
      ENDFOR
     ENDELSE ; n_offset_ang > 1 ?
    ENDIF ; i_offset_f
   ENDFOR ; frequency loop
  ENDIF ; i_offset_t
 ENDFOR ; loop over measurements

;***perform tb offset correction to plot data (if specified in par_mwr_pro.pro).

 n_meas_p = N_ELEMENTS(plot.r)

;**define offset-corrected tb
 tb_p_oc = plot.tb
 nx_p = N_ELEMENTS(plot.tb(*, 0))   
 ny_p = N_ELEMENTS(plot.tb(0, *))
 plot_oc = REPLICATE(0., nx_p, ny_p)

 n_offset_date = N_ELEMENTS(offset.date)
 n_offset_freq = N_ELEMENTS(offset.freq)
 n_offset_ang = N_ELEMENTS(offset.ang)

;**loop over measurements
 FOR i = 0l, LONG(n_meas_p-1) DO BEGIN

;*initialize offset indices
  i_offset_t = -99
  i_offset_f = -99

;*find offset correction time index (i_offset_t)
  DATE_TIME_TO_JULDAT, '20'+date, plot.time(i), jd_meas
  IF N_ELEMENTS(offset.date) EQ 1 AND offset.date(0) NE -999. THEN BEGIN
   IF jd_meas GE offset.date(0) THEN i_offset_t = 0
  ENDIF ELSE IF N_ELEMENTS(offset.date) GT 1 THEN BEGIN
   FOR j = 0, n_offset_date-2 DO BEGIN
    ii = WHERE(jd_meas GE offset.date(j) AND jd_meas LT offset.date(j+1))
    IF ii(0) NE -1 THEN i_offset_t = j
   ENDFOR
   ii = WHERE(jd_meas GE offset.date(n_offset_date-1))
   IF ii(0) NE -1 THEN i_offset_t = n_offset_date-1
  ENDIF

;*offset correct plot.tb
  IF i_offset_t GE 0 THEN BEGIN
   FOR j = 0, n_f-1 DO BEGIN
    ii = WHERE(offset.freq LT tb.f(j)+0.01 AND offset.freq GT tb.f(j)-0.01)
    IF ii(0) NE -1 THEN i_offset_f = ii(0)

    IF i_offset_f GE 0 THEN BEGIN
     IF n_offset_ang EQ 1 THEN BEGIN
      IF offset.ang(0) LT ang_high AND offset.ang(0) GT ang_low THEN BEGIN
       tb_p_oc(j, i) = plot.tb(j, i) - offset.tb(i_offset_t, i_offset_f, 0)
       plot_oc(j, i) = offset.tb(i_offset_t, i_offset_f, 0)
      ENDIF
     ENDIF ELSE BEGIN
      FOR k = 0, n_offset_ang-2 DO BEGIN
       x = [offset.ang(k), offset.ang(k+1)]
       y = x(SORT(x))
       el_p = (ang_high+ang_low)/2.
       IF el_p GE y(0) AND el_p LE y(1) THEN BEGIN
;**linearly interpolate tb offset correction to elevation angle of measurement
        steig = (offset.tb(i_offset_t, i_offset_f, k+1)-offset.tb(i_offset_t, i_offset_f, k))/(offset.ang(k+1)-offset.ang(k))
        tbx = offset.tb(i_offset_t, i_offset_f, k) + steig*(el_p-offset.ang(k))
        tb_p_oc(j, i) = plot.tb(j, i) - tbx
        plot_oc(j, i) = tbx 
       ENDIF
      ENDFOR
     ENDELSE ; n_offset_ang > 1 ?
    ENDIF ; i_offset_f
   ENDFOR ; frequency loop
  ENDIF ; i_offset_t
 ENDFOR ; loop over measurements


;***perform tb offset correction to level1c data (if specified in par_mwr_pro.pro).

 n_1c = N_ELEMENTS(l1c.time)
 n_f1c = N_ELEMENTS(l1c.f)
 n_a1c = N_ELEMENTS(l1c.el)

;**define offset-corrected tb
 tb_l1c_oc = l1c.tb
 tb_l1c = l1c.tb
 nx_l1c = N_ELEMENTS(l1c.tb(*, 0, 0))   
 ny_l1c = N_ELEMENTS(l1c.tb(0, *, 0))
 nz_l1c = N_ELEMENTS(l1c.tb(0, 0, *))
 l1c_oc = REPLICATE(0., nx_l1c, ny_l1c, nz_l1c)

 n_offset_date = N_ELEMENTS(offset.date)
 n_offset_freq = N_ELEMENTS(offset.freq)
 n_offset_ang = N_ELEMENTS(offset.ang)

 FOR i = 0, n_1c-1 DO BEGIN

;**initialize offset indices
  i_offset_t = -99
  i_offset_f = -99
  i_offset_a = -99

;**find offset correction time index (i_offset_t)
  DATE_TIME_TO_JULDAT, '20'+date, l1c.time(i), jd_meas
  IF N_ELEMENTS(offset.date) EQ 1 AND offset.date(0) NE -999. THEN BEGIN
   IF jd_meas GT offset.date(0) THEN i_offset_t = 0
  ENDIF ELSE IF N_ELEMENTS(offset.date) GT 1 THEN BEGIN
   FOR j = 0, n_offset_date-2 DO BEGIN
    ii = WHERE(jd_meas GE offset.date(j) AND jd_meas LT offset.date(j+1))
    IF ii(0) NE -1 THEN i_offset_t = j
   ENDFOR
   ii = WHERE(jd_meas GE offset.date(n_offset_date-1))
   IF ii(0) NE -1 THEN i_offset_t = n_offset_date-1
  ENDIF

;**offset correct tb
  IF i_offset_t GE 0 THEN BEGIN

   FOR j = 0, n_f1c-1 DO BEGIN
    ii = WHERE(offset.freq LT l1c.f(j)+0.01 AND offset.freq GT l1c.f(j)-0.01)
    IF ii(0) NE -1 THEN i_offset_f = ii(0)
    IF i_offset_f GE 0 THEN BEGIN
 
     FOR k = 0, n_a1c-1 DO BEGIN
      jj = WHERE(offset.ang LT l1c.el(k)+0.6 AND offset.ang GT l1c.el(k)-0.6)
      IF jj(0) NE -1 THEN BEGIN 
       i_offset_a = jj(0)
       l1c_oc(j, k, i) = offset.tb(i_offset_t, i_offset_f, i_offset_a)       
       tb_l1c_oc(j, k, i) = tb_l1c_oc(j, k, i) - offset.tb(i_offset_t, i_offset_f, i_offset_a)
      ENDIF 
     ENDFOR ; elevation loop

    ENDIF

   ENDFOR ; frequency loop
  ENDIF ; i_offset_t

 ENDFOR ; n_1c loop

;**RESCOM data
ENDIF ELSE IF rd_file EQ 'dec' THEN BEGIN
 GET_RES, raw_path, raw_path_prev, ang_low, ang_high, date, tb, plot, verbose=verbose
 az = REPLICATE(az_fix, N_ELEMENTS(tb.time))
 el = tb.el

;**RADIOMETRICS data
ENDIF ELSE IF rd_file EQ 'rad' THEN BEGIN
; GET_RAD still needs to be coded ...

;**OTHER systems ... still need to be coded

ENDIF

IF index1(0) NE -1 THEN BEGIN
 comment = ' TB offset correction applied to channels (GHz):'

 nfi1 = N_ELEMENTS(index1) 
 FOR ind1 = 0, nfi1-1 DO BEGIN
  comment = comment + ', ' + STRING(tb.f(index1(ind1)), format = '(f5.2)')
 ENDFOR

ENDIF

;**put met data (if available) on tb time
n_tb = LONG(N_ELEMENTS(tb.time))
tb_t = REPLICATE(-999., n_tb)
tb_p = REPLICATE(-999., n_tb)
tb_rh = REPLICATE(-999., n_tb)
tb_q = REPLICATE(-999., n_tb)

IF N_ELEMENTS(met.time) GT 1 THEN BEGIN 
 
 print, 'PL_MK: find matching MET values to TB times'
 FOR i = 0l, n_tb-1l, 100 DO BEGIN
  CLOSEST, tb.time(i), met.time, dmin, ind
  IF dmin LT 0.25 THEN BEGIN
   tb_t(i) = met.t(ind)
   tb_p(i) = met.p(ind)
   tb_rh(i) = met.q(ind)
  ENDIF
 ENDFOR

 RH_TO_ABSHUM, tb_t, tb_rh, tb_q 

ENDIF

;**put irt data (if available) on tb time
tb_irt = REPLICATE(-999., 1, n_tb)
el_irt = REPLICATE(-999., n_tb)

IF N_ELEMENTS(irt.time) GT 1 THEN BEGIN 

 n_wavel_irt = N_ELEMENTS(irt.wavel)
 tb_irt = REPLICATE(-999., n_wavel_irt, n_tb)
 el_irt = REPLICATE(-999., n_tb)
 
 print, 'PL_MK: find matching IRT TBs to TB times'
 FOR i = 0l, n_tb-1l DO BEGIN
  CLOSEST, tb.time(i), irt.time, dmin, ind
  IF dmin LT 10./3600. THEN BEGIN
   tb_irt(*, i) = irt.tb(*, ind)
   el_irt(i) = irt.el(ind)
  ENDIF
 ENDFOR

ENDIF

;***read ceilometer data if available

IF N_ELEMENTS(ceilo_type) GT 0 THEN BEGIN

 IF ceilo_type EQ 'ct25' THEN BEGIN

  ceilo_file_x = ceilo_path + '20'+date_m+'/*' + date_md + '*.DAT*'
  ceilo_file = FILE_SEARCH(ceilo_file_x)  

  time_ceilo = REPLICATE(-999., 2)
  bscat = REPLICATE(-999., 2, 2)
  z_ceilo = REPLICATE(-999., 2)
 
  IF ceilo_file(0) EQ '' AND verbose THEN print, 'PL_MK: no ceilometer files found' 
  IF ceilo_file(0) NE '' THEN BEGIN

   n_ceilo_files = N_ELEMENTS(ceilo_file)
   IF verbose THEN print, 'PL_MK: number of CT25K ceilometer files found: ', n_ceilo_files

   FOR i_ceilo = 0, n_ceilo_files-1 DO BEGIN

    READ_RAW_CT25K,$
     ceilo_file(i_ceilo),       $ ; name of input file
     time_ceilo1,               $ ; time as julian day, DBLARR( numprof )
     z_ceilo,                   $ ; height in m above ground, FLTARR( nz )
     detection_status, $ ; number of detected cloud levels (0..3) or flag for vertical visibility (4)
                       $ ; or something else (5), or data missing or suspect (99 - the original ascii is '/') INTARR(numprof)
     status_flag,      $ ; OK ('0'), Warning ('W') or Alarm ('A') , BYTARR(numprof)
     first_cbh, second_cbh, third_cbh,$ ; height above gnd (m) of detected cloud layers, FLTARR(numprof)
     vertical_vis, max_bscat,         $ ; height above gnd (m) of derived vertical visibility, FLTARR(numprof), NEW Nov 2007
     laser_pulse_energy,              $ ; laser pulse energy (% of optimum ?), FLTARR(numprof)
     laser_temperature,               $ ; laser temperature (°C), FLTARR(numprof)
     receiver_sensitivity,            $ ; receiver sensitivity (% of optimum ?), FLTARR(numprof)
     window_transm_est,               $ ; window contamination (0..2500mV), FLTARR(numprof)
     tilt_angle,                      $ ; tilt angle of the instrument (° against vertical) , FLTARR(numprof)
     background_light,                $ ; background light (0..2500mV), FLTARR(numprof)
     sum_backscatter,                 $ ; integrated backscatter from instrument (units ?), FLTARR(numprof)
     bscat1,                          $ ; backscatter profiles in (1e-4/km/sr), FLTARR(nz,numprof)
     measurement_parameters,          $; information about actual device setup (string[stat_len]), BYTARR(stat_len,numprof)
     status_string,                   $ ; current alarms and warnings Coded as a hex string (string[stat_len]), BYTARR(stat_len,numprof)
     scale_factor=scale_factor,       $ ; first parameter in the status line of the ceilometer (intarr(numprof))
     ceilometername=ceilometername,   $ ; every data record starts with 'CL<unit_number><software_level><message_number><message_subclass>' the procedure returns
     cloud_data_string=cloud_data,    $ ; line 2 of the last data message = cloud data in the form <N><W|A> <cld_hgt[0]> <cld_hgt[1]> <cld_hgt[2]> <statusstring>
     ceilo_para_string=ceilo_para,    $ ; line 3 of the last data message = instrument paramters in the form <scale_factor> <bksc_prof_res> <length_prof_in_sample
     verbose=0 ;

     bscat1 = bscat1/1000. ; transform to 1e-4/m/sr

     IF i_ceilo EQ 0 THEN BEGIN
      time_ceilo = time_ceilo1
      bscat = bscat1
     ENDIF ELSE BEGIN
      time_ceilo = [time_ceilo, time_ceilo1]
      bscat = [[bscat], [bscat1]] 
     ENDELSE

   ENDFOR ; n_ceilo_files

  ENDIF ; ceilo files found

 ENDIF ELSE IF ceilo_type EQ 'jen' THEN BEGIN

  print, 'PL_MK: Data processing for JENOPTIK ceilometer still needs to be coded ...' 
  stop  

 ENDIF

ENDIF ; ceilo_type specified

;***APPLY RETRIEVAL to extracted data
;**level2a data products: IWV, LWP, ZWD, ATT, ...
;(vertically integrated)

iwv = REPLICATE(-999., N_ELEMENTS(tb.time))
lwp = REPLICATE(-999., N_ELEMENTS(tb.time))
wdl = REPLICATE(-999., N_ELEMENTS(tb.time))
taa = REPLICATE(-999., N_ELEMENTS(tb.time))
diwv = REPLICATE(-999., N_ELEMENTS(tb.time))
dlwp = REPLICATE(-999., N_ELEMENTS(tb.time))
dwdl = REPLICATE(-999., N_ELEMENTS(tb.time))
dtaa = REPLICATE(-999., N_ELEMENTS(tb.time))

IF N_ELEMENTS(algo_iwv) NE 0 THEN BEGIN
 IF ret_file EQ 'RET' THEN ret_files = FILE_SEARCH(home_path + '/mwr_pro/retrievals/iwv_'+algo_iwv+'*.RET')
 IF ret_file EQ 'nc' THEN ret_files = FILE_SEARCH(home_path + '/mwr_pro/retrievals/iwv_'+algo_iwv+'*.nc')
 IF verbose THEN print, 'PL_MK: Retrieving IWV'
 IF ret_files(0) EQ '' THEN BEGIN
  print, 'PL_MK: no IWV retrieval file found'
  GOTO, SKIP_1
 ENDIF

 RET_2A_ARBANG, ret_files, ret_file=ret_file, tb.time, tb_l1b, tb_l1b_oc, tb.el, tb.f, date, iwv, diwv, err_iwv,$
                ang_ret_2a, ret_par=ret_par_iwv, tb_t=tb_t, tb_q=tb_q, tb_p=tb_p, verbose=verbose
ENDIF 
SKIP_1:

IF N_ELEMENTS(algo_lwp) NE 0 THEN BEGIN
 ret_files = FILE_SEARCH(home_path + '/mwr_pro/retrievals/lwp_'+algo_lwp+'*') 
 IF ret_files(0) EQ '' THEN BEGIN
  print, 'PL_MK: no LWP retrieval file found'
  GOTO, SKIP_2
 ENDIF
 IF verbose THEN print, 'PL_MK: Retrieving LWP'
 RET_2A_ARBANG, ret_files, ret_file=ret_file, tb.time, tb_l1b, tb_l1b_oc, tb.el, tb.f, date, lwp, dlwp, err_lwp,$
                ang_ret_2a, ret_par=ret_par_lwp, tb_t=tb_t, tb_q=tb_q, tb_p=tb_p, verbose=verbose
ENDIF
SKIP_2:

IF N_ELEMENTS(algo_wdl) NE 0 THEN BEGIN
 ret_files = FILE_SEARCH(home_path + '/mwr_pro/retrievals/wdl_'+algo_wdl+'*')  
 IF ret_files(0) EQ '' THEN BEGIN
  print, 'PL_MK: no WDL retrieval file found'
  GOTO, SKIP_3
 ENDIF
 IF verbose THEN print, 'PL_MK: Retrieving WDL'
 RET_2A_ARBANG, ret_files, ret_file=ret_file, tb.time, tb_1lb, tb_l1b_oc, tb.el, tb.f, date, wdl, dwdl, err_wdl,$
                ang_ret_2a, ret_par=ret_par_wdl, tb_t=tb_t, tb_q=tb_q, tb_p=tb_p, verbose=verbose
ENDIF
SKIP_3:

IF N_ELEMENTS(algo_taa) NE 0 THEN BEGIN
 FOR i = 0, N_ELEMENTS(algo_taa)-1 DO BEGIN
  ret_files = FILE_SEARCH(home_path + '/mwr_pro/retrievals/taa_'+algo_taa(i)+'*')
  IF ret_files(0) EQ '' THEN BEGIN
   print, 'PL_MK: no TAA retrieval file found'
   GOTO, SKIP_4
  ENDIF
  IF verbose THEN print, 'PL_MK: Retrieving TAA'
  RET_2A_ARBANG, ret_files, ret_file=ret_file, tb.time, tb_l1b, tb_l1b_oc, tb.el, tb.f, date, taa, dtaa, err_taa,$
                 ang_ret_2a, ret_par=ret_par_taa, tb_t=tb_t, tb_q=tb_q, tb_p=tb_p, verbose=verbose
  IF i EQ 0 THEN taa_a = REPLICATE(-999., N_ELEMENTS(algo_taa), N_ELEMENTS(taa))
  taa_a(i, *) = taa
 ENDFOR
SKIP_4:
 IF N_ELEMENTS(f_taa) EQ 0 THEN BEGIN
  print, 'No frequencies for TAA specified in par_mwr.pro - ABORTING'
  stop
 ENDIF
ENDIF

;**level2b data products: T-profile, q-profile
;(profiles using data from one observing direction, i.e. zenith)

IF N_ELEMENTS(algo_tze) NE 0 THEN BEGIN 
 ret_files = FILE_SEARCH(home_path + '/mwr_pro/retrievals/tze_'+algo_tze+'*')
 IF ret_files(0) EQ '' THEN BEGIN
  print, 'PL_MK: no TZE retrieval file found'
  GOTO, SKIP_5
 ENDIF
 RET_2B, ret_files, ret_file=ret_file, tb_l1b, tb_l1b_oc, tb.el, tb.f, tb.time, date, z_final, T2b, dT2b,$
         err_T2b, ang_ret_T2b, ret_par=ret_par_tze, tb_t=tb_t, tb_q=tb_q, tb_p=tb_p, verbose=verbose
                
 ENDIF
SKIP_5:
IF N_ELEMENTS(algo_hze) NE 0 THEN BEGIN
 ret_files = FILE_SEARCH(home_path + '/mwr_pro/retrievals/hze_'+algo_hze+'*')
 IF ret_files(0) EQ '' THEN BEGIN
  print, 'PL_MK: no HZE retrieval file found'
  GOTO, SKIP_6
 ENDIF
 RET_2B, ret_files, ret_file=ret_file, tb_l1b, tb_l1b_oc, tb.el, tb.f, tb.time, date, z_final, q2b, dq2b,$
         err_q2b, ang_ret_q2b, ret_par=ret_par_hze, tb_t=tb_t, tb_q=tb_q, tb_p=tb_p, verbose=verbose
ENDIF
SKIP_6:

;**level2c data products: T-profile, RH-profile
;(profiles using data from multiple observing directions under the assumption of horz. homogeneity)
IF N_ELEMENTS(algo_tel) NE 0 AND N_ELEMENTS(algo_hze) NE 0 THEN BEGIN
 ret_files = FILE_SEARCH(home_path + '/mwr_pro/retrievals/tel_'+algo_tel+'.nc')
 IF ret_files(0) EQ '' THEN BEGIN
  print, 'PL_MK: no TEL retrieval file found'
  GOTO, SKIP_7
 ENDIF
 par2c = {time:-999., z:-999., T:-999., q:-999., rh:-999., Tpot:-999., Tepot:-999.}
 IF l1c.time(0) GT 0. AND N_ELEMENTS(q2b) GT 0 THEN BEGIN
  RET_2C, ret_files(0), ret_file=ret_file, l1c.time, tb_l1c, tb_l1c_oc, l1c.el, l1c.f, met, tb.time, dq2b+q2b, q2b,$
          data_path2, file_naming_2b_hze, met.p, mwr_dir, date, mwr_meas, par2c, dpar2c, err_T2c, ret_par=ret_par_tel, verbose=verbose
  T2c = par2c.T
  dT2c = dpar2c.T
 ENDIF
 IF N_ELEMENTS(q2b) EQ 0 THEN print, 'PL_MK: no 2C retrieval calculated due to missing 2b data'
 IF l1c.time(0) LT 0. THEN print, 'PL_MK: no 2C retrieval calculated due to missing 1c data'
ENDIF 
SKIP_7:

;***COMPUTE solar zenith and azimuth angles
;*compute Julian day
DATE_TIME_TO_JULDAT, '20'+date, 0., jday2
DATE_TIME_TO_JULDAT, '20'+date_y+'0101', 0., jday1
jday = jday2 - jday1
ZENSUN, LONG(jday), tb.time, latitude, longitude, sol_zen, sol_az, sol_fac, local=0
sol_el = 90. - sol_zen
ii = WHERE(sol_az LT 0.)
IF ii(0) NE -1 THEN sol_az(ii) = 360. + sol_az(ii)

sol_el_max = MAX(sol_el)
sunrise = 0.
sunset = 24.
i_sun = WHERE(sol_el GT 0., n_sun)

IF i_sun(0) NE -1 THEN BEGIN
 i_sunrise = i_sun(0)
 sunrise = tb.time(i_sunrise)
 i_sunset = i_sun(n_sun-1)
 sunset = tb.time(i_sunset)
ENDIF

;***QUALITY CONTROL: FLAG extracted TBs and retrieval products
;--> level0, 1, 2b and 2c flags are generated

flag_1b = REPLICATE(0, N_ELEMENTS(tb.time))
flag_2a = REPLICATE(0, N_ELEMENTS(tb.time))
flag_T2b = REPLICATE(0, N_ELEMENTS(tb.time))
flag_q2b = REPLICATE(0, N_ELEMENTS(tb.time))
flag_1c = REPLICATE(0, N_ELEMENTS(l1c.time))
flag_2c = REPLICATE(0, N_ELEMENTS(l1c.time))

;- FLAGS are set as bits
;Bit1: MANUAL FILTER Band1 (user edited filter_*.dat file)
;Bit2: MANUAL FILTER Band2 (user edited filter_*.dat file)
;Bit3: MANUAL FILTER Band3 (user edited filter_*.dat file)
;Bit4: RAIN FLAG (RPG specific)
;Bit5: SANITY RECEIVER Band1 (RPG specific)
;Bit6: SANITY RECEIVER Band2 (RPG specific)
;Bit7: solar flag
;Bit8: TB THRESHOLD Band1 (set in par_mwr_pro.pro)
;Bit9: TB THRESHOLD Band2 (set in par_mwr_pro.pro)
;Bit10: TB THRESHOLD Band3 (set in par_mwr_pro.pro)
;Bit11: retrieved LWP/IWV threshold (set in par_mwr_pro.pro)
;Bit12: retrieved TEMPERATURE threshold (set in par_mwr_pro.pro)

;**ALL MWR: 
;*read FILTER FILE (edited manually)
filter_file = FILE_SEARCH(mwr_dir+'/filter_'+mwr_meas+'.dat')
IF filter_file EQ '' THEN BEGIN

 IF verbose THEN print, 'No manual filter file found'

ENDIF ELSE BEGIN

 OPENR, unit_filter, filter_file, /GET_LUN
 s = ''
 n_man = 0
 FOR i = 0, 20 DO READF, unit_filter, s
 WHILE STRMID(s, 0, 19) NE 'date of last change' DO BEGIN
  READF, unit_filter, s, FORMAT = '(a27)'
  IF STRMID(s, 0, 6) EQ date THEN BEGIN
   n_man = FIX(STRMID(s, 6, 3))
   filter_man = FLTARR(2, 3, n_man)
   f_start = FLOAT(STRMID(s, 10, 5))
   f_end = FLOAT(STRMID(s, 16, 5))
   FOR j = 0, 2 DO BEGIN
    IF STRMID(s, 22+j*2, 1) EQ '1' THEN BEGIN 
     filter_man(0, j, 0) = f_start  
     filter_man(1, j, 0) = f_end
    ENDIF
   ENDFOR

   FOR k = 0, n_man-2 DO BEGIN
    READF, unit_filter, s, FORMAT = '(a27)'
    f_start = FLOAT(STRMID(s, 10, 5))
    f_end = FLOAT(STRMID(s, 16, 5))
    FOR j = 0, 2 DO BEGIN
     IF STRMID(s, 22+j*2, 1) EQ '1' THEN BEGIN
      filter_man(0, j, k+1) = f_start
      filter_man(1, j, k+1) = f_end 
     ENDIF
    ENDFOR
   ENDFOR

  ENDIF 
 ENDWHILE

 READF, unit_filter, s, FORMAT = '(a6)'
 dolc = s
 
 FREE_LUN, unit_filter 
  
ENDELSE

;*assign FILTER FILE flags

FOR i = 0, n_man-1 DO BEGIN

 ii = WHERE(tb.time GE filter_man(0, 0, i) AND tb.time LE filter_man(1, 0, i))
 IF ii(0) NE -1 THEN BEGIN
  flag_1b(ii) = 1
  flag_2a(ii) = 1
  flag_T2b(ii) = 1
  flag_q2b(ii) = 1
 ENDIF

 ii = WHERE(l1c.time GE filter_man(0, 0, i) AND l1c.time LE filter_man(1, 0, i))
 IF ii(0) NE -1 THEN BEGIN
  flag_1c(ii) = 1
  flag_2c(ii) = 1
 ENDIF 

 ii = WHERE(tb.time GE filter_man(0, 1, i) AND tb.time LE filter_man(1, 1, i))
 IF ii(0) NE -1 THEN BEGIN
  flag_1b(ii) = flag_1b(ii) OR 2
  flag_2a(ii) = flag_2a(ii) OR 2
  flag_T2b(ii) = flag_T2b(ii) OR 2
  flag_q2b(ii) = flag_q2b(ii) OR 2
 ENDIF

 ii = WHERE(l1c.time GE filter_man(0, 1, i) AND l1c.time LE filter_man(1, 1, i))
 IF ii(0) NE -1 THEN BEGIN
  flag_1c(ii) = flag_1c(ii) OR 2
  flag_2c(ii) = flag_2c(ii) OR 2
 ENDIF 

 ii = WHERE(tb.time GE filter_man(0, 2, i) AND tb.time LE filter_man(1, 2, i))
 IF ii(0) NE -1 THEN BEGIN
  flag_1b(ii) = flag_1b(ii) OR 4
  flag_2a(ii) = flag_2a(ii) OR 4
  flag_T2b(ii) = flag_T2b(ii) OR 4
  flag_q2b(ii) = flag_q2b(ii) OR 4
 ENDIF

 ii = WHERE(l1c.time GE filter_man(0, 2, i) AND l1c.time LE filter_man(1, 2, i))
 IF ii(0) NE -1 THEN BEGIN
  flag_1c(ii) = flag_1c(ii) OR 4
  flag_2c(ii) = flag_2c(ii) OR 4
 ENDIF 

ENDFOR

;**only RPG radiometers - automatic rain flags
IF rd_file EQ 'wvl' OR rd_file EQ 'brt' THEN BEGIN

 ii = WHERE(tb.r GT 0)
 IF ii(0) NE -1 THEN BEGIN
  flag_1b(ii) = flag_1b(ii) OR 8
  flag_2a(ii) = flag_2a(ii) OR 8
  flag_T2b(ii) = flag_T2b(ii) OR 8
  flag_q2b(ii) = flag_q2b(ii) OR 8
 ENDIF

 ii = WHERE(l1c.r AND 1)
 IF ii(0) NE -1 THEN BEGIN
  flag_1c(ii) = flag_1c(ii) OR 8
  flag_2c(ii) = flag_2c(ii) OR 8
 ENDIF
ENDIF

;**only RPG radiometers - check .HKD data for sanity checks
;* time intervals of +-15 min around level0, 1, 2b and 2c data are checked for sanity

IF rd_file EQ 'wvl' OR rd_file EQ 'brt' THEN BEGIN
 IF verbose THEN print, 'PL_MK: checking HKD files' 
 n_hkd = N_ELEMENTS(hkd.time)

 IF n_hkd GT 1 THEN BEGIN
 
;* assign flags (Band1)
  ii = WHERE(hkd.rec_san(0, *) GT 0, n_ii)

  IF ii(0) NE -1 THEN BEGIN
   FOR i = 0l, n_ii-1 DO BEGIN
    
    jj = WHERE(tb.time GT hkd.time(ii(i))-0.25 AND tb.time LT hkd.time(ii(i))+0.25)
    IF jj(0) NE -1 THEN BEGIN

     flag_1b(ii) = flag_1b(ii) OR 16
     flag_2a(ii) = flag_2a(ii) OR 16
     flag_T2b(ii) = flag_T2b(ii) OR 16
     flag_q2b(ii) = flag_q2b(ii) OR 16

    ENDIF
   
    jj = WHERE(l1c.time GT hkd.time(ii(i))-0.25 AND l1c.time LT hkd.time(ii(i))+0.25)
    IF jj(0) NE -1 THEN BEGIN
     flag_1c(ii) = flag_1c(ii) OR 16
     flag_2c(ii) = flag_2c(ii) OR 16
    ENDIF
         
   ENDFOR ; number of failed sanity checks
  ENDIF ; sanity checks failed?

;* assign flags (Band2)
  ii = WHERE(hkd.rec_san(1, *) GT 0, n_ii)

  IF ii(0) NE -1 THEN BEGIN
   FOR i = 0l, n_ii-1 DO BEGIN
    
    jj = WHERE(tb.time GT hkd.time(ii(i))-0.25 AND tb.time LT hkd.time(ii(i))+0.25)
    IF jj(0) NE -1 THEN BEGIN

     flag_1b(jj) = flag_1b(jj) OR 32
     flag_2a(jj) = flag_2a(jj) OR 32
     flag_T2b(jj) = flag_T2b(jj) OR 32
     flag_q2b(jj) = flag_q2b(jj) OR 32

    ENDIF

    jj = WHERE(l1c.time GT hkd.time(ii(i))-0.25 AND l1c.time LT hkd.time(ii(i))+0.25)
    IF jj(0) NE -1 THEN BEGIN
     flag_1c(jj) = flag_1c(jj) OR 32
     flag_2c(jj) = flag_2c(jj) OR 32
    ENDIF
            
   ENDFOR ; number of failed sanity checks
  ENDIF ; sanity checks failed?

 ENDIF ; hkd file exists?

ENDIF ; rd_file?

;**ALL radiometers - check threshold defined in par_mwr_pro.pro
;Bit7: solar flag
;Bit8: TB THRESHOLD Band1 (set in par_mwr_pro.pro)
;Bit9: TB THRESHOLD Band2 (set in par_mwr_pro.pro)
;Bit10: TB THRESHOLD Band3 (set in par_mwr_pro.pro)
;Bit11: retrieved LWP/IWV threshold (set in par_mwr_pro.pro)
;Bit12: retrieved TEMPERATURE threshold (set in par_mwr_pro.pro)

;**solar flag
IF N_ELEMENTS(saf) EQ 0 THEN saf = 0.
i_check_sol = WHERE(tb.el LE sol_el_max+10. AND tb.time GT sunrise and tb.time LT sunset, n_check_sol)
IF i_check_sol(0) NE -1 THEN BEGIN
 FOR i = 0, n_check_sol-1 DO BEGIN
  j = i_check_sol(i)
  IF (tb.el(j) GE sol_el(j)-saf AND tb.el(j) LE sol_el(j)+saf) AND (az(j) GE sol_az(j)-saf AND az(j) LE sol_az(j)+saf) THEN BEGIN
   flag_1b(j) = flag_1b(j) OR 64
   flag_2a(j) = flag_2a(j) OR 64
   flag_T2b(j) = flag_T2b(j) OR 64
   flag_q2b(j) = flag_q2b(j) OR 64
  ENDIF
 ENDFOR
ENDIF
i_check_sol = WHERE(l1c.el LE sol_el_max+10. AND l1c.time GT sunrise and l1c.time LT sunset, n_check_sol)
IF i_check_sol(0) NE -1 THEN BEGIN
 FOR i = 0, n_check_sol-1 DO BEGIN
  j = i_check_sol(i)
  IF (l1c.el(j) GE sol_el(j)-saf AND l1c.el(j) LE sol_el(j)+saf) AND (l1c.az(j) GE sol_az(j)-saf AND l1c.az(j) LE sol_az(j)+saf) THEN BEGIN
   flag_1c(j) = flag_1b(j) OR 64
   flag_2c(j) = flag_2a(j) OR 64
  ENDIF
 ENDFOR
ENDIF

;**TB THRESHOLD Band1
ifreq_check = WHERE(tb.f GT 20. AND tb.f LE 40., nfreq_check)
IF ifreq_check(0) NE -1 THEN BEGIN

 IF verbose THEN print, 'PL_MK: checking TB THRESHOLD Band1'
 FOR i = 0, nfreq_check-1 DO BEGIN

  jj = WHERE(tb.tb(ifreq_check(i), *) GT TB_max OR tb.tb(ifreq_check(i), *) LT TB_min)
  IF jj(0) NE -1 THEN BEGIN
   flag_1b(jj) = flag_1b(jj) OR 128
   flag_2a(jj) = flag_2a(jj) OR 128
   flag_T2b(jj) = flag_T2b(jj) OR 128   
   flag_q2b(jj) = flag_q2b(jj) OR 128   
  ENDIF

  IF l1c.time(0) NE -999. THEN BEGIN
   FOR j = 0, N_ELEMENTS(l1c.el)-1 DO BEGIN

    jj = WHERE(l1c.tb(ifreq_check(i), j, *) GT TB_max OR l1c.tb(ifreq_check(i), j, *) LT TB_min)
    IF jj(0) NE -1 THEN BEGIN
     flag_1c(jj) = flag_1b(jj) OR 128
     flag_2c(jj) = flag_2c(jj) OR 128   
    ENDIF
   
   ENDFOR
  ENDIF

 ENDFOR

ENDIF

;**TB THRESHOLD Band2
ifreq_check = WHERE(tb.f GT 40. AND tb.f LE 70., nfreq_check)
IF ifreq_check(0) NE -1 THEN BEGIN

 IF verbose THEN print, 'PL_MK: checking TB THRESHOLD Band2'
 FOR i = 0, nfreq_check-1 DO BEGIN

  jj = WHERE(tb.tb(ifreq_check(i), *) GT TB_max OR tb.tb(ifreq_check(i), *) LT TB_min)
  IF jj(0) NE -1 THEN BEGIN
   flag_1b(jj) = flag_1b(jj) OR 256
   flag_2a(jj) = flag_2a(jj) OR 256
   flag_T2b(jj) = flag_T2b(jj) OR 256   
   flag_q2b(jj) = flag_q2b(jj) OR 256   
  ENDIF

  IF l1c.time(0) NE -999. THEN BEGIN 
   FOR j = 0, N_ELEMENTS(l1c.el)-1 DO BEGIN

    jj = WHERE(l1c.tb(ifreq_check(i), j, *) GT TB_max OR l1c.tb(ifreq_check(i), j, *) LT TB_min) 
    IF jj(0) NE -1 THEN BEGIN
     flag_1c(jj) = flag_1b(jj) OR 256
     flag_2c(jj) = flag_2c(jj) OR 256  
    ENDIF
   
   ENDFOR
  ENDIF

 ENDFOR

ENDIF

;**TB THRESHOLD Band3
ifreq_check = WHERE(tb.f LT 20. OR tb.f GT 70., nfreq_check)
IF ifreq_check(0) NE -1 THEN BEGIN

 IF verbose THEN print, 'PL_MK: checking TB THRESHOLD Band3'

 FOR i = 0, nfreq_check-1 DO BEGIN

  jj = WHERE(tb.tb(ifreq_check(i), *) GT TB_max OR tb.tb(ifreq_check(i), *) LT TB_min)
  IF jj(0) NE -1 THEN BEGIN
   flag_1b(jj) = flag_1b(jj) OR 512
   flag_2a(jj) = flag_2a(jj) OR 512
   flag_T2b(jj) = flag_Tb(jj) OR 512   
   flag_q2b(jj) = flag_q2b(jj) OR 512   
  ENDIF

  IF l1c.el(0) NE -999. THEN BEGIN 
   FOR j = 0, N_ELEMENTS(l1c.el)-1 DO BEGIN

    jj = WHERE(l1c.tb(ifreq_check(i), j, *) GT TB_max OR l1c.tb(ifreq_check(i), j, *) LT TB_min) 
    IF jj(0) NE -1 THEN BEGIN
     flag_1c(jj) = flag_1b(jj) OR 512
     flag_2c(jj) = flag_2c(jj) OR 512  
    ENDIF
   
   ENDFOR
  ENDIF
 ENDFOR

ENDIF

;**ALL radiometers - check threshold defined in par_mwr_pro.pro
;Bit11: retrieved LWP/IWV threshold (set in par_mwr_pro.pro)
IF N_ELEMENTS(iwv) GT 0 AND N_ELEMENTS(algo_iwv) NE 0 THEN BEGIN

 IF verbose THEN print, 'PL_MK: checking IWV threshold'
 ii = WHERE(iwv*SIN(el*!DTOR) GT iwv_max OR iwv*SIN(el*!DTOR) LT iwv_min AND iwv NE -999.)
 IF ii(0) NE -1 THEN flag_2a(ii) = flag_2a(ii) OR 1024
ENDIF

IF N_ELEMENTS(lwp) GT 0 AND N_ELEMENTS(algo_lwp) NE 0 THEN BEGIN

 IF verbose THEN print, 'PL_MK: checking LWP threshold'
 ii = WHERE(lwp*SIN(el*!DTOR) GT lwp_max OR lwp*SIN(el*!DTOR) LT lwp_min AND lwp NE -999.)
 IF ii(0) NE -1 THEN flag_2a(ii) = flag_2a(ii) OR 1024
ENDIF
;**ALL radiometers - check threshold defined in par_mwr_pro.pro
;Bit12: retrieved TEMPERATURE threshold (set in par_mwr_pro.pro)

IF N_ELEMENTS(T2b) GT 0 THEN BEGIN

 IF verbose THEN print, 'PL_MK: checking T threshold (2b)'
 nz = N_ELEMENTS(T2b(*, 0))
 FOR i = 0, nz-1 DO BEGIN
  ii = WHERE((T2b(i, *) GT T_max OR T2b(i, *) LT T_min) AND T2b(i, *) NE -999.)
  IF ii(0) NE -1 THEN flag_T2b(ii) = flag_T2b(ii) OR 2048
 ENDFOR

ENDIF

IF N_ELEMENTS(q2b) GT 0 THEN BEGIN

 IF verbose THEN print, 'PL_MK: checking q threshold (2b)'
 nz = N_ELEMENTS(q2b(*, 0))
 FOR i = 0, nz-1 DO BEGIN
  ii = WHERE((q2b(i, *) GT q_max OR q2b(i, *) LT q_min) AND q2b(i, *) NE -999.)
  IF ii(0) NE -1 THEN flag_q2b(ii) = flag_q2b(ii) OR 2048
 ENDFOR

ENDIF

IF N_ELEMENTS(T2c) GT 0 THEN BEGIN

 IF verbose THEN print, 'PL_MK: checking T threshold (2c)'
 nz = N_ELEMENTS(par2c.T(*, 0))
 FOR i = 0, nz-1 DO BEGIN
  ii = WHERE((par2c.T(i, *) GT T_max OR par2c.T(i, *) LT T_min) AND par2c.T(i, *) NE -999.)
  IF ii(0) NE -1 THEN flag_2c(ii) = flag_2c(ii) OR 2048
 ENDFOR

ENDIF

;***lwp offset correction (lwp_cor variable) - only carried out when flag_2a is zero
n_lwp = N_ELEMENTS(lwp)
lwp_cor = REPLICATE(-999., n_lwp)

IF set_lwp_off EQ 1 AND N_ELEMENTS(lwp) GT 0 AND N_ELEMENTS(algo_lwp) NE 0 THEN BEGIN

 OFFSET_LWP, mwr_dir, date, tb.time, tb.el, lwp, lwp_std_thres, flag_2a, lwp_cor, verbose = 0

ENDIF ; set_lwp_off eq 1

;***STORE DATA as netcdf files
;* Data are only archived as netcdf up to "date of last change" specified in filter_*.dat or
;* if offc (override filter function condition) is set to 1 in parameter file

IF N_ELEMENTS(offc) EQ 0 THEN offc = 0

IF (LONG(dolc) GE LONG(date)) OR (offc EQ 1) THEN BEGIN

;****WRITE LEVEL1 DATA

;*************
;WRITE_LEVEL1B 
;*************
;**put met data (if available) on tb time
 n_tb = LONG(N_ELEMENTS(tb.time))
 tb_t = REPLICATE(-999., n_tb)
 tb_p = REPLICATE(-999., n_tb)
 tb_rh = REPLICATE(-999., n_tb)

 IF N_ELEMENTS(met.time) GT 1 THEN BEGIN 
 
  FOR i = 0l, n_tb-1l, 100 DO BEGIN
   CLOSEST, tb.time(i), met.time, dmin, ind
   IF dmin LT 0.25 THEN BEGIN
    tb_t(i) = met.t(ind)
    tb_p(i) = met.p(ind)
    tb_rh(i) = met.q(ind)
   ENDIF
  ENDFOR

 ENDIF

;**put irt data (if available) on tb time
 tb_irt = REPLICATE(-999., 1, n_tb)
 el_irt = REPLICATE(-999., n_tb)

 IF N_ELEMENTS(irt.time) GT 1 THEN BEGIN 

  n_wavel_irt = N_ELEMENTS(irt.wavel)
  tb_irt = REPLICATE(-999., n_wavel_irt, n_tb)
  el_irt = REPLICATE(-999., n_tb)
 
  print, 'PL_MK: find matching IRT TBs to TB times'
  FOR i = 0l, n_tb-1l DO BEGIN
   CLOSEST, tb.time(i), irt.time, dmin, ind
   IF dmin LT 10./3600. THEN BEGIN
    tb_irt(*, i) = irt.tb(*, ind)
    el_irt(i) = irt.el(ind)
   ENDIF
  ENDFOR

 ENDIF

 IF verbose THEN BEGIN 
  print, 'PL_MK: writing TB-level1b data of ', date
  print, 'to ', data_path1
 ENDIF

;**get tb uncertainty estimates
 path_cal = mwr_dir+'/data/uncertainty/'

 RESTORE, file = path_cal + tb_bias_file
 RESTORE, file = path_cal + tb_cov_file

 n_f = N_ELEMENTS(tb.f)
 IF N_ELEMENTS(tb_bias) NE n_f OR N_ELEMENTS(tb_cov(*, 0)) NE n_f OR N_ELEMENTS(tb_cov(0, *)) NE n_f THEN BEGIN
  print, 'PL_MK: Dimensions of specified TB uncertainties do not agree with measured TBs!'
  print, 'PL_MK: ABORT'
  stop
 ENDIF

;**estimate integration time of measurements through histogram analysis
 t_diff_tb = FLTARR(n_tb-1)
 FOR i = 0l, n_tb-2 DO t_diff_tb(i) = tb.time(i+1)-tb.time(i)
 x = HISTOGRAM(t_diff_tb*3600., binsize=0.1, min = 0.05, max = 60., locations=l)
 i_xmax = WHERE(x EQ MAX(x))
 int_time_s = (l(i_xmax)+l(i_xmax+1))/2d
 int_time = int_time_s/3600d

;**create time_bnds array
 time_bnds = DBLARR(n_tb, 2)
 time_bnds(*, 0) = tb.time
 time_bnds(*, 1) = tb.time + int_time(0)
 time_bnds_o = time_bnds 

;***transform relative humdity go from % to [1]
 iirh = WHERE(tb_rh NE -999.)
 IF iirh(0) NE -1 THEN tb_rh(iirh) = tb_rh(iirh)/100.

 WRITE_LEVEL1B, data_path1, date, mwr_meas, latitude, longitude, altitude, tb.f, tb.time, time_bnds_o, tb_l1b_oc, l1b_oc, offset_tb_type, freq_shift,$
                flag_1b, tb_t, tb_p, tb_rh, az, el, thresholds(0:1), irt.wavel, tb_irt, el_irt, file_naming_l1b, global_atts, tb_bias, tb_cov 

;*************
;WRITE_LEVEL1C
;*************
;**put met data (if available) on l1c time

 IF N_ELEMENTS(l1c.time) GT 1 THEN BEGIN

  n_l1c = N_ELEMENTS(l1c.time)
  l1c_t = REPLICATE(-999., n_l1c)
  l1c_p = REPLICATE(-999., n_l1c)
  l1c_rh = REPLICATE(-999., n_l1c)

  IF N_ELEMENTS(met.time) GT 1 THEN BEGIN 
 
   FOR i = 0, n_l1c-1 DO BEGIN
    CLOSEST, l1c.time(i), met.time, dmin, ind
    IF dmin LT 0.25 THEN BEGIN
     l1c_t(i) = met.t(ind)
     l1c_p(i) = met.p(ind)
     l1c_rh(i) = met.q(ind)
    ENDIF
   ENDFOR

  ENDIF

  IF verbose THEN BEGIN 
   print, 'PL_MK: writing TB-level1c data of ', date
   print, 'to ', data_path1
  ENDIF

;**create time_bnds_1c array
  IF N_ELEMENTS(int_time_1c) LE 0. OR int_time_1c LE 0. THEN BEGIN
   print, 'PL_MK: You must specifiy integration time for BLB files.'
   print, 'ABORTING'
   stop
  ENDIF

  time_bnds_1c = DBLARR(N_ELEMENTS(l1c.time), 2)
  time_bnds_1c(*, 0) = l1c.time
  time_bnds_1c(*, 1) = l1c.time + int_time_1c/3600.
  time_bnds_1c_o = time_bnds_1c

;***transform relative humdity go from % to [1]
  iirh = WHERE(l1c_rh NE -999.)
  IF iirh(0) NE -1 THEN l1c_rh(iirh) = l1c_rh(iirh)/100.

  WRITE_LEVEL1C, data_path1, date, mwr_meas, latitude, longitude, altitude, l1c.f, l1c.time, time_bnds_1c_o, tb_l1c_oc, l1c_oc, offset_tb_type, freq_shift,$
                 l1c.el, l1c.az, flag_1c, l1c_t, l1c_p, l1c_rh, thresholds(0:1), file_naming_l1c, global_atts, tb_bias, tb_cov

 ENDIF

;****WRITE LEVEL 2 DATA

;****rets: will be stored in retrieval nc files
IF N_ELEMENTS(algo_iwv) EQ 0 THEN BEGIN
 algo_iwv_c = '-'
ENDIF ELSE BEGIN
 algo_iwv_c = algo_iwv
ENDELSE
IF N_ELEMENTS(algo_lwp) EQ 0 THEN BEGIN
 algo_lwp_c = '-'
ENDIF ELSE BEGIN
 algo_lwp_c = algo_lwp
ENDELSE
IF N_ELEMENTS(algo_wdl) EQ 0 THEN BEGIN
 algo_wdl_c = '-'
ENDIF ELSE BEGIN
 algo_wdl_c = algo_wdl
ENDELSE
IF N_ELEMENTS(algo_taa) EQ 0 THEN BEGIN
 algo_taa_c = '-'
ENDIF ELSE BEGIN
 algo_taa_c = algo_taa
ENDELSE
IF N_ELEMENTS(algo_tze) EQ 0 THEN BEGIN
 algo_tze_c = '-'
ENDIF ELSE BEGIN
 algo_tze_c = STRMID(algo_tze(0), 0, 8)
ENDELSE
IF N_ELEMENTS(algo_hze) EQ 0 THEN BEGIN
 algo_hze_c = '-'
ENDIF ELSE BEGIN
 algo_hze_c = STRMID(algo_hze(0), 0, 8)
ENDELSE
IF N_ELEMENTS(algo_tel) EQ 0 THEN BEGIN
 algo_tel_c = '-'
ENDIF ELSE BEGIN
 algo_tel_c = algo_tel
ENDELSE

rets = [algo_iwv_c, algo_lwp_c, algo_wdl_c, algo_taa_c, algo_tze_c, algo_hze_c, algo_tel_c]

;*************
;WRITE_LEVEL2A
;- the option for retrieving atmospheric attenuation at arbitrary frequencies must still be coded!
;*************
 IF N_ELEMENTS(algo_iwv) NE 0 OR N_ELEMENTS(algo_lwp) NE 0 OR N_ELEMENTS(algo_wdl) NE 0 OR N_ELEMENTS(algo_taa) NE 0 THEN BEGIN  

  IF verbose THEN BEGIN 
   print, 'PL_MK: writing level2a data of ', date
   print, 'to ', data_path2
  ENDIF

  IF N_ELEMENTS(lwp) NE 0 THEN lwpo = lwp
  time_bnds_o = time_bnds 


  WRITE_LEVEL2A, data_path2, date, mwr_meas, latitude, longitude, altitude, tb.time, time_bnds_o, flag_2a,$
                 tb_t, tb_p, tb_rh, az, el, thresholds(0:5), file_naming_l1b, file_naming_2a_iwv, file_naming_2a_lwp,$
                 file_naming_2a_wdl, file_naming_2a_taa, global_atts, rets, ang_ret_2a,$
                 iwv=iwv, diwv=diwv, err_iwv=err_iwv, lwpo=lwpo, lwp_cor=lwp_cor, dlwp=dlwp, err_lwp=err_lwp, wdl=wdl, dwdl=dwdl,$
                 err_wdl=err_wdl, taa_a=taa_a, dtaa=dtaa, err_taa=err_taa, f_taa=f_taa, dep_l2a_iwv=dep_l2a_iwv, dep_l2a_lwp=dep_l2a_lwp   
 ENDIF


;*************
;WRITE_LEVEL2B
;*************
 IF N_ELEMENTS(algo_tze) NE 0 OR N_ELEMENTS(algo_hze) NE 0 THEN BEGIN

  IF verbose THEN BEGIN 
   print, 'PL_MK: writing level2b data of ', date
   print, 'to ', data_path2
  ENDIF

  time_bnds_o = time_bnds

  WRITE_LEVEL2B, data_path2, date, mwr_meas, latitude, longitude, altitude, tb.time, time_bnds_o, flag_T2b, flag_q2b,$
                 tb_t, tb_p, tb_rh, az, el, thresholds(0:9), file_naming_l1b, file_naming_2b_tze, file_naming_2b_hze, global_atts, rets,$
                 z_final, T2b=T2b, dT2b=dT2b, err_T2b=err_T2b, ang_ret_T2b=ang_ret_T2b, q2b=q2b, dq2b=dq2b, err_q2b=err_q2b, ang_ret_q2b=ang_ret_q2b,$
                 dep_l2b_t=dep_l2b_t, dep_l2b_q=dep_l2b_q
 ENDIF

;*************
;WRITE_LEVEL2C
;*************

 IF N_ELEMENTS(algo_tel) AND l1c.time(0) NE -999. THEN BEGIN

  IF verbose THEN BEGIN 
   print, 'PL_MK: writing level2c data of ', date
   print, 'to ', data_path2
  ENDIF

  time_bnds_1c_o = time_bnds_1c

  WRITE_LEVEL2C, data_path2, date, mwr_meas, latitude, longitude, altitude, l1c.time, time_bnds_1c_o, flag_2c, l1c_t,$
                 l1c_p, l1c_rh, l1c.az, l1c.el, thresholds(0:7), file_naming_l1c, file_naming_2c_tel, global_atts, rets,$
                 z_final, T2c=T2c, dT2c=dT2c, err_T2c=err_T2c, dep_l2c_t=dep_l2c_t
 ENDIF


ENDIF ELSE BEGIN

 IF verbose THEN print, 'PL_MK: no netcdf files written'  

ENDELSE

;***plot level1 TB time series
n_f = N_ELEMENTS(tb.f)
n_plot_times = N_ELEMENTS(plot.time)
plot_start = plot.time(0)
plot_end = plot.time(N_ELEMENTS(plot.time)-1)

IF n_plot_times LE 2 THEN BEGIN
 print, 'PL_MK: NO PLOT CREATED!'
 GOTO, SKIP_PLOT
ENDIF

IF verbose THEN print, 'creating level1 TB-plot ...'
n_ff = n_f

IF n_f GT 7 THEN BEGIN
 x_min = [0.17, 0.57, REPLICATE(0.17, CEIL(FLOAT(n_f)/2.)), REPLICATE(0.57, (n_f)/2)]
 x_max = [0.5, 0.9, REPLICATE(0.5, CEIL(FLOAT(n_f)/2.)), REPLICATE(0.9, (n_f)/2)]

 y_min = [0.85, 0.85, 0.8 - ((0.7/(FLOAT(CEIL(FLOAT(n_f)/2.))))+(FINDGEN(CEIL(FLOAT(n_f)/2.))*(0.7/(FLOAT(CEIL(FLOAT(n_f)/2.)))))),$
          0.8 - ((0.7/(FLOAT(CEIL(FLOAT(n_f)/2.))))+(FINDGEN(FLOAT(n_f)/2.))*(0.7/(FLOAT(CEIL(FLOAT(n_f)/2.)))))]
 y_max = [0.95, 0.95, 0.8 - (FINDGEN(CEIL(FLOAT(n_f)/2.))*(0.7/(FLOAT(CEIL(FLOAT(n_f)/2.))))),$
          0.8 - (FINDGEN(FLOAT(n_f)/2.))*(0.7/(FLOAT(CEIL(FLOAT(n_f)/2.))))]

ENDIF ELSE BEGIN

 IF N_ELEMENTS(name_ap) NE 0 THEN BEGIN
  n_ff = n_f + 1

  IF N_ELEMENTS(ch2) EQ 0 THEN BEGIN
   print, 'PL_MK: need to specify ch2 in par_mwr.pro'
   stop
  ENDIF
  IF N_ELEMENTS(ch1) EQ 0 THEN BEGIN
   print, 'PL_MK: need to specify ch1 in par_mwr.pro'
   stop
  ENDIF

 ENDIF

 x_min = REPLICATE(0.17, n_ff+2)
 x_max = REPLICATE(0.9, n_ff+2)
 y_min = [0.85, 0.75, 0.7-((0.6/FLOAT(n_ff))+FINDGEN(n_ff)*(0.6/FLOAT(n_ff)))]
 y_max = [0.95, 0.85, 0.7-FINDGEN(n_ff)*(0.6/FLOAT(n_ff))]

ENDELSE

!P.font = 0
IF n_f LE 7 THEN !p.multi = [0, 1, n_ff+2]
IF n_f GT 7 THEN !p.multi = [0, 2, CEIL(FLOAT(n_ff+1.)/2.)]

!P.Charsize=1.3
!P.thick = 1.5
!X.minor = 6
SET_PLOT, 'ps'
LOADCOL, home_path+'/mwr_pro/source/col1'
DEVICE, /color, /landscape, filename = plot_path1+date+'_'+mwr_meas+'_l1.ps'

;***plotting for MWRs with less than 8 channels
IF n_f LE 7 THEN BEGIN

 a1 = STRING(ang_low, format = '(f5.1)')
 a2 = STRING(ang_high, format = '(f5.1)')
 ang_range = ', Elev. = ['+a1+', '+a2+']'

;**plot environmental T & RH values if available
 IF N_ELEMENTS(met.time) GT 1 THEN BEGIN
  tit = ''
  PLOT, met.time, met.q, yrange=[0. ,100.], ystyle=1, title = tit,$
  ytickv=1, ytitle='RH surf. [%]', ycharsize = 0.9, xticklen = 0.07,$
  xrange=[0.,24.], xstyle=1, xtitle='Time (UTC) [h]', xcharsize = 1e-5,xtickinterval=6,$
  position=[x_min(0), y_min(0), x_max(0), y_max(0)], xthick = 1.5, ythick = 1.5, /NODATA

  OPLOT, met.time, met.q, psym = 3

  dy = MAX(met.t) - MIN(met.t)
  ymin = MIN(met.t) - dy/10.
  ymax = MAX(met.t) + dy/10.

  PLOT, met.time, met.t, yrange=[ymin, ymax], ystyle=1,$
  ytickv=1, ytitle='T surf. [K]', ycharsize = 0.9, xticklen = 0.07,xtickinterval=6,$
  xrange=[0.,24.], xstyle=1, xtitle='Time (UTC) [h]', xthick = 1.5, ythick = 1.5,$
  position=[x_min(1), y_min(1), x_max(1), y_max(1)], xcharsize = 0.8, /NODATA

  OPLOT, met.time, met.t, psym = 3
 ENDIF ; end plotting of env. variables

;**brightness temperature plots  
 FOR i = 0, n_ff-1 DO BEGIN

  IF i LT n_f THEN BEGIN
   plot_tb = tb_p_oc(i, *)
   name = STRING(tb.f(i), FORMAT='(F6.2)') + ' GHz'
  ENDIF ELSE BEGIN
   plot_tb = plot.tb(ch2, *) - plot.tb(ch1, *)
   name = name_ap
  ENDELSE
 
  ii = i+2
  dy = MAX(plot_tb) - MIN(plot_tb)
  ymin = MIN(plot_tb) - dy/10.
  ymax = MAX(plot_tb) + dy/10.

;*make plot frames
  IF y_min(ii) EQ MIN(y_min) THEN BEGIN
   PLOT, plot.time, plot_tb, yrange=[ymin, ymax], ystyle = 1,$
   xrange = [0., 24.], xstyle=1, xtitle = 'Time (UTC) [h]',xtickinterval=6,$
   xcharsize = 1.5, ytickv = 2, xthick = 1.5, xticklen = -0.05,$
   position=[x_min(ii), y_min(ii), x_max(ii), y_max(ii)], /NODATA
  ENDIF ELSE BEGIN
   PLOT, plot.time, plot_tb, yrange=[ymin, ymax], ystyle = 1,$
   xrange = [0., 24.], xstyle=1, xtitle = 'Time (UTC) [h]',xtickinterval=6,$
   xcharsize = 1e-5, ytickv = 2,xticklen = -0.05, xthick = 1.5,$
   position=[x_min(ii), y_min(ii), x_max(ii), y_max(ii)], /NODATA
  ENDELSE

;*plot brightness temperatures
  OPLOT, plot.time, plot_tb, psym = 3
  
;*calc and plot mean & stddev
  m = MEAN(plot_tb)
  str = STRING(m, FORMAT='(F6.2)') + ' +/- '+$
        STRING(STDDEV(plot_tb),FORMAT='(F5.2)') + 'K'
  XYOUTS, x_min(ii)+0.6, y_max(ii)-0.02, str, /NORMAL, size=0.7, color=2

  y=[m, m]
  x =[plot_start, plot_end] 
  OPLOT, x, y, linestyle = 1

;*plot channel number in GHz
  XYOUTS, x_min(ii)+0.02, y_max(ii)-0.02, name, /normal, charsize=0.7, color=2

;*overplot flagged times

  index = WHERE(flag_1b GT 0, count)
  IF count GT 0 THEN BEGIN
   x = tb.time(index)
   y = REPLICATE(ymin, count)
   OPLOT, x, y, psym=4, color=7, thick=2, symsize = 0.7
  ENDIF
  
 ENDFOR ; loop over channels

 a1 = STRING(ang_low, format = '(f5.1)')
 a2 = STRING(ang_high, format = '(f5.1)')
 ang_range = ', Elev. = ['+a1+', '+a2+']'

 IF N_ELEMENTS(met.time) GT 1 THEN BEGIN
  XYOUTS, 0.532, 0.985, 'Measurement: '+mwr_meas + ', '+ date + ang_range, alignment=0.5, charsize=1.2, /NORMAL
 ENDIF ELSE BEGIN
  XYOUTS, 0.532, 0.785, 'Measurement: '+mwr_meas + ', '+ date + ang_range, alignment=0.5, charsize=1.2, /NORMAL
 ENDELSE
 IF count GT 1 THEN  XYOUTS, 0.532, 0.01, 'flagged data', color=7, charsize=1.0, alignment=0.5, /NORMAL
 XYOUTS, 0.1, y_min(n_ff+1)+(y_max(2)-y_min(n_ff+1))/2., 'Brightness Temperature [K]', alignment=0.5, orientation=90., /NORMAL

ENDIF ELSE BEGIN ; n_f LE 7, now n_f GT 7

;**plot environmental T & RH values if available
 IF N_ELEMENTS(met.time) GT 1 THEN BEGIN
  PLOT, met.time, met.q, yrange=[0. ,100.], ystyle=1,xtickinterval=6,$
  ytickv=1, ytitle='RH env. [%]', ycharsize = 0.9, xticklen=0.07,$
  xrange=[0.,24.], xstyle=1, xtitle='Time (UTC) [h]',$
  position=[x_min(0), y_min(0), x_max(0), y_max(0)], /NODATA

  OPLOT, met.time, met.q, psym = 3

  dy = MAX(met.t) - MIN(met.t)
  ymin = MIN(met.t) - dy/10.
  ymax = MAX(met.t) + dy/10.

  PLOT, met.time, met.t, yrange=[ymin, ymax], ystyle=1,xtickinterval=6,$
  ytickv=1, ytitle='T env. [K]', ycharsize = 0.9, xticklen=0.07,$
  xrange=[0.,24.], xstyle=1, xtitle='Time (UTC) [h]',$
  position=[x_min(1), y_min(1), x_max(1), y_max(1)], /NODATA

  OPLOT, met.time, met.t, psym = 3
 ENDIF ; end plotting of env. variables

;**brightness temperature plots  
 FOR i = 0, n_f-1 DO BEGIN

  ii = i+2
  dy = MAX(tb_p_oc(i, *)) - MIN(tb_p_oc(i, *))
  ymin = MIN(tb_p_oc(i, *)) - dy/10.
  ymax = MAX(tb_p_oc(i, *)) + dy/10.

;*make plot frames
  IF y_min(ii) EQ MIN(y_min) OR i EQ n_f-1 THEN BEGIN
   PLOT, plot.time, tb_p_oc(i, *), yrange=[ymin, ymax], ystyle = 1,$
   xrange = [0., 24.], xstyle=1, xtitle = 'Time (UTC) [h]',xtickinterval=6,$
   xcharsize = 1.5, ytickv = 2, xthick = 1.5, xticklen = -0.05,$
   position=[x_min(ii), y_min(ii), x_max(ii), y_max(ii)], /NODATA
  ENDIF ELSE BEGIN
   PLOT, plot.time, tb_p_oc(i, *), yrange=[ymin, ymax], ystyle = 1,$
   xrange = [0., 24.], xstyle=1, xtitle = 'Time (UTC) [h]',xtickinterval=6,$
   xcharsize = 1e-5, ytickv = 2,xticklen = -0.05, xthick = 1.5,$
   position=[x_min(ii), y_min(ii), x_max(ii), y_max(ii)], /NODATA
  ENDELSE

;*plot brightness temperatures
  OPLOT, plot.time, tb_p_oc(i, *), psym = 3
  
;*calc and plot mean & stddev

  m = MEAN(tb_p_oc(i, *))
  str = STRING(m, FORMAT='(F6.2)') + ' +/- '+$
        STRING(STDDEV(tb_p_oc(i, *)),FORMAT='(F5.2)') + 'K'
  XYOUTS, x_min(ii)+0.2, y_max(ii)-0.02, str, /NORMAL, size=0.7, color=2

  y=[m, m]
  x =[plot_start, plot_end] 
  OPLOT, x, y, linestyle = 1

;*plot channel number in GHz
  XYOUTS, x_min(ii)+0.02, y_max(ii)-0.02, STRING(tb.f(i), FORMAT='(F6.2)')+ ' GHz',/normal, charsize=0.7, color=2

;*overplot flagged times
  index = WHERE(flag_1b GT 0, count)
  IF count GT 0 THEN BEGIN
   x = tb.time(index)
   y = REPLICATE(ymin, count)
   OPLOT, x, y, psym=4, color=7, thick=2, symsize = 0.7
  ENDIF

 ENDFOR ; loop over channels

 a1 = STRING(ang_low, format = '(f5.1)')
 a2 = STRING(ang_high, format = '(f5.1)')
 ang_range = ', Elev. = ['+a1+', '+a2+']'

 XYOUTS, 0.532, 0.985, 'Measurement: '+mwr_meas + ', '+ date + ang_range, alignment=0.5, charsize=1.2, /NORMAL
 IF count GT 1 THEN  XYOUTS, 0.532, 0.03, 'flagged data', color=7, charsize=1.0, alignment=0.5, /NORMAL
 XYOUTS, 0.1, y_min(n_f+1)+(y_max(2)-y_min(n_f+1))/2., 'Brightness Temperature [K]', alignment=0.5, orientation=90., /NORMAL

ENDELSE ; n_f GT 7

DEVICE, /CLOSE

;***plot tb spectra
if verbose then print,"creating '",date+'_'+mwr_meas+'_l1_sp.ps','"'
!p.font =0
!P.Charsize=1.2
!P.thick = 1.3
!X.minor = 5

SET_PLOT, 'PS'
;DEVICE, /color, /landscape, filename = plot_path1+date+'_'+mwr_meas+'_l1_sp.ps'
DEVICE, /color, /portrait, filename = plot_path1+date+'_'+mwr_meas+'_l1_sp.ps'

x_min = 0.17
x_max = 0.9
y_min = 0.1  
y_max = 0.8
!P.charsize = 1.8

tb_mean = FLTARR(n_f)
tb_std = FLTARR(n_f)
FOR i= 0, n_f-1 DO BEGIN
 tb_mean(i) = MEAN(tb_p_oc(i, *))
 tb_std(i) = STDDEV(tb_p_oc(i, *))
ENDFOR

i_f_plot = SORT(tb.f)
f_plot = tb.f(i_f_plot)

tit = 'Measurement: '+mwr_meas+' '+date+': TB daily means +/-STDEV' + ang_range
PLOT, f_plot, tb_mean, yrange = [MIN(tb_mean)-MAX(tb_std), MAX(tb_mean)+MAX(tb_std)], ystyle=0,$
      xrange=[f_plot(0), f_plot(n_f-1)], xstyle=0, xtitle='Frequency [GHz]', ytitle = 'Brightness temperature [K]',$
      ytickv=2, position=[x_min, y_min, x_max, y_max], psym=4, xcharsize = 1.5, ycharsize = 1.5, title = tit,$
      xthick = 1.5, ythick = 1.5

OPLOT, f_plot, tb_mean+tb_std, psym = 1
OPLOT, f_plot, tb_mean-tb_std, psym = 1
FOR i = 0, N_ELEMENTS(f_plot)-1 DO OPLOT, [f_plot(i), f_plot(i)], [tb_mean(i)-tb_std(i), tb_mean(i)+tb_std(i)] 

; close ps-file
DEVICE, /CLOSE

;***plot angles & flags
IF verbose then print,"creating '",date+'_'+mwr_meas+'_flag.ps','"'
!p.font = 0
!P.Charsize = 1.8
!P.thick = 1.5
!P.multi = 14
!X.minor = 6

SET_PLOT, 'PS'
DEVICE, /color, /landscape, filename = plot_path1+date+'_'+mwr_meas+'_flag.ps'

x_min1 = REPLICATE(0.17, 12)
x_max1 = REPLICATE(0.9, 12)
x_min2 = REPLICATE(0.17, 2)
x_max2 = REPLICATE(0.9, 2)

start_y1 = 0.05
end_y1 = 0.7
start_y2 = end_y1
end_y2 = 0.95

y_min1 = start_y1 + FINDGEN(12)*(end_y1-start_y1)/12.
y_max1 = start_y1 + (FINDGEN(12)*(end_y1-start_y1)/12.) + (end_y1-start_y1)/12.
y_min2 = start_y2 + (FINDGEN(2)*(end_y2-start_y2)/2.)
y_max2 = start_y2 + (FINDGEN(2)*(end_y2-start_y2)/2.) + (end_y2-start_y2)/2.

;**plot elevation
PLOT, tb.time, tb.el, yrange=[0., 99.], ystyle=1,$
ytickv=1, ytitle='Elevation [deg]', ycharsize = 0.8, xticklen = 0.07,$
xrange=[0.,24.], xstyle=1, xcharsize = 1e-5, xtickinterval = 6,$
position=[x_min2(0), y_min2(0), x_max2(0), y_max2(0)], xthick = 1.5, ythick = 1.5, /NODATA
OPLOT, tb.time, el, color = 3, thick = 1
;**plot azimuth
PLOT, tb.time, az, yrange=[0., 360.], ystyle=1, title = 'Angles and flags: ' + mwr_meas + ', ' + date,$
ytickv=1, ytitle='Azimuth [deg]', ycharsize = 0.8, xticklen = 0.07,$
xrange=[0.,24.], xstyle=1, xcharsize = 1e-5, xtickinterval = 6,$
position=[x_min2(1), y_min2(1), x_max2(1), y_max2(1)], xthick = 1.5, ythick = 1.5, /NODATA
OPLOT, tb.time, az, color = 5, thick = 1
 
;**data flag
X = [-1, 1, 1, -1]
Y = [-1, -1, 1, 1]
USERSYM, X, Y, /FILL

ytit = ['Bit1', 'Bit2', 'Bit3', 'Bit4', 'Bit5', 'Bit6', 'Bit7', 'Bit8', 'Bit9', 'Bit10', 'Bit11', 'Bit12']
bit_char = ['Manual filter Band1', 'Manual filter Band2', 'Manual filter Band3', 'Rain', 'Sanity Receiver1', 'Sanity Receiver2', 'Sun', 'TB thres. Band1',$
            'TB thres. Band2', 'TB thres. Band3', 'LWP/IWV retrieval thres.', 'Temperature retrieval thres.']

FOR i = 0, 10 DO BEGIN ; loop over 12 flags
 PLOT, [0, 0], [24, 1], /NODATA, yrange=[0.9,1.1], ystyle = 1, ytitle = ytit(i),$
       xrange = [0., 24.], xstyle = 1, xcharsize = 0.0000001, ycharsize = 0.9,$
       position = [x_min1(11-i), y_min1(11-i), x_max1(11-i), y_max1(11-i)], yticks = 1, ytickname = [' ', ' '],$
       xtickinterval=6, xticklen = 0.07

 index = WHERE((flag_2a AND 2^i) EQ 2^i, count)
 IF count GT 0 THEN BEGIN
  IF count EQ 1 THEN BEGIN
   count = 2
   index = [index(0), index(0)]
  ENDIF
  h = REPLICATE(1, count)
  OPLOT, tb.time(index), h, psym=8, thick=5, color=2
 ENDIF

 XYOUTS, x_min1(11-i)+(x_max1(11-i)-x_min1(11-i))*0.08, y_max1(11-i)-(y_max1(11-i)-y_min1(11-i))*0.25, bit_char(i), charsize = 0.6, /normal
ENDFOR

i = 11
PLOT, [0, 0], [24, 1], /NODATA, yrange=[0.9,1.1], ystyle = 1, ytitle = ytit(i),$
      xrange = [0., 24.], xstyle = 1, xcharsize = 1., ycharsize=0.9, ytickname = [' ', ' '], xtitle = 'Time [UTC]',$
      position = [x_min1(11-i), y_min1(11-i), x_max1(11-i), y_max1(11-i)], yticks = 1, xticklen = 0.2, xtickinterval=6

index = WHERE((flag_2a AND 2^i) EQ 2^i, count)
IF count GT 0 THEN BEGIN
 IF count EQ 1 THEN BEGIN
  count = 2
  index = [index(0), index(0)]
 ENDIF
 h = REPLICATE(1, count)
 OPLOT, tb.time(index), h, psym=8, thick=5, color=2
ENDIF

XYOUTS, x_min1(11-i)+(x_max1(11-i)-x_min1(11-i))*0.08, y_max1(11-i)-(y_max1(11-i)-y_min1(11-i))*0.25, bit_char(i), charsize = 0.6, /normal

; close ps-file
DEVICE, /CLOSE

;***plot level2a (iwv, lwp, wdl)
algos_2a = REPLICATE('dummy', 3)
data_2a = REPLICATE(-999., 3, N_ELEMENTS(plot.time))
units_2a = [' kgm!e-2!N', ' gm!e-2!N', ' mm']
ytit_2a = ['IWV / kgm!e-2!N', 'LWP / gm!e-2!N', 'Wet delay / mm'] 
IF set_lwp_off EQ 1 THEN ytit_2a = ['IWV / kgm!e-2!N', 'LWP (offset corrected) / gm!e-2!N', 'Wet delay / mm']
i_ang = WHERE(tb.el GT ang_low AND tb.el LT ang_high AND iwv GE 0., n_i_ang)

IF N_ELEMENTS(algo_iwv) EQ 0 AND N_ELEMENTS(algo_lwp) EQ 0 AND N_ELEMENTS(algo_wdl) EQ 0 THEN BEGIN
 print, 'no IWV, LWP, WDL retrieval specified - nothing to plot'
 GOTO, SKIP_PLOT1
ENDIF

IF n_i_ang GT 1 THEN BEGIN

 data_2a = REPLICATE(-999., 3, N_ELEMENTS(plot.time))

 IF N_ELEMENTS(algo_iwv) NE 0 THEN BEGIN
  algos_2a(0) = 'iwv'
  data_2a(0, *) = iwv(i_ang)
 ENDIF

 IF N_ELEMENTS(algo_lwp) NE 0 THEN BEGIN
  algos_2a(1) = 'lwp'
  data_2a(1, *) = lwp(i_ang)*1000.
  IF set_lwp_off EQ 1 THEN data_2a(1, *) = lwp_cor(i_ang)*1000.
 ENDIF
 IF N_ELEMENTS(algo_wdl) NE 0 THEN BEGIN
  algos_2a(2) = 'wdl'
  data_2a(2, *) = wdl(i_ang)
 ENDIF

 ii = WHERE(algos_2a NE 'dummy', n_algos_2a)
 IF n_algos_2a(0) NE -1 THEN BEGIN
  algos_2a = algos_2a(ii)
  data_2a = data_2a(ii, *)
  units_2a = units_2a(ii)
  ytit_2a = ytit_2a(ii)
  user_range_2a_min = user_range_min(ii)
  user_range_2a_max = user_range_max(ii)
 
  IF verbose then print,"creating '",date+'_'+mwr_meas+'_l2a.ps','"'
  !p.font = 0
  !P.Charsize = 1.8
  !P.thick = 1.3
  !X.minor = 6

  SET_PLOT, 'PS'
  DEVICE, /color, /landscape, filename = plot_path2+date+'_'+mwr_meas+'_l2a.ps'

  x_min = REPLICATE(0.17, n_algos_2a+2)
  x_max = REPLICATE(0.9, n_algos_2a+2)
  y_min = [0.85, 0.75, 0.7-((0.6/FLOAT(n_algos_2a))+FINDGEN(n_algos_2a)*(0.6/FLOAT(n_algos_2a)))]
  y_max = [0.95, 0.85, 0.7-FINDGEN(n_algos_2a)*(0.6/FLOAT(n_algos_2a))]

;**plot environmental T & RH values if available
  IF N_ELEMENTS(met.time) GT 1 THEN BEGIN
   tit = 'Measurement: '+mwr_meas + ', '+ date + ang_range
   PLOT, met.time, met.q, yrange=[0. ,100.], ystyle=1, title = tit,$
   ytitle='RH env. [%]', ycharsize = 0.9, xticklen = 0.07, yticks=4,$
   xrange=[0.,24.], xstyle=1, xtitle='Time (UTC) [h]', xcharsize = 1e-5,xtickinterval=6,$
   position=[x_min(0), y_min(0), x_max(0), y_max(0)], xthick = 1.5, ythick = 1.5, /NODATA

   OPLOT, met.time, met.q, psym = 3

   dy = MAX(met.t) - MIN(met.t)
   ymin = MIN(met.t) - dy/10.
   ymax = MAX(met.t) + dy/10.

   PLOT, met.time, met.t, yrange=[ymin, ymax], ystyle=1,xtickinterval=6,$
   ytitle='T env. [K]', ycharsize = 0.9, xticklen = 0.07, yticks=4,$
   xrange=[0.,24.], xstyle=1, xtitle='Time (UTC) [h]', xthick = 1.5, ythick = 1.5,$
   position=[x_min(1), y_min(1), x_max(1), y_max(1)], xcharsize = 0.8, /NODATA

   OPLOT, met.time, met.t, psym = 3
  ENDIF ; end plotting of env. variables

  FOR i = 0, n_algos_2a-1 DO BEGIN

   ii = i+2
   index = WHERE(flag_2a(i_ang) EQ 0, count)
   IF plot_range EQ 'flag' AND count GT 1 THEN BEGIN
    dy = MAX(data_2a(i, index)) - MIN(data_2a(i, index))
    ymin = MIN(data_2a(i, index)) - dy/10.
    ymax = MAX(data_2a(i, index)) + dy/10.
   ENDIF ELSE BEGIN
    dy = MAX(data_2a(i, *)) - MIN(data_2a(i, *))
    ymin = MIN(data_2a(i, *)) - dy/10.
    ymax = MAX(data_2a(i, *)) + dy/10.
   ENDELSE

   IF user_range_2a_min(i) NE -999. THEN ymin = user_range_2a_min(i)
   IF user_range_2a_max(i) NE -999. THEN ymax = user_range_2a_max(i)

;*make plot frames
   xtl = -0.05
   IF y_min(ii) EQ MIN(y_min) OR i EQ n_algos_2a-1 THEN BEGIN
    tit = ''
    IF i EQ 0 THEN BEGIN
     tit = 'Measurement: '+mwr_meas + ', '+ date + ang_range
     xtl = 0.02
    ENDIF

    PLOT, plot.time, data_2a(i, *), yrange=[ymin, ymax], ystyle = 1,xtickinterval=6,$
    xrange = [0., 24.], xstyle=1, xtitle = 'Time (UTC) [h]', ytitle = ytit_2a(i),$
    xcharsize = 1.5, ytickv = 2, xthick = 1.5, xticklen = xtl, ythick= 1.5,$
    position=[x_min(ii), y_min(ii), x_max(ii), y_max(ii)], title = tit, /NODATA, /NOERASE
   ENDIF ELSE BEGIN
    tit = ''
    IF i EQ 0 THEN BEGIN
     IF N_ELEMENTS(met.time) EQ 0 THEN tit = 'Measurement: '+mwr_meas + ', '+ date + ang_range
     xtl = 0.02
    ENDIF
    PLOT, plot.time, data_2a(i, *), yrange=[ymin, ymax], ystyle = 1,xtickinterval=6,$
    xrange = [0., 24.], xstyle=1, xtitle = 'Time (UTC) [h]', ytitle = ytit_2a(i),$
    xcharsize = 1e-5, ytickv = 2, xticklen = xtl, xthick = 1.5, ythick= 1.5,$
    position=[x_min(ii), y_min(ii), x_max(ii), y_max(ii)], title = tit, /NODATA, /NOERASE
   ENDELSE

;*plot level2a data
   IF plot_range EQ 'flag' AND count GT 1 THEN BEGIN
    OPLOT, plot.time, data_2a(i, index), psym = 3
   ENDIF ELSE IF plot_range EQ 'all' THEN BEGIN
    OPLOT, plot.time, data_2a(i, *), psym = 3
   ENDIF
   OPLOT, [plot_start, plot_end], [0., 0.], linestyle = 1

;*calc and plot mean & stddev

   m = MEAN(data_2a(i, *))
   str = STRING(m, FORMAT='(F6.2)') + ' +/- '+$
         STRING(STDDEV(data_2a(i, *)),FORMAT='(F6.2)') + units_2a(i)
   XYOUTS, x_min(ii)+0.2, y_max(ii)-0.02, str, /NORMAL, size=0.7, color=2

   y=[m, m]
   x =[plot_start, plot_end] 
   OPLOT, x, y, linestyle = 2, color = 2

;*overplot flagged times
   index = WHERE(flag_2a GT 0, count)
   IF count GT 1 THEN BEGIN
    x = tb.time(index)
    y = REPLICATE(ymin, count)
    OPLOT, x, y, psym=4, color=7, thick=2, symsize = 0.7
    XYOUTS, 0.532, -0.02, 'flagged data', color=7, charsize=1.0, alignment=0.5, /NORMAL
   ENDIF

  ENDFOR ; loop over level2a products
  
 ENDIF ; n_algos

ENDIF ; i_ang --> ploting range

SKIP_PLOT1:
DEVICE, /CLOSE

;***plot level2a azimuth-time contours @ one or two elevation angles 
;iwv:
IF N_ELEMENTS(aztp_iwv) EQ 0 THEN aztp_iwv = 0
IF N_ELEMENTS(algo_iwv) GT 0 AND aztp_iwv EQ 1 THEN BEGIN
 filename = plot_path2+date+'_'+mwr_meas+'_l2a_iwv_aztp.ps'
 ytit = 'IWV / kgm!e-2!N'

 AZTP, par='IWV', filename=filename, date=date, mwr_meas=mwr_meas, time=tb.time, data=iwv, ele=tb.el, azi=az, ytit=ytit, el_aztp=el_aztp,$
       az_scan_thres=az_scan_thres, color_tab=5, flag_2a=flag_2a, verbose=1  
ENDIF

;lwp:
IF N_ELEMENTS(aztp_lwp) EQ 0 THEN aztp_lwp = 0
IF N_ELEMENTS(algo_lwp) GT 0 AND aztp_lwp EQ 1 THEN BEGIN
 filename = plot_path2+date+'_'+mwr_meas+'_l2a_lwp_aztp.ps'
 ytit = 'LWP / gm!e-2!N' 
 lwpx = lwp
 IF set_lwp_off EQ 1 THEN BEGIN
  ytit = 'LWP / gm!e-2!N'
  lwpx = lwp
 ENDIF
 AZTP, par='LWP', filename=filename, date=date, mwr_meas=mwr_meas, time=tb.time, data=lwpx*1000., ele=tb.el, azi=az, ytit=ytit, el_aztp=el_aztp,$
       az_scan_thres=az_scan_thres, color_tab=1, flag_2a=flag_2a, verbose=1  
ENDIF

;***plot level2b (T profile, q profile) data
i_ang = WHERE(tb.el GT ang_low AND tb.el LT ang_high, n_i_ang)

algos_2b = REPLICATE('dummy', 2)
datas_2b = REPLICATE(-999., 2, N_ELEMENTS(z_final), N_ELEMENTS(plot.time))
units_2b = [' K', ' gm!e-3!N']
ytits_2b = ['Temperature', 'Abs. humidity'] 

IF N_ELEMENTS(algo_tze) EQ 0 AND N_ELEMENTS(algo_hze) EQ 0 THEN BEGIN
 print, 'no TZE or HZE retrieval specified - nothing to plot'
 GOTO, SKIP_PLOT2
ENDIF

IF n_i_ang GT 1 THEN BEGIN

 datas_2b = REPLICATE(-999., 2, N_ELEMENTS(z_final), N_ELEMENTS(plot.time))

 IF N_ELEMENTS(algo_tze) NE 0 THEN BEGIN
  algos_2b(0) = 'tze'
  datas_2b(0, *, *) = T2b(*, i_ang)
 ENDIF
 IF N_ELEMENTS(algo_hze) NE 0 THEN BEGIN
  algos_2b(1) = 'hze'
  datas_2b(1, *, *) = q2b(*, i_ang)*1000.
 ENDIF

 ii = WHERE(algos_2b NE 'dummy', n_algos_2b)
 IF n_algos_2b(0) NE -1 THEN BEGIN
  algo_2b = algos_2b(ii)
  data_2b = REFORM(datas_2b(ii, *, *))
  unit_2b = units_2b(ii)
  ytit_2b = ytits_2b(ii)

  FOR ii = 0, n_algos_2b-1 DO BEGIN
 
   IF verbose then print, 'creating '+date+'_'+mwr_meas+'_l2b_'+algo_2b(ii)+'.ps'
   SET_PLOT, 'ps'
   LOADCT,39
   !X.minor = 6

   DEVICE, /color, /landscape, filename = plot_path2+date+'_'+mwr_meas+'_l2b_'+algo_2b(ii)+'.ps'

   IF algo_2b(ii) EQ 'tze' THEN BEGIN
    index = WHERE(flag_T2b(i_ang) EQ 0, count)
    flagx = flag_T2b(i_ang)
   ENDIF
   IF algo_2b(ii) EQ 'hze' THEN BEGIN
    index = WHERE(flag_q2b(i_ang) EQ 0, count)
    flagx = flag_q2b(i_ang)
   ENDIF
   index_h = WHERE(z_final LE h_max)

   IF plot_range EQ 'flag' AND count GT 1 THEN BEGIN
    dy = MAX(data_2b(ii, index_h, index)) - MIN(data_2b(ii, index_h, index))
    min_p = MIN(data_2b(ii, index_h, index)) - dy/10.
    max_p = MAX(data_2b(ii, index_h, index)) + dy/10.
   ENDIF ELSE BEGIN
    dy = MAX(data_2b(ii, index_h, *)) - MIN(data_2b(ii, index_h, *))
    min_p = MIN(data_2b(ii, index_h, *)) - dy/10.
    max_p = MAX(data_2b(ii, index_h, *)) + dy/10.
   ENDELSE

   IF min_p LT 0. THEN min_p = 0.
   max_p = CEIL(max_p)

   IF user_range_min(ii+2) NE -999. THEN min_p = user_range_min(ii+2)
   IF user_range_max(ii+2) NE -999. THEN max_p = user_range_max(ii+2)
  
   pos_cont = [0.1, 0.1, 0.88, 0.85]
   pos_flag = [0.1, 0.85, 0.88, 0.9]
   pos_bar = [0.93, 0.5, 0.97, 0.85]

   tit = ytit_2b(ii) + ' ('+ algo_2b(ii) +'), ' + mwr_meas

   PLOT_CONT, REFORM(data_2b(ii, *, *)), max_p, min_p, plot.time, flagx, z_final, h_min, h_max,$
              plot_start, plot_end, cont_dist(ii), units_2b(ii), date, pos_cont, pos_bar, pos_flag,$
              tit, deltat2b, home_path

  ENDFOR ; loop over level2b products
  
 ENDIF ; n_algos

ENDIF ; n_i_ang

SKIP_PLOT2:
DEVICE, /CLOSE

;***plot level2c (T&Tpot profiles, q profile&IWV, rh_profile&LWP) data
IF N_ELEMENTS(algo_tel) EQ 0 OR l1c.time(0) LT 0 THEN BEGIN
 print, 'no TEL retrieval specified - nothing to plot'
 GOTO, SKIP_PLOT3
ENDIF

plot_start = l1c.time(0)
plot_end = l1c.time(N_ELEMENTS(l1c.time)-1)
data_2c = REPLICATE(-999., 4, N_ELEMENTS(z_final), N_ELEMENTS(par2c.time))
units_2c = ['K', 'K', 'kgm!e-2!N', 'gm!e-3!N', 'gm!e-2!N', ' %']
ytits_2c = ['Temperature', 'Eq. pot. temperature', 'IWV / kgm!e-2!N', 'Abs. humidity', 'LWP gm!e-3!N', 'Relative humidity' ] 
plots_2c = ['t', 'tpot', 'qiwv', 'rhlwp']

data_2c(0, *, *) = par2c.T(*, *)
data_2c(1, *, *) = par2c.Tepot(*, *) 
data_2c(2, *, *) = par2c.q(*, *)*1000. 
data_2c(3, *, *) = par2c.rh(*, *) 

index = WHERE(flag_2c EQ 0, count)
index_h = WHERE(z_final LE h_max)

FOR i = 0, 3 DO BEGIN 
 
 IF verbose then print, 'creating '+date+'_'+mwr_meas+'_l2c_'+plots_2c(i)+'.ps'
 SET_PLOT, 'ps'
 LOADCT,39
 !X.minor = 6

 DEVICE, /color, /landscape, filename = plot_path2+date+'_'+mwr_meas+'_l2c_'+plots_2c(i)+'.ps'

 IF plot_range EQ 'flag' AND count GT 1 THEN BEGIN
  dy = MAX(data_2c(i, index_h, index)) - MIN(data_2c(i, *, index))
  min_p = MIN(data_2c(i, index_h, index)) - dy/10.
  max_p = MAX(data_2c(i, index_h, index)) + dy/10.
 ENDIF ELSE BEGIN
  dy = MAX(data_2c(i, index_h, *)) - MIN(data_2c(i, *, *))
  min_p = MIN(data_2c(i, index_h, *)) - dy/10.
  max_p = MAX(data_2c(i, index_h, *)) + dy/10.
 ENDELSE

 IF min_p LT 0. THEN min_p = 0.
 max_p = CEIL(max_p)
  
 pos_cont = [0.1, 0.1, 0.88, 0.85]
 pos_flag = [0.1, 0.85, 0.88, 0.9]
 pos_bar = [0.93, 0.5, 0.97, 0.85]

 IF i LT 2 THEN BEGIN 

  IF i EQ 0 THEN tit = 'Temperature (level2c), '+mwr_meas 
  IF i EQ 1 THEN tit = 'Eq. pot. temperature (level2c), '+mwr_meas

  IF user_range_min(3) NE -999. THEN min_p = user_range_min(3)
  IF user_range_max(3) NE -999. THEN max_p = user_range_max(3)

  PLOT_CONT, REFORM(data_2c(i, *, *)), max_p, min_p, par2c.time, flag_2c, z_final, h_min, h_max,$
             plot_start, plot_end, cont_dist(0), 'K', date, pos_cont, pos_bar, pos_flag,$
             tit, deltat2c, home_path

 ENDIF ELSE IF i EQ 2 THEN BEGIN

  pos_cont = [0.1, 0.1, 0.88, 0.7]
  pos_flag = [0.1, 0.85, 0.88, 0.9]
  pos_bar = [0.93, 0.5, 0.97, 0.85]
  pos_ts = [0.1, 0.7, 0.88, 0.85]
  tit = 'Abs. humidity and IWV (level2c), '+mwr_meas

  index = WHERE(flag_2a EQ 0, count)
  IF plot_range EQ 'flag' AND count GT 1 THEN BEGIN
   dy = MAX(data_2a(0, index)) - MIN(data_2a(0, index))
   ymin = MIN(data_2a(0, index)) - dy/10.
   ymax = MAX(data_2a(0, index)) + dy/10.
  ENDIF ELSE BEGIN
   dy = MAX(data_2a(0, *)) - MIN(data_2a(0, *))
   ymin = MIN(data_2a(0, *)) - dy/10.
   ymax = MAX(data_2a(0, *)) + dy/10.  
  ENDELSE

  IF user_range_min(0) NE -999. THEN ymin = user_range_min(0)
  IF user_range_max(0) NE -999. THEN ymax = user_range_max(0)

  IF user_range_min(4) NE -999. THEN min_p = user_range_min(4)
  IF user_range_max(4) NE -999. THEN max_p = user_range_max(4)

  PLOT_CONT_TS, REFORM(data_2c(i, *, *)), data_2a(0, *), max_p, min_p, ymax, ymin, par2c.time, plot.time, flag_2c, z_final, h_min, h_max,$
                plot_start, plot_end, cont_dist(1), 'gm!e-3!N', 'IWV kgm!e-2!N', date, pos_cont, pos_bar, pos_flag,$
                pos_ts, tit, deltat2c, home_path

 ENDIF ELSE IF i EQ 3 THEN BEGIN

  pos_cont = [0.1, 0.1, 0.88, 0.7] 
  pos_flag = [0.1, 0.85, 0.88, 0.9]
  pos_bar = [0.93, 0.5, 0.97, 0.85]
  pos_ts = [0.1, 0.7, 0.88, 0.85]
  tit = 'Rel. humidity and LWP (level2c), '+mwr_meas

  index = WHERE(flag_2a EQ 0, count)
  IF plot_range EQ 'flag' AND count GT 1 THEN BEGIN
   dy = MAX(data_2a(1, index)) - MIN(data_2a(1, index))
   ymin = MIN(data_2a(1, index)) - dy/10.
   ymax = MAX(data_2a(1, index)) + dy/10.
  ENDIF ELSE BEGIN
   dy = MAX(data_2a(1, *)) - MIN(data_2a(1, *))
   ymin = MIN(data_2a(1, *)) - dy/10.
   ymax = MAX(data_2a(1, *)) + dy/10.  
  ENDELSE

  IF user_range_min(1) NE -999. THEN ymin = user_range_min(1)
  IF user_range_max(1) NE -999. THEN ymax = user_range_max(1)

  IF user_range_min(5) NE -999. THEN min_p = user_range_min(5) 
  IF user_range_max(5) NE -999. THEN max_p = user_range_max(5)

  PLOT_CONT_TS, REFORM(data_2c(i, *, *)), data_2a(1, *), max_p, min_p, ymax, ymin, par2c.time, plot.time, flag_2c, z_final, h_min, h_max,$
                plot_start, plot_end, cont_dist(2), '%', 'LWP gm!e-2!N', date, pos_cont, pos_bar, pos_flag,$
                pos_ts, tit, deltat2c, home_path

 ENDIF
  
ENDFOR ; loop over level2c plots
  
SKIP_PLOT3:
DEVICE, /CLOSE

;***plot cloud radar and level2a LWP data
IF N_ELEMENTS(radar) EQ 0 THEN radar = 0

IF radar EQ 1 THEN BEGIN
;**plot all radar measurements (also off-zenith measurements, azimuth and elevation angle)

 filename1_plot_radar = plot_path2 + date + '_' + mwr_meas + '_l2_radar_lwp.ps'
 filename2_plot_radar = plot_path2 + date + '_' + mwr_meas + '_l2_radar.ps'

 print, 'PL_MK: creating radar quicklooks: '
 print, filename1_plot_radar
 print, filename2_plot_radar

 lwp_min_rad = MIN(data_2a(1,*))
 lwp_max_rad = MAX(data_2a(1,*))
 
 IF user_range_min(1) NE -999. THEN lwp_min_rad = user_range_min(1)
 IF user_range_max(1) NE -999. THEN lwp_max_rad = user_range_max(1)

 PLOT_RADAR_ALL, data_path_radar, filter_file_radar, filename1_plot_radar,$
                 filename2_plot_radar, date, mwr_meas, plot.time, REFORM(data_2a(1,*)), lwp_min_rad, lwp_max_rad

ENDIF

;**create CILI (Ceilometer, IRT, LWP, IWV)
IF N_ELEMENTS(ceilo_type) GT 0 AND N_ELEMENTS(algo_iwv) GT 0 AND N_ELEMENTS(algo_lwp) GT 0 THEN BEGIN

 SET_PLOT, 'ps'
 DEVICE, /color, BITS_PER_PIXEL = 8, file=plot_path2+date+'_'+mwr_meas+'_l2_cili.ps'
 tit = date + ', ' + mwr_meas
 !P.multi = [0, 1, 4]
 !P.charsize = 1.0

;*ceilometer backscatter plot
 fac = 0.9 ; (factor for scaling color bars)
 mu='!9b'
 bar_tit1 = 'LOG(!9b!3), CT25K' 
 bar_tit2 = '10!e-4!N (srad m)!e-1!N'

 LOADCT, 39

 y_min = [0.1, 0.3, 0.5, 0.7]
 y_max = [0.3, 0.5, 0.7, 0.9]

 x_min = 0.07
 x_max = 0.85

;bscat(height, time)

;*max/min bscat
 bscat_max = -999.
 n_time_ceilo = N_ELEMENTS(bscat(0, *))
 n_z_ceilo = N_ELEMENTS(z_ceilo)
 bscat_min = 1e-3

 FOR i = 0, n_z_ceilo-1 DO BEGIN
  ii = WHERE(bscat(i, *) LE 0.)
  IF ii(0) NE -1 THEN bscat(i, ii) = bscat_min
  bscat(i, *) = ALOG10(bscat(i, *))  

  bscat_max_x = MAX(bscat(i, *))
  IF bscat_max_x GT bscat_max THEN bscat_max = bscat_max_x
 ENDFOR

;*create homogeneous time array of 15s
 n_res = 15.
 n_time_c_h = 3600.*24./n_res
 time_c_h = FINDGEN(n_time_c_h)*24./n_time_c_h
 DATE_TIME_TO_JULDAT, '20'+date, time_c_h, time_c_h_jd

;*create homogeneous array of bscat on 15s
 bscat_h = REPLICATE(ALOG(bscat_min), n_z_ceilo, n_time_c_h)

;*put bscat on bscat_h array

 FOR i = 0, n_time_c_h-1 DO BEGIN
 
  ii = WHERE(ABS(time_c_h_jd(i)-time_ceilo) LE n_res/(3600.*24.)) 
  IF ii(0) NE -1 THEN BEGIN
   CLOSEST, time_c_h(i), time_ceilo(ii), dmin, ind
   bscat_h(*, i) = bscat(*, ii(ind))
  ENDIF
 ENDFOR

 arr = DBZ_2_BYTE(bscat_h, bscat_max + bscat_max*0.0001, ALOG10(bscat_min), fac)
 arr = TRANSPOSE(arr)

 TV, arr, x_min, y_min(3), xsize = 0.78, ysize = 0.2, /normal

 PLOT, bscat_h(0, *), xcharsize = 0.01, xrange = [1., n_time_ceilo],$
       YRANGE = [z_ceilo(0), z_ceilo(n_z_ceilo-1)], color = 0, title = tit, xticks = 12, xminor = 2,$
       charsize = 1.2, xstyle = 1, ystyle = 1, yminor = -1, ycharsize = 1.0,$
       Position = [x_min, y_min(3), x_max, y_max(3)], /Nodata, /Normal, font = 3, ytitle = 'Height [km]',$
       yticklen = -0.007, xticklen = 0.05, /NOERASE

 x1 = DBZ_2_BYTE(bscat_max, bscat_max + bscat_max*0.0001, bscat_min, fac)
 x2 = DBZ_2_BYTE(bscat_min, bscat_max + bscat_max*0.0001, bscat_min, fac)
 ncolor = ABS(FIX(x2-x1))

 BarColor = BINDGEN(ncolor) + FIX(MIN([x1, x2]))
 BarLabel = [bscat_min, bscat_max]
 !P.Color = 0.

 COLLABEL2, BarColor, BarLabel, CHARSIZE = 1.,$
            DIRECTION = 1, FONT = 0,$
            LabFormat = '(F5.2)',LABOFFSET = 35,$
            POSITION = [x_max+0.05, y_min(3), x_max+0.08, y_max(3)]

 XYOUTS, x_max+0.065, y_max(3)+0.02, bar_tit1, color = 0, font = 3, charsize = 0.8, alignment = 0.5, /Normal
 XYOUTS, x_max+0.065, y_min(3)-0.03, bar_tit2, color = 0, font = 3, charsize = 0.8, alignment = 0.5, /Normal

;*IRT plot (if data available)

 LOADCOL, home_path+'/mwr_pro/source/col1'

 ii = WHERE(tb_irt(0, *) NE -999.)

 IF ii(0) NE -1 THEN BEGIN

  tb_irt_plot = REPLICATE(-999., n_wavel_irt, N_ELEMENTS(plot.time))
  FOR i = 0, n_wavel_irt-1 DO BEGIN
   tb_irt_plot(i, *) = tb_irt(i, i_ang)
  ENDFOR 
   
  PLOT, plot.time, tb_irt_plot(0, *), xcharsize = 0.01, xrange = [0., 24.],$
        YRANGE = [irt_range_min, irt_range_max], color = 0, xticks = 12, xminor=2,$
        charsize = 1.2, xstyle = 1, ystyle = 1, yminor = -1, ycharsize = 1.0,$
        Position = [x_min, y_min(2), x_max, y_max(2)], /Nodata, /Normal, ytitle = 'TB [K]',$
        /NOERASE

  char = ['IRT@'+STRING(irt.wavel(0), format='(f4.1)')+'!9m!3m', 'IRT@'+STRING(irt.wavel(1), format='(f4.1)')+'!9m!3m']

  FOR i = 0, n_wavel_irt-1 DO BEGIN
   OPLOT, plot.time, tb_irt_plot(i, *), psym = 3, color= i+1, thick = 3
   XYOUTS, 1.,  irt_range_max - (i+1)*(irt_range_max-irt_range_min)/10., char(i), color = i+1, charsize = 0.5
  ENDFOR 

 ENDIF

;*LWP plot
 ii = WHERE(data_2a(1, *) NE -999.)
 IF ii(0) NE -1 THEN BEGIN

  index = WHERE(flag_2a(i_ang) EQ 0, count)
  IF plot_range EQ 'flag' AND count GT 1 THEN BEGIN
   dy = MAX(data_2a(1, index)) - MIN(data_2a(1, index))
   ymin = MIN(data_2a(1, index)) - dy/10.
   ymax = MAX(data_2a(1, index)) + dy/10.
  ENDIF ELSE BEGIN
   dy = MAX(data_2a(1, *)) - MIN(data_2a(1, *))
   ymin = MIN(data_2a(1, *)) - dy/10.
   ymax = MAX(data_2a(1, *)) + dy/10.
  ENDELSE
   
  PLOT, plot.time, data_2a(1, *), xcharsize = 0.01, xrange = [0., 24.],$
        YRANGE = [ymin, ymax], color = 0, xticks = 12, xminor = 2,$
        charsize = 1.2, xstyle = 1, ystyle = 1, ycharsize = 1.0, yminor = 0, $
        Position = [x_min, y_min(1), x_max, y_max(1)], /Nodata, /Normal, font = 3, ytitle = 'LWP [gm!e-2!N]',$
        /NOERASE

  OPLOT, plot.time, data_2a(1, *), psym = 3, color= 5, thick = 3
  OPLOT, [0., 24.], [0., 0.], linestyle = 2, color = 0

 ENDIF

;*IWV plot
 ii = WHERE(data_2a(0, *) NE -999.)
 IF ii(0) NE -1 THEN BEGIN

  index = WHERE(flag_2a(i_ang) EQ 0, count)
  IF plot_range EQ 'flag' AND count GT 1 THEN BEGIN
   dy = MAX(data_2a(0, index)) - MIN(data_2a(0, index))
   ymin = MIN(data_2a(0, index)) - dy/10.
   ymax = MAX(data_2a(0, index)) + dy/10.
  ENDIF ELSE BEGIN
   dy = MAX(data_2a(0, *)) - MIN(data_2a(0, *))
   ymin = MIN(data_2a(0, *)) - dy/10.
   ymax = MAX(data_2a(0, *)) + dy/10.
  ENDELSE
   
  PLOT, plot.time, data_2a(0, *), xrange = [0., 24.],$
        YRANGE = [ymin, ymax], color = 0,$
        charsize = 1.2, xstyle = 1, ystyle = 1, ycharsize = 1.0, xticks = 12, xminor = 2, $
        Position = [x_min, y_min(0), x_max, y_max(0)], /Nodata, /Normal, font = 3, ytitle = 'IWV [kgm!e-2!N]',$
        xticklen = 0.05, xtit = 'Time [UTC]', /NOERASE

  OPLOT, plot.time, data_2a(0, *), psym = 3, color= 6, thick = 3

 ENDIF

 DEVICE, /close

ENDIF ; ceilo_type specified

SKIP_PLOT:
END

