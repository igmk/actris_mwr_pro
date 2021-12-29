;+
;***********
PRO RET_2C,$
;***********
;INPUT
ret_files,$            ;specified retrieval file (one only)
ret_file=ret_file,$    ;'ret'--> use RPG-formatted files (default), 'sav' --> use IDL binary files (not coded yet)
time,$                 ;time array in hours
tb_org,$               ;(frequency x angles x time)
tb,$                   ;(frequency x angles x time) offset corrected
elev,$                 ;(angles x time) --> elevation angles of measured TBs
freq,$                 ;frequency of TB measurements
met,$                  ;environmental variables (structure): {time:time_m, p:pres_m, t:temp_m, q:humi_m, r:rain_m}
time2b,$               ;level2b retrieval times
q2b_org,$              ;level2b absolute humidity data
q2b,$                  ;level2b absolute humidity data offset corrected
hze_path,$             ;path of level2 data
fn_hze,$               ;file naming convention of humidity profiles
pres2b,$               ;surface pressure on level2b time grid
mwr_data_path,$        ;
date,$                 ;yymmdd       
mwr_meas,$             ;MWR measurement (user-defined three letter code)
;OUTPUT
predictand,$           ;retrieved parameters (structure): {time:time, z:z, T:T, q:q, rh:rh, Tpot:Tpot, Tepot:Tepot}
dpredictand,$          ;retrieved parameters (structure) : predictand_org - predictand_oc
predictand_err,$       ;standard error of retrieved temperature
;KEYWORDS:
ret_par=ret_par,$      ;retrieval parameters contained in nc files
verbose=verbose
; Abstract:
;* read specified *.RET file (ret_files)
;  and calculate par profiles: 
;  - temperature
;  - relative humidity
;  - potential temperature
;  - equivalent potential temperature
; * routine needs existing 2b data to work correctly!
; --> either q2b needs to be passed on from prior RET_2B call or
; --> this routine searches for existing level2b netcdf data
; Author:
; U. Loehnert
; Date:
; 2011-02-11
; Dependencies:
; -
; Changes:
; *20120518 (UL): BUG FIXED
;  The mutli-angle, multi-freqeuncy TB-array-rearangement for mulitplication with the regression coefficients was wrong.
;  This lead to faulty wave-like structures in the derived temperature profile. TB-array-rearangement now OK.  
; *20120725 (UL)
;  Calculation of potential temperature is now skipped if no surface pressure information is available. 
; *2014-10-27 (UL):
;  added the possibility to read netcdf retrieval files, containing also relevant retrieval parameters
; *2015-01-08 (UL):
;  Error in Tv and eq. pot. temp. calculation detected - division of q by 1000. was wrong; corrected.
;- 

IF N_ELEMENTS(time) EQ 0 OR time(0) EQ -999. THEN BEGIN
 print, 'RET_2C: No level1c data found - aborting level2c retrieval'
ENDIF

g = 9.81
Rl = 287.
Rw = 462.
cp = 1004.

hsb = 300./3600. ; averaging window for level2b data

IF N_ELEMENTS(time2b) EQ 0 THEN BEGIN
;***search for existing level2b netcdf data

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
 file_hze = hze_path + fn_hze.kkk + '_' + fn_hze.sss + '_' + fn_hze.inst + fn_hze.instnn + '_'+ fn_hze.lll2 + '_' + fn_hze.var + '_' + fn_hze.vnn + '_' + date_string + '.nc'

 filename = FILE_SEARCH(file_hze)
 IF filename EQ '' THEN BEGIN
  print, 'RET_2C: No 2b data found - aborting level2c retrieval'
  GOTO, ABORT
 ENDIF
 READ_LEVEL2B_NC, filename, algo, comment, time2b, z2b, T2b, q2b, flag2b, temp2b, pres2b, relh2b
 n_time2b = N_ELEMENTS(time2b)
ENDIF

IF verbose THEN print, 'RET_2C: Number of total time steps level2b: ', N_ELEMENTS(time2b)

dummy = -999.
n2c = N_ELEMENTS(time)
nz2b = N_ELEMENTS(q2b(*, 0))
T = REPLICATE(dummy, nz2b, n2c)
rh = REPLICATE(dummy, nz2b, n2c)
Theta = REPLICATE(dummy, nz2b, n2c)
Theta_e = REPLICATE(dummy, nz2b, n2c)
rh = REPLICATE(dummy, nz2b, n2c)
q = REPLICATE(dummy, nz2b, n2c)
es = REPLICATE(dummy, nz2b, n2c)
e = REPLICATE(dummy, nz2b, n2c)
Tv = REPLICATE(dummy, nz2b, n2c)
mr = REPLICATE(dummy, nz2b, n2c)

T_org = REPLICATE(dummy, nz2b, n2c)
rh_org = REPLICATE(dummy, nz2b, n2c)
Theta_org = REPLICATE(dummy, nz2b, n2c)
Theta_e_org = REPLICATE(dummy, nz2b, n2c)
rh_org = REPLICATE(dummy, nz2b, n2c)
q_org = REPLICATE(dummy, nz2b, n2c)
es_org = REPLICATE(dummy, nz2b, n2c)
e_org = REPLICATE(dummy, nz2b, n2c)
Tv_org = REPLICATE(dummy, nz2b, n2c)
mr_org = REPLICATE(dummy, nz2b, n2c)

;***read coefficents

IF ret_file EQ 'RET' OR N_ELEMENTS(ret_file) EQ 0 THEN BEGIN

 GET_COEFF_LEVEL2C, ret_files, angle, f, z, offset, coeff_tel, coeff_sen, f_bl, verbose=0

ENDIF ELSE IF ret_file EQ 'nc' THEN BEGIN

 GET_COEFF_LEVEL2C_NC, ret_files, ret_par, angle, f, z, offset, coeff_tel, f_bl, predictand_err, verbose=0

ENDIF

nf = N_ELEMENTS(f)
nz = N_ELEMENTS(z)
na = N_ELEMENTS(angle)
nbl = N_ELEMENTS(f_bl)
nzen = nf - nbl

;*check if chosen retrieval file is compatible with the measured frequency channels and angles

go = 0
ibl = -1
izen = -1
FOR i = 0, nf-1 DO BEGIN 

 im = WHERE(freq EQ f(i))
 im = im(0)
 IF im EQ -1 THEN BEGIN
  print, 'RET_2C ABORT: Frequencies in retrieval file are not compatible with measurements'
  print, 'retrieval file frequencies: ', f
  print, 'measurement frequencies: ', freq
  GOTO, ABORT
 ENDIF
 IF go EQ 0 THEN BEGIN
  imf = im
  go = 1
 ENDIF ELSE BEGIN
  imf = [imf, im]
 ENDELSE

ENDFOR

FOR i = 0, nbl-1 DO BEGIN
 ibb = WHERE(f_bl(i) EQ f)
 IF ibb(0) NE -1 THEN BEGIN
  IF ibl(0) EQ -1 THEN BEGIN
   ibl = ibb
  ENDIF ELSE BEGIN
   ibl = [ibl, ibb]
  ENDELSE
 ENDIF 

ENDFOR

go = 0
FOR j = 0, na-1 DO BEGIN
 im = WHERE(elev GT angle(j)-0.6 AND elev LT angle(j)+0.6) 
 im = im(0)
 IF im EQ -1 THEN BEGIN

  print, 'RET_2C ABORT: Elevation angles in retrieval file are not compatible with measurements'
  print, 'retrieval file angles: ', angle
  print, 'measurement angle: ', elev
 
  GOTO, ABORT

 ENDIF ELSE BEGIN

  IF go EQ 0 THEN BEGIN
   ima = im
   go = 1
  ENDIF ELSE BEGIN
   ima = [ima, im]
  ENDELSE
   
 ENDELSE
ENDFOR

;***brightness temperature for retrieval
tb_algo = tb(imf, ima, *)
tb_algo_org = tb_org(imf, ima, *)

FOR i = 0, n2c-1 DO BEGIN

 b_match = WHERE(time2b GE time(i)-hsb AND time2b LT time(i)+hsb and q2b(0, *) GE 0)

 IF b_match(0) NE -1 THEN BEGIN
  ivalpres = WHERE(pres2b(b_match) GT 0.)
  IF ivalpres(0) NE -1 THEN pres_surf = MEAN(pres2b(b_match(ivalpres)))
  IF ivalpres(0) EQ -1 THEN print, 'RET_2C: no surface pressure data available!!'
 ENDIF

;**rearange tbs

 IF nzen GE 1 THEN BEGIN 
  tb_algo1 = REFORM(tb_algo(0:nzen-1, 0, i))
  FOR j = 0, nbl-1 DO tb_algo1 = [tb_algo1, REFORM(tb_algo(ibl(j), *, i))] 
 ENDIF ELSE IF nzen EQ 0 THEN BEGIN
  tb_algo1 = REFORM(tb_algo(ibl(0), *, i))
  FOR j = 1, nbl-1 DO tb_algo1 = [tb_algo1, REFORM(tb_algo(ibl(j), *, i))]
 ENDIF 

 IF nzen GE 1 THEN BEGIN 
  tb_algo1_org = REFORM(tb_algo_org(0:nzen-1, 0, i))
  FOR j = 0, nbl-1 DO tb_algo1_org = [tb_algo1_org, REFORM(tb_algo_org(ibl(j), *, i))] 
 ENDIF ELSE IF nzen EQ 0 THEN BEGIN
  tb_algo1_org = REFORM(tb_algo_org(ibl(0), *, i))
  FOR j = 1, nbl-1 DO tb_algo1_org = [tb_algo1_org, REFORM(tb_algo_org(ibl(j), *, i))]
 ENDIF 

 FOR j = 0, nz-1 DO T(j, i) = offset(j) + coeff_tel(j, *)#tb_algo1  
 FOR j = 0, nz-1 DO T_org(j, i) = offset(j) + coeff_tel(j, *)#tb_algo1_org  

;**calculate other 2c variables using 2b data
 IF b_match(0) NE -1 THEN BEGIN

  FOR j = 0, nz-1 DO BEGIN
   q(j, i) = MEAN(q2b(j, b_match)) ; absolute humidity in kg/m^-3
   q_org(j, i) = MEAN(q2b_org(j, b_match)) ; absolute humidity in kg/m^-3

   ABSHUM_TO_RH, T(j, i), q(j, i), rhx
   rh(j, i) = rhx
   Tv(j, i) = T(j, i)*(1+0.608*q(j, i))  ; virtual temperature in K

   ABSHUM_TO_RH, T_org(j, i), q_org(j, i), rhx_org
   rh_org(j, i) = rhx_org
   Tv_org(j, i) = T_org(j, i)*(1+0.608*q_org(j, i))  ; virtual temperature in K

  ENDFOR

;*calculate pressure in each level using barometric height formula
  IF N_ELEMENTS(pres_surf) GT 0 THEN BEGIN
   p_baro = FLTARR(nz)
   p_baro(0) = pres_surf*100.
   p_baro_org = FLTARR(nz)
   p_baro_org(0) = pres_surf*100.

   FOR jj = 1, nz-1 DO p_baro(jj) = p_baro(jj-1)*exp(-g*(z(jj)-z(jj-1))/(Rl*(Tv(jj, i)+Tv(jj-1, i))/2))  ; pressure on each level in Pa

   FOR j = 0, nz-1 DO BEGIN
    Theta(j, i) = T(j, i)*(100000./p_baro(j))^(Rl/cp) ; potential temperature in K
    e = 462.*T(j, i)*q(j, i)
    mr(j, i) = 0.622*e/(p_baro(j)-e) ; water vapor mixing ratio in kg/kg
    lv = (2500.-2.42*(T(j, i)-273.15))*1000. ; latent heat of vaporization in J/kg
    Theta_e(j, i) = Theta(j, i)+(lv*mr(j, i)/cp)*(100000./p_baro(j))^(Rl/cp) ; equivalent potential temperature in K
   ENDFOR

   FOR jj = 1, nz-1 DO p_baro_org(jj) = p_baro_org(jj-1)*exp(-g*(z(jj)-z(jj-1))/(Rl*(Tv_org(jj, i)+Tv_org(jj-1, i))/2))  ; pressure on each level in Pa

   FOR j = 0, nz-1 DO BEGIN
    Theta_org(j, i) = T_org(j, i)*(100000./p_baro_org(j))^(Rl/cp) ; potential temperature in K
    e_org = 462.*T_org(j, i)*q_org(j, i)
    mr_org(j, i) = 0.622*e_org/(p_baro_org(j)-e_org) ; water vapor mixing ratio in kg/kg
    lv_org = (2500.-2.42*(T_org(j, i)-273.15))*1000. ; latent heat of vaporization in J/kg
    Theta_e_org(j, i) = Theta_org(j, i)+(lv_org*mr_org(j, i)/cp)*(100000./p_baro_org(j))^(Rl/cp) ; equivalent potential temperature in K
   ENDFOR

  ENDIF ELSE BEGIN
    print, 'RET_2C: no calculation of potential temperature because no surface pressure data available'
  ENDELSE
 ENDIF

ENDFOR ; loop over n2c

;***calculate potential and equivalent potential temperature

ABORT:
predictand = {time:time, z:z, T:T, q:q, rh:rh, Tpot:Theta, Tepot:Theta_e}
dpredictand = {time:time, z:z, T:T_org-T, q:q_org-q, rh:rh_org-rh, Tpot:Theta_org-Theta, Tepot:Theta_e_org-Theta_e}

END
