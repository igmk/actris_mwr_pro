;+
;***********
PRO RET_2B,$
;***********
;INPUT
ret_files,$            ;set of retrieval files
ret_file=ret_file,$    ;'ret'--> use RPG-formatted files (default), 'sav' --> use IDL binary files (not coded yet)
tb_org,$               ;(n_freq x n_meas)
tb,$                   ;(n_freq x n_meas) offset corrected
elev,$                 ;(n_meas)
freq,$                 ;frequency of tb measurements
time,$                 ;time of measurements
date,$                 ;yymmdd of measurements
;OUTPUT
z,$                    ;height grid in m
predictand,$           ;T or q (nz x n_meas)
dpredictand,$          ;predictand_org - predictand_oc
predictand_err_a,$     ;standard error of predictand
angle_a,$              ;retrieval angles
;KEYWORDS
ret_par=ret_par,$      ;retrieval parameters contained in nc files
tb_t=tb_t,$            ;T / K on TB time grid
tb_q=tb_q,$            ;q / kgm-3 on TB time grid
tb_p=tb_p,$            ;p / Pa on TB time grid
verbose=verbose
; Abstract:
;* read set of .RET or .nc files specified in ret_files
;  and calculate par (=q or T-profile) at defined 
;  elevation angles (elev). 
; Author:
; U. Loehnert
; Date:
; 2011-02-11
; Dependencies:
; -
; Changes:
; 20121001 (UL)
; * parameter par_all was redundant and lead to errors in case of multiple 
;   retrieval files; par_all now removed
; 2014-10-27 (UL):
; * added the possibility to read netcdf retrieval files, containing also relevant retrieval parameters

;-

n_freq = 100
n_ret = N_ELEMENTS(ret_files)
n_meas = N_ELEMENTS(time)
angle_a = FLTARR(n_ret)

IF N_ELEMENTS(ret_files) EQ 0 OR ret_files(0) EQ '' THEN BEGIN
 print, 'RET_2B: No retrieval files found!'
 stop
ENDIF

IF ret_file EQ 'RET' OR N_ELEMENTS(ret_file) EQ 0 THEN BEGIN

 GET_COEFF_LEVEL2B, ret_files(0), 1, angle, f, z, offset, c_l, c_q, verbose=0

ENDIF ELSE IF ret_file EQ 'nc' THEN BEGIN

 GET_COEFF_LEVEL2B_NC, ret_files(0), ret_par, angle, f, z, offset, c_l, c_q, c_s, predictand_err, verbose=1

ENDIF

;****fill missing surface values with nearest neighbor values
tb_p = tb_p*100. 
T_help = -999.
q_help = -999.
p_help = -999.
time_help = -999.

IF ret_par.surf EQ 'in_surface' THEN BEGIN
 ii = WHERE(tb_t GT 250. AND tb_t LT 330. AND tb_q GE 0. AND tb_q LT 0.02 AND tb_p GT 0. AND tb_p LT 120000.)
 T_help = tb_t(ii)
 q_help = tb_q(ii) 
 p_help = tb_p(ii)
 time_help = time(ii)

 FOR i = 0l, N_ELEMENTS(time)-1l DO BEGIN
  IF tb_t(i) GT 250. AND tb_t(i) LT 330. AND tb_q(i) GE 0. AND tb_q(i) LT 0.02 AND $
     tb_p(i) GT 0. AND tb_p(i) LT 120000. THEN BEGIN
;  OK!!
  ENDIF ELSE BEGIN
   CLOSEST, time(i), time_help, dmin, index
   tb_t(i) = T_help(index)
   tb_q(i) = q_help(index)
   tb_p(i) = p_help(index)
  ENDELSE
 ENDFOR

 surf = TRANSPOSE([[tb_t], [tb_q], [tb_p]])
ENDIF

nz = N_ELEMENTS(z)
predictand = REPLICATE(-999., nz, n_meas)
dpredictand = REPLICATE(0., nz, n_meas)

;***read coefficents
predictand_err_a = REPLICATE(-999., n_ret, nz)

FOR i_ret = 0, n_ret-1 DO BEGIN

 IF ret_file EQ 'RET' OR N_ELEMENTS(ret_file) EQ 0 THEN BEGIN

  GET_COEFF_LEVEL2B, ret_files(i_ret), 1, angle, f, z, offset, c_l, c_q, verbose=0
  n_f = N_ELEMENTS(f)
  nz = N_ELEMENTS(z)
  angle_a(i_ret) = angle

 ENDIF ELSE IF ret_file EQ 'nc' THEN BEGIN

  GET_COEFF_LEVEL2B_NC, ret_files(i_ret), ret_par, angle, f, z, offset, c_l, c_q, c_s, predictand_err, verbose=1
  n_f = N_ELEMENTS(f)
  nz = N_ELEMENTS(z)
  predictand_err_a(i_ret, *) = predictand_err
  angle_a(i_ret) = angle

 ENDIF

;*check if chosen retrieval file is compatible with the frequency channels and elevation angles
; extracted from the measurements
 i_ang = WHERE(elev GT angle-0.6 AND elev LT angle+0.6) ; 0.6 deg is the typical pointing accuracy of HATPRO
 IF i_ang(0) EQ -1 THEN BEGIN
 
  print, 'RET_2B: No retrieval at elevation angle '+STRING(angle, format='(f6.2)')+ ' carried out'
  print, 'because no measurements at this angle were found!'
 
 ENDIF ELSE BEGIN
  
  n_f = N_ELEMENTS(f) 
  go = 0
  FOR i = 0, n_f-1 DO BEGIN
   im = WHERE(freq EQ f(i))
   IF im(0) EQ -1 THEN BEGIN
    print, 'ABORTING in RET_2B: Retrieval file is not compatible with measurements'
    GOTO, ABORT
   ENDIF
   IF go EQ 0 THEN BEGIN
    ima = im
    go = 1
   ENDIF ELSE BEGIN
    ima = [ima, im] 
   ENDELSE
  ENDFOR

;***brightness temperature for retrieval
  tb_algo = tb(ima, *)
  tb_algo = tb_algo(*, i_ang)
  tb_algo_org = tb_org(ima, *)
  tb_algo_org = tb_algo_org(*, i_ang)
  time_algo = time(i_ang) 
  n_ang = N_ELEMENTS(i_ang)

;***calculate predictand
  FOR i = 0l, LONG(n_ang-1l) DO BEGIN

   FOR j = 0, nz-1 DO BEGIN
    IF ret_par.reg EQ 'linear' THEN BEGIN  
     predictand(j, i_ang(i)) = offset(j) + c_l(j, *)#tb_algo(0:n_f-1, i)
     dpredictand(j, i_ang(i)) = offset(j) + c_l(j, *)#tb_algo_org(0:n_f-1, i) - predictand(j, i_ang(i)) 
     IF ret_par.surf EQ 'in_surface' THEN BEGIN 
      predictand(j, i_ang(i)) = offset(j) + c_l(j, *)#tb_algo(0:n_f-1, i) + c_s(j, *)#surf(*, i) 
      dpredictand(j, i_ang(i)) = offset(j) + c_l(j, *)#tb_algo_org(0:n_f-1, i) + c_s(j, *)#surf(*, i) - predictand(j, i_ang(i))   
     ENDIF
        
    ENDIF ELSE IF ret_par.reg EQ 'quadratic' THEN BEGIN
     predictand(j, i_ang(i)) = offset(j) + c_l(j, *)#tb_algo(0:n_f-1, i) + c_q(j, *)#tb_algo(0:n_f-1, i)^2
     dpredictand(j, i_ang(i)) = offset(j) + c_l(j, *)#tb_algo_org(0:n_f-1, i) + c_q(j, *)#tb_algo_org(0:n_f-1, i)^2 - predictand(j, i_ang(i))
     IF ret_par.surf EQ 'in_surface' THEN BEGIN
      predictand(j, i_ang(i)) = offset(j) + c_l(j, *)#tb_algo(0:n_f-1, i) + c_q(j, *)#tb_algo(0:n_f-1, i)^2 + c_s(j, *)#surf(*, i) 
      dpredictand(j, i_ang(i)) = offset(j) + c_l(j, *)#tb_algo_org(0:n_f-1, i) + c_q(j, *)#tb_algo_org(0:n_f-1, i)^2 + c_s(j, *)#surf(*, i) - predictand(j, i_ang(i))
     ENDIF
    ENDIF 
   ENDFOR
  ENDFOR

 ENDELSE

ENDFOR ; ret_files loop

ABORT:
END
