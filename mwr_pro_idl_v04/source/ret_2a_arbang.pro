;+
;******************
PRO RET_2A_ARBANG,$
;******************
;INPUT
ret_files,$            ;set of retrieval files
ret_file=ret_file,$    ;'ret'--> use RPG-formatted files (default), 'ncdf' --> use netcdf files
time,$                 ;(n_meas) in decimal hours
tb_org,$               ;(n_freq x n_meas)
tb,$                   ;(n_freq x n_meas) offset corrected
elev,$                 ;(n_meas)
freq,$                 ;frequency of tb measurements
date,$                 ;yymmdd of measurements
;OUTPUT
predictand,$           ;IWV, LWP or WDL (n_meas) offset corrected
dpredictand,$          ;IWV-IWV_oc, LWP-LWP_oc or WDL-WDL_oc (n_meas)
predictand_err_a,$     ;standard error of predictand due to MVR contained in different (angle-dependent) retrieval files
angle_a,$              ;retrieval angles
ret_par=ret_par,$      ;retrieval parameters contained in nc files
tb_t=tb_t,$            ;T / K on TB time grid             
tb_q=tb_q,$            ;q / kgm-3 on TB time grid
tb_p=tb_p,$            ;p / Pa on TB time grid 
;KEYWORDS
verbose=verbose

; Abstract:
;* read set of .RET or .NC files specified in ret_files
;  and calculate predictand (=IWV, LWP or WDL) at arbitrary
;  elevation angles (elev). Note that if elev is not
;  equal to angle (see below), the quality of the retrieved
;  product will be dependent on ret_files because
;  the applied retrieval coefficients are calculated through
;  linear interpolation.
; Author:
; U. Loehnert
; Date:
; 2008-03-30
; Dependencies:
; -
; Changes:
; 2011-02-11 (UL): 
; * freq added as further INPUT variable
; * f removed as OUTPUT variable
; * renamed from RET_IWVLWP_MA to RET_2A_ARBANG
; * consistency check for frequencies in retrieval file and from measurements added
; 2014-10-27 (UL):
; * added the possibility to read netcdf retrieval files, containing also relevant retrieval parameters
; 2014-12-01 (UL):
; * possible tb offset correction now alreay calculated in pl_mk.pro; calculation taken out of ret_2a_arbang
;-

IF N_ELEMENTS(ret_files) EQ 0 OR ret_files(0) EQ '' THEN BEGIN
 print, 'RET_2A_ARBANG: No retrieval files found!'
 stop
ENDIF

IF ret_file EQ 'RET' OR N_ELEMENTS(ret_file) EQ 0 THEN BEGIN

 GET_COEFF_LEVEL2A, 0, ret_files(0), lin_qua, angle, f, offset, c_l, c_q

ENDIF ELSE IF ret_file EQ 'nc' THEN BEGIN

 GET_COEFF_LEVEL2A_NC, 0, ret_files(0), ret_par, angle, f, offset, c_l, c_q, c_s, predictand_err

ENDIF

n_f = N_ELEMENTS(f)
n_ret = N_ELEMENTS(ret_files)
offset_a = FLTARR(n_ret)
c_l_a = FLTARR(n_f, n_ret)
c_q_a = FLTARR(n_f, n_ret)
c_s_a = FLTARR(3, n_ret)
angle_a = FLTARR(n_ret)
surf = -999.

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

;***read coefficents
predictand_err_a = REPLICATE(-999., n_ret)

FOR j = 0, n_ret-1 DO BEGIN

 IF ret_file EQ 'RET' OR N_ELEMENTS(ret_file) EQ 0 THEN BEGIN

  GET_COEFF_LEVEL2A, 0, ret_files(j), lin_qua, angle, f, offset, c_l, c_q

 ENDIF ELSE IF ret_file EQ 'nc' THEN BEGIN

  GET_COEFF_LEVEL2A_NC, 0, ret_files(j), ret_par, angle, f, offset, c_l, c_q, c_s, predictand_err
  predictand_err_a(j) = predictand_err

 ENDIF

 angle_a(j) = angle
 offset_a(j) = offset

 IF N_ELEMENTS(c_l) NE N_ELEMENTS(c_l_a(*, j)) THEN BEGIN
  print, 'RET_2A_ARBANG: Inconsistency in specified retrieval files'
  stop
 ENDIF
 c_l_a(*, j) = c_l
 c_q_a(*, j) = c_q
 c_s_a(*, j) = c_s

;*check if chosen retrieval file is compatible with the freqeuncy channels extracted from the measurements
 go = 0
 FOR i = 0, n_f-1 DO BEGIN
  im = WHERE(freq EQ f(i))
  IF im(0) EQ -1 THEN BEGIN
   print, 'ABORTING RET_2A_ARBANG: Retrieval file is not compatible with measurements'
   print, 'measurements @ ', freq
   print, 'retrieval file frequency: ', f(i)
   print, 'Frequency mismatch!' 
   GOTO, NO_RETRIEVAL
  ENDIF
  IF j GT 0 THEN BEGIN
   IF TOTAL(f) NE TOTAL(f_old) THEN BEGIN
    print, 'ABORTING RET_2A_ARBANG: Retrieval file are not compatible with each other'
    print, 'Frequency mismatch!'
    GOTO, NO_RETRIEVAL
   ENDIF
  ENDIF
  IF go EQ 0 THEN BEGIN
   ima = im
   go = 1
  ENDIF ELSE BEGIN
   ima = [ima, im] 
  ENDELSE
 ENDFOR
 f_old = f

ENDFOR

;***brightness temperature for retrieval
tb_algo = tb(ima, *)
tb_algo_org = tb_org(ima, *)

i_sort = SORT(angle_a)
angle_a = angle_a(i_sort)
offset_a = offset_a(i_sort)
c_l_a = c_l_a(*, i_sort)
c_q_a = c_q_a(*, i_sort)
c_s_a = c_s_a(*, i_sort)

angle_c = angle_a
angle_c(N_ELEMENTS(angle_c)-1) = angle_c(N_ELEMENTS(angle_c)-1) + 0.6

n_meas = N_ELEMENTS(elev)
predictand = REPLICATE(-999., n_meas)
dpredictand = REPLICATE(0., n_meas)

;***calculate predictand
go = 0
FOR i = 0l, LONG(n_meas-1) DO BEGIN

 IF n_ret GT 1 THEN BEGIN

;**linearly interpolate regression coefficients to elevation angle of measurement
  FOR j = 0, n_ret-2 DO BEGIN

   IF elev(i) GE angle_c(j)-0.6 AND elev(i) LE angle_c(j+1)+0.6 THEN BEGIN

    steig_o = (offset_a(j+1)-offset_a(j))/(angle_a(j+1)-angle_a(j))
    offset_x = offset_a(j) + steig_o*(elev(i)-angle_a(j))

    c_l_x = FLTARR(n_f)
    c_q_x = FLTARR(n_f)
    c_s_x = FLTARR(3)

    FOR k = 0, n_f-1 DO BEGIN
     steig_l = (c_l_a(k, j+1)-c_l_a(k, j))/(angle_a(j+1)-angle_a(j))
     c_l_x(k) = c_l_a(k, j) + steig_l*(elev(i)-angle_a(j))
     steig_q = (c_q_a(k, j+1)-c_q_a(k, j))/(angle_a(j+1)-angle_a(j))
     c_q_x(k) = c_q_a(k, j) + steig_q*(elev(i)-angle_a(j))
    ENDFOR
    IF ret_par.surf EQ 'in_surface' THEN BEGIN
     FOR k = 0, 2 DO BEGIN
      steig_s = (c_s_a(k, j+1)-c_s_a(k, j))/(angle_a(j+1)-angle_a(j))
      c_s_x(k) = c_s_a(k, j) + steig_s*(elev(i)-angle_a(j))     
     ENDFOR
    ENDIF
    IF ret_par.reg EQ 'linear' THEN BEGIN
     predictand(i) = offset_x + TOTAL(c_l_x*tb_algo(0:n_f-1, i))
     dpredictand(i) = -1.*(predictand(i) - (offset_x + TOTAL(c_l_x*tb_algo_org(0:n_f-1, i))))
     IF ret_par.surf EQ 'in_surface' THEN BEGIN 
      predictand(i) = offset_x + TOTAL(c_l_x*tb_algo(0:n_f-1, i)) + TOTAL(c_s_x*surf(*, i))                
      dpredictand(i) = -1.*(predictand(i) + (offset_x + TOTAL(c_l_x*tb_algo_org(0:n_f-1, i)) + TOTAL(c_s_x*surf(*, i))))
     ENDIF
    ENDIF
    
    IF ret_par.reg EQ 'quadratic' THEN BEGIN
     predictand(i) = offset_x + TOTAL(c_l_x*tb_algo(0:n_f-1, i)) + TOTAL(c_q_x*(tb_algo(0:n_f-1, i)^2))
     dpredictand(i) = -1.*(predictand(i) - (offset_x + TOTAL(c_l_x*tb_algo_org(0:n_f-1, i)) + TOTAL(c_q_x*(tb_algo_org(0:n_f-1, i)^2))))
     IF ret_par.surf EQ 'in_surface' THEN BEGIN
      predictand(i) = offset_x + TOTAL(c_l_x*tb_algo(0:n_f-1, i)) + TOTAL(c_q_x*(tb_algo(0:n_f-1, i)^2)) + TOTAL(c_s_x*surf(*, i))                
      dpredictand(i) = -1.*(predictand(i) - (offset_x + TOTAL(c_l_x*tb_algo_org(0:n_f-1, i)) + TOTAL(c_q_x*(tb_algo_org(0:n_f-1, i)^2)) + TOTAL(c_s_x*surf(*, i))))
     ENDIF
    ENDIF

   ENDIF

  ENDFOR
 ENDIF ELSE BEGIN

  IF n_ret EQ 1 AND elev(i) GT angle_a(0)-0.6 AND elev(i) LT angle_a(0)+0.6 THEN BEGIN
   IF ret_par.reg EQ 'linear' THEN BEGIN
    predictand(i) = offset_a(0) + TOTAL(c_l_a(*, 0)*tb_algo(0:n_f-1, i))
    dpredictand(i) = -1.*(predictand(i) - (offset_a(0) + TOTAL(c_l_a(*, 0)*tb_algo_org(0:n_f-1, i)))) 
    IF ret_par.surf EQ 'in_surface' THEN BEGIN
     predictand(i) = offset_a(0) + TOTAL(c_l_a(*, 0)*tb_algo(0:n_f-1, i)) + TOTAL(c_s_a(*, 0)*surf(*, i))
     dpredictand(i) = -1.*(predictand(i) - (offset_a(0) + TOTAL(c_l_a(*, 0)*tb_algo_org(0:n_f-1, i)) + TOTAL(c_s_a(*, 0)*surf(*, i))))
    ENDIF
   ENDIF
   IF ret_par.reg EQ 'quadratic' THEN BEGIN
    predictand(i) = offset_a(0) + TOTAL(c_l_a(*, 0)*tb_algo(0:n_f-1, i)) + TOTAL(c_q_a(*, 0)*(tb_algo(0:n_f-1, i)^2))
    dpredictand(i) = -1.*(predictand(i) - (offset_a(0) + TOTAL(c_l_a(*, 0)*tb_algo_org(0:n_f-1, i)) + TOTAL(c_q_a(*, 0)*(tb_algo_org(0:n_f-1, i)^2))))
    IF ret_par.surf EQ 'in_surface' THEN BEGIN
     predictand(i) = offset_a(0) + TOTAL(c_l_a(*, 0)*tb_algo(0:n_f-1, i)) + TOTAL(c_q_a(*, 0)*(tb_algo(0:n_f-1, i)^2)) + TOTAL(c_s_a(*, 0)*surf(*, i))
     dpredictand(i) = -1.*(predictand(i) - (offset_a(0) + TOTAL(c_l_a(*, 0)*tb_algo_org(0:n_f-1, i)) + TOTAL(c_q_a(*, 0)*(tb_algo_org(0:n_f-1, i)^2)) + TOTAL(c_s_a(*, 0)*surf(*, i))))
    ENDIF
   ENDIF
  ENDIF ELSE BEGIN
   IF verbose AND go EQ 0 THEN BEGIN
    print, 'Warning RET_2A_ARBANG: retrievals not applied to all cases due to angle mismatch'
    go = 1
   ENDIF
  ENDELSE
 ENDELSE

ENDFOR
NO_RETRIEVAL:

END
