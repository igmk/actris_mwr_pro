;+
;************
PRO GET_RPG,$
;************
;INPUT
raw_path,$         ;path to raw data this day
raw_path_prev,$    ;path to raw data previous day
ang_low,$          ;plotting range elevation (low)
ang_high,$         ;plotting range elevation (high)
rd_file,$          ;radiometer data type specifier
date,$             ;date STRING (yyyymmdd)
;OUTPUT
tb_struc,$         ;tb variables (structure): {time:time_tb, tb:tb, el:el_tb, az:az_tb, f:freq, r:rain_tb}
plot_struc,$       ;variables to be plotetd (structure): {time:time_tb_p, tb:tb_p, r:rain_tb_p}
l1c_struc,$        ;elevation scanning tb variable (structure): {time:time_l1c, tb:tb_l1c, el:el_l1c, az:az_l1c, f:freq_blb, r:rain_0}
met_struc,$        ;environmental variables (structure): {time:time_m, p:pres_m, t:temp_m, q:humi_m, r:rain_m}
hkd_struc,$        ;receiver sanity checks for HKD files (structure): {time:time_hkd, rec_san:rec_sanity}  
irt_struc,$        ;IRT TBs (structure): {wavel:wavel_irt, time:time_irt, tb:tb_irt, el:el_irt, az:az_irt}
;KEYWORDS
verbose=verbose

; Abstract:
;* read raw RPG radiometer data and put them into daily structures
; Author:
; U. Loehnert
; Date:
; 2011-02-10
; Dependencies:
; -
; Changes:
; *20120518 (UL) BUG FIXED
;  Up to now, the data from the first file of each raw-data file type was stored twice in the output structure.
; *20120522 (UL) BUG FIXES
;  Time stamps for IRT and HKD files were wrong. Time from MET files was erroneously used.
; *20120725 (UL)
;  If no data is avaialble from certain file types MET, HKD or IRT for "date", the structure entries are set back to -999.
; *20120930 (UL)
;  dimension of tb_l1c set to dummy values if no BLB files available
; *20130415 (UL) BUG FIXED
;  Indices for BLB-TBs were set incorrectly; Up to now only BLB files with one BLB scan were used correctly
; *20130806 (UL) BUG FIXED
;  Indices for BLB-TBs were still set incorrectly...
;-

;**extract day from date string
date_d = FIX(STRMID(date, 4, 2))

IF rd_file EQ 'wvl' THEN BEGIN
;**search and count wvl files (wvl files are RPG specific)
 in_name_wvl_today = FILE_SEARCH(raw_path + '*.WVL')
 IF in_name_wvl_today(0) EQ '' THEN in_name_wvl_today = FILE_SEARCH(raw_path + '*.wvl')

 in_name_wvl_prev = FILE_SEARCH(raw_path_prev + '*.WVL')
 IF in_name_wvl_prev(0) EQ '' THEN in_name_wvl_prev = FILE_SEARCH(raw_path_prev + '*.wvl')

 in_name_wvl = [in_name_wvl_prev, in_name_wvl_today]
 IF in_name_wvl_prev(0) EQ '' THEN in_name_wvl = [in_name_wvl_today]
 IF in_name_wvl_today(0) EQ '' THEN in_name_wvl = [in_name_wvl_prev]
 n_in_wvl = N_ELEMENTS(in_name_wvl)
 IF in_name_wvl(0) EQ '' THEN n_in_wvl = 0

 IF verbose then print, 'total number of WVL files found: ', n_in_wvl

;**extract data from wvl files
 go = 0
 FOR ifile = 0, n_in_wvl-1 DO BEGIN

  GET_SPEC, verbose, in_name_wvl(ifile), n_wvl, time_wvl, rain_wvl, f_wvl, tb_wvl, az_wvl, el_wvl

;*angle range for quicklooks
  i_ang = WHERE(el_wvl GT ang_low AND el_wvl LT ang_high)

  IF go EQ 0 AND n_wvl GT 1 THEN BEGIN

   time_tbw = time_wvl
   tbw = tb_wvl
   rain_tbw = rain_wvl
   az_tbw = az_wvl
   el_tbw = el_wvl

   IF i_ang(0) NE -1 THEN BEGIN
    time_tbw_p = time_wvl(i_ang)
    tbw_p = tb_wvl(*, i_ang)
    rain_tbw_p = rain_wvl(i_ang)
   ENDIF

  ENDIF ELSE IF go EQ 1 AND n_wvl GT 1 THEN BEGIN

   time_tbw = [time_tbw, time_wvl]
   tbw = [[tbw], [tb_wvl]]
   rain_tbw = [rain_tbw, rain_wvl]
   az_tbw = [az_tbw, az_wvl] 
   el_tbw = [el_tbw, el_wvl]

   IF i_ang(0) NE -1 THEN BEGIN  
    time_tbw_p = [time_tbw_p, time_wvl(i_ang)]
    tbw_p = [[tbw_p], [tb_wvl(*, i_ang)]]
    rain_tbw_p = [rain_tbw_p, rain_wvl(i_ang)]   
   ENDIF

  ENDIF

  IF n_wvl GT 1 THEN go = 1
 ENDFOR ; loop over n_in_wvl

;**convert times
 time_tbw = time_tbw/86400d
 time_tbw = time_tbw + JULDAY(1,1,2001,0,0,0)
 i_sort = SORT(time_tbw)
 time_tbw = time_tbw(i_sort)
 rain_tbw = rain_tbw(i_sort)
 tbw = tbw(*, i_sort)
 az_tbw = az_tbw(i_sort)
 el_tbw = el_tbw(i_sort)

 CALDAT, time_tbw, month_tbw, day_tbw, year_tbw, hour_tbw, min_tbw, sec_tbw

 i_day = WHERE(day_tbw EQ date_d)
 time_tbw = hour_tbw(i_day) + min_tbw(i_day)/60d + sec_tbw(i_day)/3600d
 rain_tbw = rain_tbw(i_day)
 tbw = tbw(*, i_day)
 az_tbw = az_tbw(i_day)
 el_tbw = el_tbw(i_day)
 nw = LONG(N_ELEMENTS(time_tbw))

;**convert times for arrays designated for plotting
 time_tbw_p = time_tbw_p/86400d
 time_tbw_p = time_tbw_p + JULDAY(1,1,2001,0,0,0)
 i_sort = SORT(time_tbw_p)
 time_tbw_p = time_tbw_p(i_sort)
 rain_tbw_p = rain_tbw_p(i_sort)
 tbw_p = tbw_p(*, i_sort)

 CALDAT, time_tbw_p, month_tbw_p, day_tbw_p, year_tbw_p, hour_tbw_p, min_tbw_p, sec_tbw_p

 i_day = WHERE(day_tbw_p EQ date_d)
 time_tbw_p = hour_tbw_p(i_day) + min_tbw_p(i_day)/60d + sec_tbw_p(i_day)/3600d
 rain_tbw_p = rain_tbw_p(i_day)
 tbw_p = tbw_p(*, i_day)
 nw_p = LONG(N_ELEMENTS(time_tbw_p))

;**search and count olc files (olc files are RPG specific)
 in_name_olc_today = FILE_SEARCH(raw_path + '*.OLC')
 IF in_name_olc_today(0) EQ '' THEN in_name_olc_today = FILE_SEARCH(raw_path + '*.olc')

 in_name_olc_prev = FILE_SEARCH(raw_path_prev + '*.OLC')
 IF in_name_olc_prev(0) EQ '' THEN in_name_olc_prev = FILE_SEARCH(raw_path_prev + '*.olc')

 in_name_olc = [in_name_olc_prev, in_name_olc_today]
 IF in_name_olc_prev(0) EQ '' THEN in_name_olc = [in_name_olc_today]
 IF in_name_olc_today(0) EQ '' THEN in_name_olc = [in_name_olc_prev]
 n_in_olc = N_ELEMENTS(in_name_olc)
 IF in_name_olc(0) EQ '' THEN n_in_olc = 0
 
 IF verbose then print, 'total number of OLC files found: ', n_in_olc

;**extract data from olc files

 go = 0
 FOR ifile = 0, n_in_olc-1 DO BEGIN
  GET_SPEC, verbose, in_name_olc(ifile), n_olc, time_olc, rain_olc, f_olc, tb_olc, az_olc, el_olc

;*angle range for quicklooks
  i_ang = WHERE(el_olc GT ang_low AND el_olc LT ang_high)

  IF go EQ 0 AND n_olc GT 1 THEN BEGIN
   time_tbo = time_olc
   tbo = tb_olc
   rain_tbo = rain_olc
   az_tbo = az_olc
   el_tbo = el_olc

   IF i_ang(0) NE -1 THEN BEGIN
    time_tbo_p = time_olc(i_ang)
    tbo_p = tb_olc(*, i_ang)
    rain_tbo_p = rain_olc(i_ang)
   ENDIF

  ENDIF ELSE IF go EQ 1 AND n_olc GT 1 THEN BEGIN
   time_tbo = [time_tbo, time_olc]
   tbo = [[tbo], [tb_olc]]
   rain_tbo = [rain_tbo, rain_olc]
   az_tbo = [az_tbo, az_olc] 
   el_tbo = [el_tbo, el_olc]

   IF i_ang(0) NE -1 THEN BEGIN  
    time_tbo_p = [time_tbo_p, time_olc(i_ang)]
    tbo_p = [[tbo_p], [tb_olc(*, i_ang)]]
    rain_tbo_p = [rain_tbo_p, rain_olc(i_ang)]   
   ENDIF

  ENDIF
  IF n_olc GT 1 THEN go = 1
 ENDFOR ; loop over n_in_olc

;**convert times
 time_tbo = time_tbo/86400d
 time_tbo = time_tbo + JULDAY(1,1,2001,0,0,0)
 i_sort = SORT(time_tbo)
 time_tbo = time_tbo(i_sort)
 rain_tbo = rain_tbo(i_sort)
 tbo = tbo(*, i_sort)
 az_tbo = az_tbo(i_sort)
 el_tbo = el_tbo(i_sort)

 CALDAT, time_tbo, month_tbo, day_tbo, year_tbo, hour_tbo, min_tbo, sec_tbo

 i_day = WHERE(day_tbo EQ date_d)
 time_tbo = hour_tbo(i_day) + min_tbo(i_day)/60d + sec_tbo(i_day)/3600d
 rain_tbo = rain_tbo(i_day)
 tbo = tbo(*, i_day)
 az_tbo = az_tbo(i_day)
 el_tbo = el_tbo(i_day)
 no = N_ELEMENTS(time_tbo)

;**convert times for arrays designated for plotting
 time_tbo_p = time_tbo_p/86400d
 time_tbo_p = time_tbo_p + JULDAY(1,1,2001,0,0,0)
 i_sort = SORT(time_tbo_p)
 time_tbo_p = time_tbo_p(i_sort)
 rain_tbo_p = rain_tbo_p(i_sort)
 tbo_p = tbo_p(*, i_sort)

 CALDAT, time_tbo_p, month_tbo_p, day_tbo_p, year_tbo_p, hour_tbo_p, min_tbo_p, sec_tbo_p

 i_day = WHERE(day_tbo_p EQ date_d)
 time_tbo_p = hour_tbo_p(i_day) + min_tbo_p(i_day)/60d + sec_tbo_p(i_day)/3600d
 rain_tbo_p = rain_tbo_p(i_day)
 tbo_p = tbo_p(*, i_day)
 no_p = N_ELEMENTS(time_tbo_p)

;**put olc data on time_tbw
 time_tb = time_tbw 
 nf_wvl = N_ELEMENTS(f_wvl)
 nf_olc = N_ELEMENTS(f_olc)
 tb = REPLICATE(-999., nf_wvl+nf_olc, nw)
 rain_tb = rain_tbw
 el_tb = el_tbw
 az_tb = az_tbw

;*look for nearest neighbor tbw/tbo
 FOR iw = 0l, nw-1l DO BEGIN
  tb(0:nf_wvl-1, iw) = tbw(*, iw)
  CLOSEST, time_tb(iw), time_tbo, dmin, ind
  IF dmin LT 2d/3600d THEN tb(nf_wvl:nf_wvl+nf_olc-1, iw) = tbo(*, ind)
 ENDFOR 

;**put olc data on time_tbw (plotting designated)
 time_tb_p = time_tbw_p 
 tb_p = REPLICATE(-999., nf_wvl+nf_olc, nw_p)
 rain_tb_p = rain_tbw_p

;*look for nearest neighbor tbw/tbo
 FOR iw = 0l, LONG(nw_p)-1l DO BEGIN
  tb_p(0:nf_wvl-1, iw) = tbw_p(*, iw)
  CLOSEST, time_tb_p(iw), time_tbo_p, dmin, ind
  IF dmin LT 2d/3600d THEN tb_p(nf_wvl:nf_wvl+nf_olc-1, iw) = tbo_p(*, ind)
 ENDFOR 
 plot_struc = {time:time_tb_p, tb:tb_p, r:rain_tb_p}
 
 freq = [f_wvl, f_olc]

ENDIF ; wvl/olc

IF rd_file EQ 'brt' THEN BEGIN

;**search and count brt files
;(brt files are RPG specific)

 in_name_brt_today = FILE_SEARCH(raw_path + '*.BRT')
 IF in_name_brt_today(0) EQ '' THEN in_name_brt_today = FILE_SEARCH(raw_path + '*.brt')

 in_name_brt_prev = FILE_SEARCH(raw_path_prev + '*.BRT')
 IF in_name_brt_prev(0) EQ '' THEN in_name_brt_prev = FILE_SEARCH(raw_path_prev + '*.brt')

 in_name_brt = [in_name_brt_prev, in_name_brt_today]
 IF in_name_brt_prev(0) EQ '' THEN in_name_brt = [in_name_brt_today]
 IF in_name_brt_today(0) EQ '' THEN in_name_brt = [in_name_brt_prev]
 n_in_brt = N_ELEMENTS(in_name_brt)
 IF in_name_brt(0) EQ '' THEN n_in_brt = 0

 IF verbose THEN print, 'total number of BRT files found: ', n_in_brt

;**extract data from brt files

 go = 0
 FOR ifile = 0, n_in_brt-1 DO BEGIN
  GET_BRT, in_name_brt(ifile), n_brt, time_brt, rain_brt, freq, tb_brt, az_brt, el_brt, verbose=0
;*angle range for quicklooks
  i_ang = WHERE(el_brt GT ang_low AND el_brt LT ang_high)

  IF go EQ 0 AND n_brt GT 1 THEN BEGIN

   time_tb = time_brt
   tb = tb_brt
   rain_tb = rain_brt
   az_tb = az_brt
   el_tb = el_brt

   IF i_ang(0) NE -1 THEN BEGIN
    time_tb_p = time_brt(i_ang)
    tb_p = tb_brt(*, i_ang)
    rain_tb_p = rain_brt(i_ang)
   ENDIF

   go = 1

  ENDIF ELSE IF go EQ 1 AND n_brt GT 1 THEN BEGIN

   time_tb = [time_tb, time_brt]
   tb = [[tb], [tb_brt]]
   rain_tb = [rain_tb, rain_brt]
   az_tb = [az_tb, az_brt] 
   el_tb = [el_tb, el_brt]
   IF i_ang(0) NE -1 THEN BEGIN  
    IF N_ELEMENTS(time_tb_p) GT 0 THEN BEGIN
     time_tb_p = [time_tb_p, time_brt(i_ang)]
     tb_p = [[tb_p], [tb_brt(*, i_ang)]]
     rain_tb_p = [rain_tb_p, rain_brt(i_ang)]   
    ENDIF ELSE BEGIN
     time_tb_p = time_brt(i_ang)
     tb_p = tb_brt(*, i_ang)
     rain_tb_p = rain_brt(i_ang)
    ENDELSE
   ENDIF

  ENDIF

 ENDFOR

;**convert times
 time_tb = time_tb/86400d
 time_tb = time_tb + JULDAY(1,1,2001,0,0,0)
 i_sort = SORT(time_tb)
 time_tb = time_tb(i_sort)
 rain_tb = rain_tb(i_sort)
 tb = tb(*, i_sort)
 az_tb = az_tb(i_sort)
 el_tb = el_tb(i_sort)
 CALDAT, time_tb, month_tb, day_tb, year_tb, hour_tb, min_tb, sec_tb
 i_day = WHERE(day_tb EQ date_d)

 time_tb = hour_tb(i_day) + min_tb(i_day)/60d + sec_tb(i_day)/3600d
 rain_tb = rain_tb(i_day)
 tb = tb(*, i_day)
 az_tb = az_tb(i_day)
 el_tb = el_tb(i_day)

;**convert times for arrays designated for plotting
 IF N_ELEMENTS(time_tb_p) GT 1 THEN BEGIN
  time_tb_p = time_tb_p/86400d
  time_tb_p = time_tb_p + JULDAY(1,1,2001,0,0,0)
  i_sort = SORT(time_tb_p)
  time_tb_p = time_tb_p(i_sort)
  rain_tb_p = rain_tb_p(i_sort)
  tb_p = tb_p(*, i_sort)

  CALDAT, time_tb_p, month_tb_p, day_tb_p, year_tb_p, hour_tb_p, min_tb_p, sec_tb_p

  i_day = WHERE(day_tb_p EQ date_d)
  time_tb_p = hour_tb_p(i_day) + min_tb_p(i_day)/60d + sec_tb_p(i_day)/3600d
  rain_tb_p = rain_tb_p(i_day)
  tb_p = tb_p(*, i_day)
  plot_struc = {time:time_tb_p, tb:tb_p, r:rain_tb_p}
 ENDIF
ENDIF ; brt
tb_struc = {time:time_tb, tb:tb, el:el_tb, az:az_tb, f:freq, r:rain_tb}

IF rd_file EQ 'brt' OR rd_file EQ 'wvl' THEN BEGIN

;**search and count blb files
;(blb files are RPG specific)

 in_name_blb_today = FILE_SEARCH(raw_path + '*.BLB')
 IF in_name_blb_today(0) EQ '' THEN in_name_blb_today = FILE_SEARCH(raw_path + '*.blb')

 in_name_blb_prev = FILE_SEARCH(raw_path_prev + '*.BLB')
 IF in_name_blb_prev(0) EQ '' THEN in_name_blb_prev = FILE_SEARCH(raw_path_prev + '*.blb')

 in_name_blb = [in_name_blb_prev, in_name_blb_today]
 IF in_name_blb_prev(0) EQ '' THEN in_name_blb = [in_name_blb_today]
 IF in_name_blb_today(0) EQ '' THEN in_name_blb = [in_name_blb_prev]
 n_in_blb = N_ELEMENTS(in_name_blb)
 IF in_name_blb(0) EQ '' THEN n_in_blb = 0

 IF verbose THEN print,'total number of BLB files found: ', n_in_blb

;**extract data from blb files

 go = 0
 freq_blb = -999.
 ang_el_blb = -999. 
 n_all_blb = 0

 FOR ifile = 0, n_in_blb-1 DO BEGIN
  GET_BLB, in_name_blb(ifile), n_blb, time_blb, rain_blb, freq_blb, na_blb, tb_blb, el_blb, az_blb, verbose=0

  n_all_blb = n_blb + n_all_blb
;*find maximum dimension of el_blb for one day
  IF ifile EQ 0 THEN BEGIN
   ang_el_blb = el_blb
   ang_az_blb = az_blb
  ENDIF ELSE BEGIN
   n_ang_el_blb = N_ELEMENTS(el_blb)
   FOR k = 0, n_ang_el_blb-1 DO BEGIN
    k_ind = WHERE(el_blb(k) EQ ang_el_blb)
    IF k_ind(0) EQ -1 THEN BEGIN
     ang_el_blb = [ang_el_blb, el_blb(k)]     
     ang_az_blb = [ang_az_blb, az_blb(k)]     
     isort = SORT(ang_el_blb)
     ang_el_blb = ang_el_blb(isort)
     isort = SORT(ang_az_blb)
     ang_az_blb = ang_az_blb(isort)
    ENDIF
   ENDFOR
  ENDELSE 
 ENDFOR

 n_all_blb_x = n_all_blb
 IF n_all_blb EQ 0 THEN n_all_blb_x = 1
 tb_l1c = REPLICATE(-999., N_ELEMENTS(freq_blb), N_ELEMENTS(ang_el_blb), n_all_blb_x)

 FOR ifile = 0, n_in_blb-1 DO BEGIN

  GET_BLB, in_name_blb(ifile), n_blb, time_blb, rain_blb, freq_blb, na_blb, tb_blb, el_blb, az_blb, verbose=0

  IF in_name_blb(ifile) EQ '' THEN n_blb = 0
;tb_blb: (freq x ang x time)

  IF go EQ 0 AND n_blb GT 0 THEN BEGIN
   time_l1c = time_blb
   n_el_blb = N_ELEMENTS(el_blb)

   FOR k = 0, n_el_blb-1 DO BEGIN
    k_ind = WHERE(el_blb(k) EQ ang_el_blb)
    tb_l1c(*, k_ind(0), 0:n_blb-1) = tb_blb(*, k, *)        
   ENDFOR
   n_blb_old = n_blb
 
   rain_l1c = rain_blb
   az_l1c = ang_az_blb ; assuming that az_blb stay constant over the day
   el_l1c = ang_el_blb

   go = 1

  ENDIF ELSE IF go EQ 1 AND n_blb GT 0 THEN BEGIN

   time_l1c = [time_l1c, time_blb]
   n_el_blb = N_ELEMENTS(el_blb)
   rain_l1c = [rain_l1c, rain_blb]

   FOR k = 0, n_el_blb-1 DO BEGIN
    k_ind = WHERE(el_blb(k) EQ ang_el_blb)

    tb_l1c(*, k_ind(0), n_blb_old:n_blb+n_blb_old-1) = tb_blb(*, k, *)        
   ENDFOR
   n_blb_old = n_blb + n_blb_old

  ENDIF

 ENDFOR ; loop over ifile

;**convert times
 IF n_in_blb GT 0 THEN BEGIN
  time_l1c = time_l1c/86400d
  time_l1c = time_l1c + JULDAY(1,1,2001,0,0,0)
  i_sort = SORT(time_l1c)
  time_l1c = time_l1c(i_sort)
  rain_l1c = rain_l1c(i_sort)
  tb_l1c = tb_l1c(*, *, i_sort)
; az_l1c = az_l1c(*, i_sort)
; el_l1c = el_l1c(*, i_sort)

  CALDAT, time_l1c, month_l1c, day_l1c, year_l1c, hour_l1c, min_l1c, sec_l1c

  i_day = WHERE(day_l1c EQ date_d)
  IF i_day(0) NE -1 THEN BEGIN
   time_l1c = hour_l1c(i_day) + min_l1c(i_day)/60d + sec_l1c(i_day)/3600d
   rain_l1c = rain_l1c(i_day)
   tb_l1c = tb_l1c(*, *, i_day)
;  az_l1c = az_l1c(*, i_day)
;  el_l1c = el_l1c(*, i_day)
  ENDIF ELSE BEGIN
   time_l1c = -999.
   rain_l1c = -999.
   tb_l1c = -999.  
  ENDELSE

  l1c_struc = {time:time_l1c, tb:tb_l1c, el:el_l1c, az:az_l1c, f:freq_blb, r:rain_l1c}
 ENDIF

;**search and count met files
;(met files are RPG specific)

 in_name_met_today = FILE_SEARCH(raw_path + '*.MET')
 IF in_name_met_today(0) EQ '' THEN in_name_met_today = FILE_SEARCH(raw_path + '*.met')

 in_name_met_prev = FILE_SEARCH(raw_path_prev + '*.MET')
 IF in_name_met_prev(0) EQ '' THEN in_name_met_prev = FILE_SEARCH(raw_path_prev + '*.met')

 in_name_met = [in_name_met_prev, in_name_met_today]
 IF in_name_met_prev(0) EQ '' THEN in_name_met = [in_name_met_today]
 IF in_name_met_today(0) EQ '' THEN in_name_met = [in_name_met_prev]
 n_in_met = N_ELEMENTS(in_name_met)
 IF in_name_met(0) EQ '' THEN n_in_met = 0

 IF verbose then print, 'total number of MET files found: ', n_in_met

 go = 0
 FOR ifile = 0, n_in_met-1 DO BEGIN
  GET_MET, in_name_met(ifile), time_met, rain_met, temp_met, pres_met, humi_met, n_met, verbose=0

  IF n_met GT 1 AND go EQ 0 THEN BEGIN
   time_m = time_met
   temp_m = temp_met
   pres_m = pres_met
   humi_m = humi_met
   rain_m = rain_met
  ENDIF ELSE IF n_met GT 1 AND go EQ 1 THEN BEGIN
   time_m = [time_m, time_met]
   temp_m = [temp_m, temp_met]
   pres_m = [pres_m, pres_met]
   humi_m = [humi_m, humi_met]
   rain_m = [rain_m, rain_met]
  ENDIF
  IF n_met GT 1 THEN go = 1
 ENDFOR

;**convert times
 IF n_in_met GT 0 THEN BEGIN
  time_m = time_m/86400d
  time_m = time_m + JULDAY(1,1,2001,0,0,0)
  i_sort = SORT(time_m)
  time_m = time_m(i_sort)
  pres_m = pres_m(i_sort)
  temp_m = temp_m(i_sort)
  humi_m = humi_m(i_sort)
  rain_m = rain_m(i_sort)

  CALDAT, time_m, month_m, day_m, year_m, hour_m, min_m, sec_m
  i_day = WHERE(day_m EQ date_d, n_day)

  IF n_day GE 1 THEN BEGIN
   time_m = hour_m(i_day) + min_m(i_day)/60d + sec_m(i_day)/3600d
   pres_m = pres_m(i_day)
   temp_m = temp_m(i_day)
   humi_m = humi_m(i_day)
   rain_m = rain_m(i_day)
  ENDIF ELSE BEGIN
   time_m = -999.
   pres_m = -999. 
   temp_m = -999. 
   humi_m = -999. 
   rain_m = -999.
  ENDELSE

  met_struc = {time:time_m, p:pres_m, t:temp_m, q:humi_m, r:rain_m}

 ENDIF

;**search and count hkd files
;(hkd files are RPG specific)

 in_name_hkd_today = FILE_SEARCH(raw_path + '*.HKD')
 IF in_name_hkd_today(0) EQ '' THEN in_name_hkd_today = FILE_SEARCH(raw_path + '*.hkd')

 in_name_hkd_prev = FILE_SEARCH(raw_path_prev + '*.HKD')
 IF in_name_hkd_prev(0) EQ '' THEN in_name_hkd_prev = FILE_SEARCH(raw_path_prev + '*.hkd')

 in_name_hkd = [in_name_hkd_prev, in_name_hkd_today]
 IF in_name_hkd_prev(0) EQ '' THEN in_name_hkd = [in_name_hkd_today]
 IF in_name_hkd_today(0) EQ '' THEN in_name_hkd = [in_name_hkd_prev]
 n_in_hkd = N_ELEMENTS(in_name_hkd)
 IF in_name_hkd(0) EQ '' THEN n_in_hkd = 0

 IF verbose then print, 'total number of HKD files found: ', n_in_hkd

 go = 0
 FOR ifile = 0, n_in_hkd-1 DO BEGIN
  GET_HKD, in_name_hkd(ifile), tb_struc.f, time_hkd, rec_sanity, verbose=0
  n_hkd = N_ELEMENTS(time_hkd)

  IF n_in_hkd GT 0 AND go EQ 0 THEN BEGIN
   time_s = time_hkd
   rec_s = rec_sanity
  ENDIF ELSE IF n_in_hkd GT 0 AND go EQ 1 THEN BEGIN
   time_s = [time_s, time_hkd]
   rec_s = [[rec_s], [rec_sanity]]
  ENDIF
  IF n_hkd GT 1 THEN go = 1
 ENDFOR

;**convert times
 IF n_in_hkd GT 0 THEN BEGIN
  time_s = time_s/86400d
  time_s = time_s + JULDAY(1,1,2001,0,0,0)
  i_sort = SORT(time_s)
  time_s = time_s(i_sort)
  rec_s = rec_s(*, i_sort)

  CALDAT, time_s, month_s, day_s, year_s, hour_s, min_s, sec_s
  i_day = WHERE(day_s EQ date_d, n_day)

  IF n_day GE 1 THEN BEGIN
   time_s = hour_s(i_day) + min_s(i_day)/60d + sec_s(i_day)/3600d
   rec_s = rec_s(*, i_day)
  ENDIF ELSE BEGIN
   time_s = -999.
   rec_s = -999.
  ENDELSE

  hkd_struc = {time:time_s, rec_san:rec_s}

 ENDIF

;**search and count irt files
;(irt files are RPG specific)

 in_name_irt_today = FILE_SEARCH(raw_path + '*.IRT')
 IF in_name_irt_today(0) EQ '' THEN in_name_irt_today = FILE_SEARCH(raw_path + '*.irt')

 in_name_irt_prev = FILE_SEARCH(raw_path_prev + '*.IRT')
 IF in_name_irt_prev(0) EQ '' THEN in_name_irt_prev = FILE_SEARCH(raw_path_prev + '*.irt')

 in_name_irt = [in_name_irt_prev, in_name_irt_today]
 IF in_name_irt_prev(0) EQ '' THEN in_name_irt = [in_name_irt_today]
 IF in_name_irt_today(0) EQ '' THEN in_name_irt = [in_name_irt_prev]
 n_in_irt = N_ELEMENTS(in_name_irt)
 IF in_name_irt(0) EQ '' THEN n_in_irt = 0

 IF verbose then print, 'total number of IRT files found: ', n_in_irt

 go = 0
 FOR ifile = 0, n_in_irt-1 DO BEGIN
  GET_IRT, in_name_irt(ifile), time_irt, wavel_irt, tb_irt, el_irt, az_irt, verbose=0
  n_irt = N_ELEMENTS(time_irt)

  IF n_in_irt GT 0 AND go EQ 0 THEN BEGIN
   time_s = time_irt
   tb_irt_s = tb_irt
   el_irt_s = el_irt
   az_irt_s = az_irt
  ENDIF ELSE IF n_in_irt GT 0 AND go EQ 1 THEN BEGIN
   time_s = [time_s, time_irt]
   tb_irt_s = [[tb_irt_s], [tb_irt]]
   el_irt_s = [el_irt_s, el_irt]
   az_irt_s = [az_irt_s, az_irt]
  ENDIF
  IF n_irt GT 1 THEN go = 1
 ENDFOR

;**convert times
 IF n_in_irt GT 0 THEN BEGIN
  time_s = time_s/86400d
  time_s = time_s + JULDAY(1,1,2001,0,0,0)
  i_sort = SORT(time_s)
  time_s = time_s(i_sort)
  tb_irt_s = tb_irt_s(*, i_sort)
  el_irt_s = el_irt_s(i_sort)
  az_irt_s = az_irt_s(i_sort)

  CALDAT, time_s, month_s, day_s, year_s, hour_s, min_s, sec_s
  i_day = WHERE(day_s EQ date_d, n_day)

  IF n_day GE 1 THEN BEGIN  
   time_s = hour_s(i_day) + min_s(i_day)/60d + sec_s(i_day)/3600d
   tb_irt_s = tb_irt_s(*, i_day)
   el_irt_s = el_irt_s(i_day)
   az_irt_s = az_irt_s(i_day) 
  ENDIF ELSE BEGIN
   time_s = -999.
   tb_irt_s = -999.
   el_irt_s = -999.
   az_irt_s = -999.
  ENDELSE

  irt_struc = {wavel:wavel_irt, time:time_s, tb:tb_irt_s, el:el_irt_s, az:az_irt_s}

 ENDIF

ENDIF

END
