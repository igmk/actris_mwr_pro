;+
;************
PRO GET_RES,$
;************
;INPUT
raw_path,$         ;path to raw data this day
raw_path_prev,$    ;path to raw data previous day
ang_low,$          ;plotting range elevation (low)
ang_high,$         ;plotting range elevation (high)
date,$             ;yyyymmdd (STRING)
;OUTPUT
tb,$               ;tb variables (structure): {time:time_tb, tb:tb, el:el_tb}
plot,$             ;variables to be plotetd (structure): {time:time_tb_p, tb:tb_p}
;KEYWORDS
verbose = verbose

; Abstract:
;* read raw RESCOM radiometer data and put them into daily structures
; Author:
; U. Loehnert
; Date:
; 2011-02-10
; Dependencies:
; -
; Changes:
; *20120518 (UL) BUG FIXED
;  Up to now, the data from the first file of each raw-data file type was stored twice in the output structure.
;-

;**initialize output with dummy values
tb = {time:-999., tb:-999., el:-999.}
plot = {time:-999., tb:-999., r:-999.}

;**set day number of month
date_d = FIX(STRMID(date, 4, 2))

;**search and count dec files
;(dec files are RESCOM specific)

in_name_dec = FILE_SEARCH(raw_path +'*'+date+'*.dec', /FOLD_CASE)
n_in_dec = N_ELEMENTS(in_name_dec)
IF in_name_dec(0) EQ '' THEN n_in_dec = 0

IF verbose THEN BEGIN
 print, 'searching in: ', raw_path
 print, 'total number of DEC files found: ', n_in_dec
ENDIF

time_tb = -999.
tb = REPLICATE(-999., 1, 1)
el_tb = -999.
time_tb_p = REPLICATE(-999., 1, 1)
tb_p = -999.

go = 0
FOR ifile = 0, n_in_dec-1 DO BEGIN
 GET_DEC, verbose, in_name_dec(ifile), time_dec, tb_dec, el_dec, freq

;*sort frequencies and TBs to increasing values
 isort_freq = SORT(freq)
 freq = freq(isort_freq)
 tb_dec = tb_dec(isort_freq, *)

;*angle range for quicklooks
 i_ang = WHERE(el_dec GT ang_low AND el_dec LT ang_high)
 IF i_ang(0) EQ -1 THEN BEGIN
  print, 'GET_RES: no elevation angles within plotting range!'
 ENDIF 

 IF go EQ 0 AND n_in_dec GE 1 THEN BEGIN
  time_tb = time_dec
  tb = tb_dec
  el_tb = el_dec

  IF i_ang(0) NE -1 THEN BEGIN
   time_tb_p = time_dec(i_ang)
   tb_p = tb_dec(*, i_ang)
  ENDIF

  go = 1
 ENDIF ELSE IF go EQ 1 AND n_in_dec GE 1 THEN BEGIN

  time_tb = [time_tb, time_dec]
  tb = [[tb], [tb_dec]]
  el_tb = [el_tb, el_dec]

  IF i_ang(0) NE -1 THEN BEGIN
   time_tb_p = [time_tb_p, time_dec(i_ang)]
   tb_p = [[tb], [tb_dec(*, i_ang)]]
  ENDIF

 ENDIF

ENDFOR

;**convert times
i_sort = SORT(time_tb)
time_tb = time_tb(i_sort)
tb = tb(*, i_sort)
el_tb = el_tb(i_sort)

CALDAT, time_tb, month_tb, day_tb, year_tb, hour_tb, min_tb, sec_tb

i_day = WHERE(day_tb EQ date_d)

time_tb = hour_tb(i_day) + min_tb(i_day)/60d + sec_tb(i_day)/3600d
tb = tb(*, i_day)
el_tb = el_tb(i_day)

;**convert times for arrays designated for plotting
i_sort = SORT(time_tb_p)
time_tb_p = time_tb_p(i_sort)
tb_p = tb_p(*, i_sort)

CALDAT, time_tb_p, month_tb_p, day_tb_p, year_tb_p, hour_tb_p, min_tb_p, sec_tb_p

i_day = WHERE(day_tb_p EQ date_d)
IF i_day(0) NE -1 THEN BEGIN
 time_tb_p = hour_tb_p(i_day) + min_tb_p(i_day)/60d + sec_tb_p(i_day)/3600d
 tb_p = tb_p(*, i_day)
ENDIF

tb = {time:time_tb, tb:tb, el:el_tb, f:freq}
plot = {time:time_tb_p, tb:tb_p, r:-999.}

END