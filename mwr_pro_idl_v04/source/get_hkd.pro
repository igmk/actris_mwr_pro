;+
;************
PRO GET_HKD,$
;************
;INPUT:
filename,$                  ;*.hkd file 
freq,$                      ;frequencies in GHz 
;OUPUT:
time,$                      ;seconds since 1.1.1970
rec_sanity,$                ;QC variable (flag)
;KEYWORDS
verbose=verbose
; Abstract:
; * reads binary *.hkd files and extracts a "receiver sanity check" for
; * radiometer quality controll. 
; * rec_sanity = 0 --> OK, rec_sanity = 1 --> internal radiometer malfunction 
; * The dimension of rec_sanity is (2xn) --> discrimination between two bands
; * It is based on "status" variable:
; * Bit1-7: Band1 channel check --> 1: OK, 0: malfunction
; * Bit9-15: Band2 channel check --> 1: OK, 0: malfunction
; * Bits25&26: --> 2: no thermal stability in Band1 receiver
; * Bits27&28: --> 2: no thermal stability in Band2 receiver
; * Bit30: --> 1: hot load sensors disagree more than 0.3K
; * ---> see RPG manual (*.HKD appendix) for further details
; Authors:
; U. Loehnert
; Date:
; 2011-02-08
; Dependencies:
; -
; Changes:
; *2012-08-27 (UL)
;  freq variable added as INPUT; non-existing channels are no-longer flagged
;-

ON_ERROR, 2

OPENR, unit, filename, /GET_LUN, ERROR=err

n_freq = N_ELEMENTS(freq)

IF N_ELEMENTS(verbose) EQ 0 THEN verbose = 0
rec_sanity = REPLICATE(-999., 2)
;verbose = 1

IF err NE 0 THEN BEGIN
  n = 0
  GOTO, END_HKD
ENDIF
ON_ERROR, 0

code = 1l
READU, unit, code
IF code NE 837854832 THEN BEGIN 
 PRINT, 'Error in HKD file code'    
 GOTO, END_HKD 
ENDIF

n = 0l
READU, unit, n
IF verbose THEN print, 'Number of samples=', n
IF n LT 1 THEN GOTO, HKD_END

; Time Reference - 1 : UTC , 0 : Local Time

time_ref = 0l
READU, unit, time_ref
IF verbose THEN BEGIN
CASE time_ref of
  1: print, 'UTC'
  0: print, 'Local Time'
  ELSE: print, 'No valid time'
  ENDCASE
ENDIF

hkdsel = 0l
READU, unit, hkdsel
IF verbose THEN BEGIN
 print, 'HKDSelect: ', hkdsel 
ENDIF

time = LONARR(n)
alarm = REPLICATE(1b, n)
lon = REPLICATE(FLOAT(-999.), n)
lat = REPLICATE(FLOAT(-999.), n)
temps = REPLICATE(FLOAT(-999.), 4, n)
stab = REPLICATE(FLOAT(-999.), 2, n)
flash = REPLICATE(-99l, n)
quali = REPLICATE(-99l, n)
status = REPLICATE(-99l, n)
rec_sanity = REPLICATE(-99l, 2, n)

FOR i = 0l, LONG(n-1l) DO BEGIN
 
 l1 = 0l
 READU, unit, l1
 time(i) = l1
 IF verbose THEN print, 'time: ', l1

 s1 = 1b
 READU, unit, s1
 alarm(i) = s1
 IF verbose THEN print, 'alarm: ', s1

 x1 = FLOAT(0.)
 x2 = FLOAT(0.)

 IF (hkdsel AND 1) NE 0 THEN BEGIN
  READU, unit, x1, x2
  lon(i) = x1
  lat(i) = x2
  IF verbose THEN print, 'lon/lat: ', x1, x2
 ENDIF

 IF (hkdsel AND 2) NE 0 THEN BEGIN
  x3 = FLOAT(0.)
  x4 = FLOAT(0.)
  READU, unit, x1, x2, x3, x4
  temps(0, i) = x1
  temps(1, i) = x2
  temps(2, i) = x3
  temps(3, i) = x4   
  IF verbose THEN print, 'Temps. amb1/amb2/humpro/tempro: ', temps(*, i)
 ENDIF

 IF (hkdsel AND 4) NE 0 THEN BEGIN
  READU, unit, x1, x2
  stab(0, i) = x1
  stab(1, i) = x2
  IF verbose THEN print, 'temp. stability rec. 1/2: ', stab(*, i)
 ENDIF

 IF (hkdsel AND 8) NE 0 THEN BEGIN
  READU, unit, l1
  flash(i) = l1
  IF verbose THEN print, 'flash mem. remaining: ', l1
 ENDIF

 IF (hkdsel AND 16) NE 0 THEN BEGIN
  READU, unit, l1
  quali(i) = l1
  IF verbose THEN print, 'ret. quality: ', quali(i)
 ENDIF

 IF (hkdsel AND 32) NE 0 THEN BEGIN
  READU, unit, l1
  status(i) = l1
  IF verbose THEN print, 'status (see manual): ', status(i)
 ENDIF

;***receiver sanity check generated 
;*based on status variable
;*Bit1-7: Band1 channel check --> 1: OK, 0: malfunction
;*Bit9-15: Band2 channel check --> 1: OK, 0: malfunction
;*Bits25&26: --> 2: no thermal stability in Band1 receiver
;*Bits27&28: --> 2: no thermal stability in Band2 receiver
;*Bit30: --> 1: hot load sensors disagree more than 0.3K

 rec_sanity(*, i) = 0l
 f_counter = 0

 IF verbose THEN print, filename
 FOR j = 1l, 7l DO BEGIN
  f_counter = f_counter + 1
  IF f_counter GT n_freq THEN GOTO, SKIP
  IF (status(i) AND 2l^(j-1l)) EQ 0 THEN BEGIN
   rec_sanity(0, i) = 1l
   IF verbose THEN print, 'malfunction channel: ',j
  ENDIF
 ENDFOR

 FOR j = 9l, 15l DO BEGIN
  f_counter = f_counter + 1
  IF f_counter GE n_freq THEN GOTO, SKIP
  IF (status(i) AND 2l^(j-1l)) EQ 0 THEN BEGIN
   rec_sanity(1, i) = 1l
   IF verbose THEN print, 'malfunction channel: ',j
  ENDIF
 ENDFOR

SKIP:
 IF (status(i) AND 2l^(25l)) NE 0 THEN BEGIN
  rec_sanity(0, i) = 1l
  IF verbose THEN print, 'problem thermal stab. Band1'
 ENDIF
 IF (status(i) AND 2l^(27l)) NE 0 THEN BEGIN
  rec_sanity(1, i) = 1l
  IF verbose THEN print, 'problem thermal stab. Band2'
 ENDIF
 IF (status(i) AND 2l^(29l)) NE 0 THEN BEGIN
  rec_sanity(*, i) = 1l
  IF verbose THEN print, 'hot load targets disagree'
 ENDIF
ENDFOR

HKD_END:

CLOSE,unit
FREE_LUN,unit

END_HKD:
END
