;+
;************
PRO GET_MET,$
;************
;INPUT
filename,$
;OUPUT
time,$
rain,$
temp,$
pres,$
humi,$
n,$
;KEYWORDS
verbose=verbose                         
; Abstract:
;* read raw RPG .MET files
; Author:
; U. Loehnert
; Date:
; 2012-07-25
; Dependencies:
; -
; Changes:
; XXXX-XX-XX: ???
;-

IF N_ELEMENTS(verbose) EQ 0 THEN verbose = 0

IF verbose THEN PRINT, 'GET_MET; filename=', filename
IF filename EQ '' THEN BEGIN
  n = 0
  RETURN
ENDIF
ON_IOERROR, SPEC_END
OPENR,unit,filename,/GET_LUN
code=1l
READU,unit,code

IF code NE 599658943 THEN BEGIN
 IF code NE 599658944 THEN BEGIN
  PRINT,'Error in MET file code'
  GOTO, SPEC_END
 ENDIF
ENDIF


n=0l
READU,unit,n

IF n GT 86400 THEN BEGIN
 n = 0
 GOTO, SPEC_END
ENDIF
IF verbose THEN print,'No of samples=',n
IF n LT 1 THEN GOTO,SPEC_END

n_add = 0b ; number of additional sensors
IF code EQ 599658944 THEN BEGIN
 READU,unit,n_add
 IF verbose THEN print,n_add
ENDIF

xmin = -999.e0 ; a 4Byte floating point / single precision default value
xmax = -999.e0 ; a 4Byte floating point / single precision default value
FOR i =0,2 DO BEGIN
  READU,unit,xmin,xmax
  IF verbose THEN print,xmin,xmax
ENDFOR

time_ref = 0l
READU,unit,time_ref
IF verbose THEN BEGIN
CASE time_ref of
  1: print,'UTC'
  0: print,'Local Time'
  ELSE: print,'No valid time'
  ENDCASE
ENDIF

time=LONARR(n)
rain=BYTARR(n)
temp=FLTARR(n)
pres=FLTARR(n)
humi=FLTARR(n)
x1 = 0l ; a 4Byte integer / Long default value
x2 = 0b ; a 1Byte default value
x3 = 1E1 ; single precision / 4Byte floating point
x4 = 1E1 ; single precision / 4Byte floating point
x5 = 1E1 ; single precision / 4Byte floating point
FOR i =0l,LONG(n-1l) DO BEGIN
  READU,unit,x1,x2,x3,x4,x5
  time(i) = x1
  rain(i) = x2
  temp(i) = x4
  pres(i) = x3
  humi(i) = x5
  IF verbose THEN PRINT,i,x1,x2,x3,x4,x5
ENDFOR

SPEC_END:
;print,n
CLOSE,unit
FREE_LUN,unit


END
