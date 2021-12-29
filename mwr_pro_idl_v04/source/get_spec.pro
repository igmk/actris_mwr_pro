;+
;*************
PRO GET_SPEC,$
;*************
;INPUT:
druck,$                     ;talkative=1, quiet=0
filename,$                  ;*.wvl or *.olc file
;OUPUT:
n,$                         ;number of measurements
time,$                      ;seconds since 1.1.1970
rain,$                      ;rain flag (0/1)
f,$                         ;array of frequencies
tb,$                        ;TB arry (f x time)
az,$                        ;array of azimuth angles
el                          ;array of elevation angles
; Abstract:
; * reads binary *.wvl or *.olc files
; Authors:
; S. Crewell, U. Loehnert
; Date:
; 2008-09-12
; Dependencies:
; -
; Changes:
; XXXX-XX-XX:
; changed program ...
;-

ON_ERROR, 2
ON_IOERROR, END_SPEC

;*override verbose
druck = 0

OPENR,unit,filename,/GET_LUN,ERROR=err
IF err NE 0 THEN BEGIN
  n = 0
  GOTO,END_SPEC
ENDIF
ON_ERROR,0

code=1l
READU,unit,code
IF code NE 955874342 AND code NE 456783953 THEN PRINT,'Error in SPEC file code'

n=0l
READU,unit,n
IF druck THEN print,'Number of samples=',n
IF n LT 1 THEN GOTO,SPEC_END

xmin = -999.
xmax = -999.
READU,unit,xmin,xmax
IF druck THEN print,xmin,xmax

time_ref = 0l
READU,unit,time_ref
IF druck THEN BEGIN
CASE time_ref of
  1: print,'UTC'
  0: print,'Local Time'
  ELSE: print,'No valid time'
  ENDCASE
ENDIF

n_f = 0l
READU,unit,n_f
IF druck THEN PRINT,'Number of frequencies=',n_f

f = FLTARR(n_f)
READU,unit,f
IF druck THEN PRINT,'Frequencies=',f

time= LONARR(n)
rain= BYTARR(n)
tb  = FLTARR(n_f,n)
az  = FLTARR(n)
el  = FLTARR(n)
xxx = FLTARR(n)

x1 = 0l
x2 = 0b
help_var = FLTARR(n_f)
FOR i =0l,LONG(n-1l) DO BEGIN
  READU,unit,x1,x2,help_var,x
  time(i) = x1
  rain(i) = x2
  tb(*,i) = help_var
  xxx(i) = x
  sign = 1
  IF x LT 0. THEN sign = -1.
  az(i)   = sign*FLOOR(x/100.)/10.
  el(i)   = x - (sign*az(i)*1000.)
;  time_x  = STRING(time(i)/(86400.d)+JULDAY(1,1,2001,0),format='(C(Cdi02,":",Cmoi02,":",Cyi04,":",Chi02,":",Cmi02,":",Csi02))')
;  IF druck THEN print, i, time(i), rain(i), tb(1,i), x, az(i), el(i)
ENDFOR

SPEC_END:
CLOSE, unit
FREE_LUN, unit

END_SPEC:

RETURN

END

