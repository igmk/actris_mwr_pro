;+
;************
PRO GET_IRT,$
;************
;INPUT:
filename,$                  ;*.blb file 
;OUPUT:
time,$                      ;seconds since 1.1.2001
wavel,$                     ;center wavelength (in um)
tir,$                       ;TB (IR)
el,$                        ;array of elevation angles
az,$                        ;array azimuth angle
;KEYWORDS
verbose=verbose
; Abstract:
; * reads binary *.irt (RPG infrared thermometer) files
; Authors:
; U. Loehnert
; Date:
; 2012-02-10
; Dependencies:
; -
; Changes:
; XXXX-XX-XX:
; changed program ...
;-

ON_IOERROR, SPEC_END
ON_ERROR, 2

OPENR,unit,filename,/GET_LUN,ERROR=err

IF err NE 0 THEN BEGIN
  n = 0
  GOTO,END_IRT
ENDIF
ON_ERROR,0

code=1l
READU,unit,code
IF code GT 671112496 OR code LT 671112495 THEN PRINT,'Error in SPEC file code'

IF code EQ 671112496 THEN ir_type = 'kt2'
IF code EQ 671112495 THEN ir_type = 'kt1'

IF N_ELEMENTS(verbose) EQ 0 THEN verbose = 0

n=0l
READU, unit, n
IF verbose THEN print,'Number of samples=',n
IF n LT 1 THEN GOTO,SPEC_END

xmin = 1e1 ; predefine as single precision (4byte float)
xmax = 1e1 ; predefine as single precision (4byte float)
READU,unit,xmin,xmax
IF verbose THEN print,xmin,xmax

time_ref = 0l ; predefine as long int = 4byte integer = 32bit int...
READU,unit,time_ref
IF verbose THEN BEGIN
CASE time_ref of
  1: print,'UTC'
  0: print,'Local Time'
  ELSE: print,'No valid time'
  ENDCASE
ENDIF

n_wavel = 1 ; default
wavel = 11.1
IF ir_type EQ 'kt2' THEN BEGIN
 n_wavel = 0l
 READU, unit, n_wavel
 wavel = FLTARR(n_wavel)
 READU, unit, wavel
ENDIF

time = LONARR(n)
rain = BYTARR(n)
tir  = FLTARR(n_wavel, n)
az  = REPLICATE(-999., n)
el  = REPLICATE(-999., n)

T = 0l ; predefine as long int (4byte/32bit)
RF = 0b ; predefine as byte (1byte/8bit)
IF ir_type EQ 'kt2' THEN IRT = REPLICATE(-999., n_wavel) ; predefine as single (4byte float)
IF ir_type EQ 'kt1' THEN IRT = -999.; predefine as single (4byte float)

ang_irt = -999.
FOR i = 0l,LONG(n-1l) DO BEGIN
  IF ir_type EQ 'kt2' THEN READU,unit,T,RF,IRT,ang_irt
  IF ir_type EQ 'kt1' THEN READU,unit,T,RF,IRT
  time(i) = T
  rain(i) = RF
  tir(0:n_wavel-1, i) = IRT
  sign = 1
  IF ang_irt LT 0. THEN sign = -1.
  IF ir_type EQ 'kt2' THEN BEGIN
   az(i)   = sign*FLOOR(ang_irt/100.)/10.
   el(i)   = ang_irt - (sign*az(i)*1000.)
  ENDIF
 
  IF verbose THEN print,i,time(i),rain(i),tir(*, i),az(i),el(i)
ENDFOR


SPEC_END:

CLOSE,unit
FREE_LUN,unit

END_IRT:

END

