;+
;************
PRO GET_BLB,$
;************
;INPUT:
filename,$                  ;*.blb file 
;OUPUT:
n,$                         ;number of measurements
time,$                      ;seconds since 1.1.1970
rain,$                      ;rain flag (0/1)
f,$                         ;array of frequencies
n_a,$                       ;number of elevation angles
tb,$                        ;TB arry (f x angles x time)
el,$                        ;array of elevation angles (angles)
az,$                        ;array azimuth angle (angles)
;KEYWORDS
verbose=verbose,$
angle_calc = angle_calc     ;default is 0, 1 means calculate angles in old manner (see below)
; Abstract:
; * reads binary *.blb (boundary layer elevation scanning) files
; Authors:
; S. Crewell, U. Loehnert
; Date:
; 2008-09-12
; Dependencies:
; -
; Changes:
; 2014-10-24 (UL)
; - bux fix to to read azimuth and elevation in a consistent manner, not sure when
;   RPG changed there angle coding, so option is built in to switch back to the old
;   type of angle calculations.
; changed program ...
;-

ON_IOERROR, BLB_END_1
OPENR,unit,filename,/GET_LUN
code=1l
READU,unit,code
IF verbose THEN print, 'file code: ', code

IF code LT 567845847 OR code GT 567845848 THEN BEGIN
 PRINT,'Error in BLB file code'
 print, 'code should be: 567845847 or 567845848'
 print, 'code is: ', code
 GOTO, BLB_END
ENDIF

IF N_ELEMENTS(angle_calc) EQ 0 THEN angle_calc = 0

n=0l
READU,unit,n
IF verbose THEN print,'Number of samples=',n
IF n LT 1 THEN GOTO,BLB_END

nf=14l
IF code EQ 567845848 THEN READU,unit,nf
IF verbose THEN print,'Number of frequencies=',nf

xmin = FLTARR(nf)
xmax = FLTARR(nf)
READU,unit,xmin,xmax
IF verbose THEN BEGIN
  print,'Min=',xmin
  print,'Max=',xmax
ENDIF

time_ref = 0l
READU,unit,time_ref
IF verbose THEN BEGIN
  CASE time_ref of
  1: print,'UTC'
  0: print,'Local Time'
  ELSE: print,'No valid time'
  ENDCASE
ENDIF

IF code EQ 567845847 THEN READU,unit,nf

f = FLTARR(nf)
READU,unit,f
IF verbose THEN PRINT,'Frequencies=',f

n_a = 0l
READU,unit,n_a
IF verbose THEN PRINT,'Number of angles=', n_a

ang = FLTARR(n_a)
READU, unit, ang
IF verbose THEN PRINT,'Angles=', ang

time = LONARR(n)
rain = BYTARR(n)
tb = FLTARR(nf,n_a,n)

x1 = 0l
x2 = 0b
x = -2.
help = FLTARR(n_a)
FOR i = 0, n-1 DO BEGIN
  READU, unit, x1, x2
  time(i) = x1
  rain(i) = x2
  FOR i_f = 0, nf-1 DO BEGIN 
   READU, unit, help
   tb(i_f,*,i) = help
   IF verbose THEN print, i, time(i), rain(i), i_f, tb(i_f,0,i)
   READU, unit, x
  ENDFOR
ENDFOR

el = FLTARR(n_a)
az = FLTARR(n_a)
sign = 1
FOR iang = 0, n_a-1 DO BEGIN

;*the following lines (from ## to ##) were created on 24.10.14 to correct the angle calculation from RPG .BRT files
;##
 x = ang(iang)
 IF angle_calc EQ 0 THEN BEGIN
  els = x - 100.*FIX(x/100.)
  azs = (x-els)/1000.

  IF azs LE 360. THEN BEGIN ; first quadrant
   el(iang) = els
   az(iang) = azs
  ENDIF ELSE IF azs GT 360. AND azs LT 1000. THEN BEGIN
   print, 'Inconsistency in angle calculation!'
   print, 'ABORT in get_brt.pro'
   stop
  ENDIF ELSE IF azs GT 1000. THEN BEGIN
   az(iang) = azs - 1000.
   el(iang) = 100. + els
  ENDIF
 ENDIF
;##

;*the following lines (from ### to ###) state the angle calculation until 24.10.14
;###
 IF angle_calc EQ 0 THEN BEGIN 
  IF el(iang) LT 0. THEN sign = -1.
  az(iang) = sign*FLOOR(el(iang)/100.)/10.
  el(iang) = el(iang) - (sign*az(iang)*1000.)
 ENDIF
;###
ENDFOR

BLB_END:

FREE_LUN, unit

BLB_END_1:

END

