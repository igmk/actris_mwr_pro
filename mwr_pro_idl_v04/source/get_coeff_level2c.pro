;+
;**********************
PRO GET_COEFF_LEVEL2C,$
;**********************
;INPUT
filename,$             ;retrieval file
;OUTPUT
angles,$               ;elevation angles
f,$                    ;frequency array in GHz
z,$                    ;height grid in m
offset,$               ;MLR offset (n_z)
coeff_tel,$            ;linear retrieval coefficients (n_z x (n_f + n_bl*n_ang-1))
coeff_sen,$            ;retrieval coefficients for additional sensors (T_gr, q_gr, p_gr, TB_ir) (n_z)
f_bl,$                 ;freqeuncies used for BL scanning
;KEYWORDS
verbose=verbose
; Abstract:
;* read set of single .RET file for level2c retrieval
;  and return MLR retrieval coefficients
; Author:
; U. Loehnert, S. Crewell
; Date:
; 2011-02-25
; Dependencies:
; -
; Changes:
;-

s = ''
str = ''
dummy = -999.

OPENR, unit, filename, /GET_LUN

;***read retrieval type
WHILE s NE '# Retrieval Type' DO BEGIN
   READF, unit, s
   s = STRMID(s, 0, 16)
ENDWHILE
READF, unit, s 
lin_qua = FIX(STRMID(s,3,2))
IF verbose THEN PRINT,'Retrieval type =', lin_qua

;***check if meteorological/infrared sensor data is used in retrieval
;--> if yes, then TS, HS, PS or IS must be set to 1
;--> if no, then TS, HS, PS or IS must be set to 0
TS=0
PS=0
HS=0
IS=0
WHILE s NE '# >>>>>>>>>>>>>>>>>>>Meteorological Sensors' DO BEGIN
 READF, unit, s
 s = STRMID(s, 0, 43)
ENDWHILE

READF, unit, s
READF, unit, s, TS, format = '(a3, i1)'
READF, unit, s
READF, unit, s, HS, format = '(a3, i1)'
READF, unit, s
READF, unit, s, PS, format = '(a3, i1)'
READF, unit, s
READF, unit, s, IS, format = '(a3, i1)'

;***read frequencies
WHILE s NE '# >>>>>>>>Brightness' DO BEGIN
   READF, unit, s
   s = STRMID(s,0,20)
ENDWHILE
READF,unit,s        ; frequencies
s = STRMID(s,3,STRLEN(s)-3)
words = STRSPLIT(s,' ',/extract)
n_f = N_ELEMENTS(words)
f   = FLOAT(words)
IF verbose THEN PRINT,'Frequencies=',f

;***read angles
READF,unit,s
angles = FLTARR(6)
READF,unit,s ; angles
s = STRMID(s, 3, STRLEN(s)-3)
words = STRSPLIT(s,' ',/extract)
n_ang = N_ELEMENTS(words)
angles = FLOAT(words)
IF verbose THEN PRINT,'Angle=',angles

;***read altitudes
WHILE str NE 'AL=' DO BEGIN
   READF,unit,s
   str = STRMID(s,0,3)
ENDWHILE
s = STRMID(s,3,STRLEN(s)-3)
words = STRSPLIT(s,' ',/extract)
n_z = N_ELEMENTS(words)
z   = FLOAT(words)
IF verbose THEN BEGIN
   PRINT,'# Altitudes =',n_z
   ;PRINT,'Altitudes=',z
ENDIF

;***read offset
offset = FLTARR(n_z)
WHILE str NE '# >>>>>>>>>>>>>>>>>>Offset' DO BEGIN
   READF,unit,s
   str = STRMID(s,0,26)
ENDWHILE
FOR i = 0,n_z-1 DO BEGIN
  READF,unit,s
  ind = STRPOS(s, '#')
  s = STRMID(s, 3, ind-3)
  help_me = STRSPLIT(s,' ',/extract)
  offset(i) = FLOAT(help_me)
;  offset(i) = FLOAT(STRMID(s,3,12))
  IF verbose THEN print,'Offset=',z(i),offset(i)
ENDFOR

;***read meteorological/IR sensor coefficients
sum_coeff_sen = TOTAL(TS+HS+PS+IS)
coeff_sen = REPLICATE(dummy, n_z, 4)

IF sum_coeff_sen GT 0 THEN BEGIN 
 WHILE str NE '#     (linear Coefficients!)' DO BEGIN
  READF, unit, s
  str = STRMID(s, 0, 28)
 ENDWHILE 
 FOR i = 0, n_z-1 DO BEGIN
  READF, unit, s
  ind = STRPOS(s, '#')   
  s = STRMID(s, 3, ind-3)
  help_me = STRSPLIT(s,' ',/extract)
  FOR j = 0, sum_coeff_sen-1 DO BEGIN
   coeff_sen(i, j) = help_me(j)
  ENDFOR
  IF verbose THEN PRINT,'coeff_sen=', z(i), coeff_sen(i, *) 
 ENDFOR
ENDIF

;***read linear coefficients
coeff_tel = REPLICATE(dummy, n_z, n_f*n_ang)
help_me = FLTARR(n_ang)

mem = 0
f_bl = dummy

FOR i_f = 0, n_f-1 DO BEGIN 
 WHILE str NE 'TL=' DO BEGIN
  READF, unit, s
  str = STRMID(s, 0, 3)
 ENDWHILE

 ind = STRPOS(s, '#')
 s = STRMID(s, 3, ind-3)
 help_me = STRSPLIT(s,' ',/extract)
 nn = N_ELEMENTS(help_me)
 coeff_tel(0, mem:mem+nn-1) = FLOAT(help_me)

 FOR iz = 1, n_z-1 DO BEGIN
  READF, unit, s
  ind = STRPOS(s, '#')
  s = STRMID(s, 3, ind-3)
  help_me =  STRSPLIT(s,' ',/extract)
  nn = N_ELEMENTS(help_me)
  coeff_tel(iz, mem:mem+nn-1) = FLOAT(help_me) 
 ENDFOR
 IF N_ELEMENTS(help_me) GT 1 THEN BEGIN
  IF f_bl(0) EQ dummy THEN BEGIN
   f_bl = f(i_f) 
  ENDIF ELSE BEGIN
   f_bl = [f_bl, f(i_f)]
  ENDELSE
 ENDIF
 mem = N_ELEMENTS(help_me)+mem
 str = ''
ENDFOR

ii = WHERE(coeff_tel(0, *) NE dummy)
coeff_tel = coeff_tel(*, ii)

IF verbose THEN BEGIN
 print, 'height 0: coeff_sen=', TRANSPOSE(coeff_sen(0, *))
 print, 'height n_z: coeff_sen=', TRANSPOSE(coeff_sen(n_z-1, *))
 print,'height 0: coeff_tel=', TRANSPOSE(coeff_tel(0, *))
 print,'height n_z: coeff_tel=', TRANSPOSE(coeff_tel(n_z-1, *))
ENDIF

CLOSE,unit
FREE_LUN,unit

END

