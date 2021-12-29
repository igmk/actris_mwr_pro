;+
;**********************
PRO GET_COEFF_LEVEL2A,$
;**********************
;INPUT:
verbose,$
filename,$
;OUTPUT:
lin_qua,$
angle,$               ;elevation angle at which retrtieval is valid
f,$                   ;frequenc(y)ies at which retrieval is valid
offset,$              ;offset MLR
c_l,$                 ;linear coefficients (X frequencies)
c_q                   ;quadratic coefficients (X frequencies)
; Abstract:
; * read retrieval coefficients for IWV, LWP or ZWD
; Author:
; S. Crewell, U. Loehnert
; Date:
; 2006-02-14
; Dependencies:
; -
; Changes:
; XXXX-XX-XX: ???
;-

s=''
str = ''
OPENR,unit,filename,/GET_LUN, ERROR = err

IF (err EQ 0) THEN BEGIN
  ; read file
  ;..............................................read retrieval type
  WHILE str NE '# Retrieval Type' DO BEGIN
    READF,unit,s
    str = STRMID(s,0,16)
  ENDWHILE
  READF,unit,s        ; frequencies
  lin_qua = FIX(STRMID(s,3,2))
  IF verbose THEN PRINT,'Retrieval type =',lin_qua
  
  ;..............................................read frequencies
  WHILE s NE '# >>>>>>>>Brightness' DO BEGIN
    s=''
    READF,unit,s
    s = STRMID(s,0,20)
  ENDWHILE
  READF, unit, s, format = '(a200)'        ; frequencies
  s = STRMID(s,3,STRLEN(s)-3)
  words = STRSPLIT(s,' ',/extract)
  n_f = N_ELEMENTS(words)
  f   = FLOAT(words)
  IF verbose THEN PRINT,'Frequencies: ',f

  ;..............................................read angles
  FOR i = 0,1 DO READF,unit,s
  angle = FLOAT(STRMID(s,3,11))
  IF verbose THEN PRINT,'Angle=',angle
  
  ;..............................................read offset
  FOR i = 0,5 DO READF,unit,s
  offset = FLOAT(STRMID(s,3,25))
  IF verbose THEN PRINT,'Offset=',offset
  
  ;..............................................read linear coefficients

  FOR i = 0,7 DO READF,unit,s
  s = STRMID(s,3,STRLEN(s)-3)
  words = STRSPLIT(s,' ',/extract)
  n_cl = N_ELEMENTS(words)
  c_l   = FLOAT(words)
  IF verbose THEN PRINT,'C_l=',c_l
  
  ;..............................................read quadratic coefficients
  FOR i = 0,1 DO READF,unit,s
  s = STRMID(s,3,STRLEN(s)-3)
  words = STRSPLIT(s,' ',/extract)
  n_cq = N_ELEMENTS(words)
  c_q   = FLOAT(words)
  IF verbose THEN PRINT,'C_q=',c_q
  
 FREE_LUN,unit

ENDIF ELSE BEGIN ; (err NE 0)
  PRINTF, -2, 'Error while opening "',filename,'"'
  PRINTF, -2, !ERROR_STATE.MSG
ENDELSE



END

