;+
;**********************
PRO GET_COEFF_LEVEL2B,$
;**********************
;INPUT
filename,$             ;retrieval file
;OUTPUT
lin_qua,$              ;lin_qua = 0 -> linear; 1 -> quadratic
angle,$                ;elevation angle
freq,$                 ;frequency array in GHz
z,$                    ;height grid in m
offset,$               ;MLR offset (n_z)
c_l,$                  ;linear retrieval coefficients (n_z x n_f)
c_q,$                  ;quadratic retrieval coefficients (n_z x n_f) ... if lin_qua=1
;KEYWORDS
verbose=verbose
; Abstract:
;* read set of single .RET file for level2b data
;  and return MLR retrieval coefficients
; Author:
; U. Loehnert, S. Crewell
; Date:
; 2011-02-25
; Dependencies:
; -
; Changes:
;-


s=''
str =''
OPENR, unit, filename, /GET_LUN

;..............................................read retrieval type
WHILE s NE '# Retrieval Type' DO BEGIN
   READF,unit,s
   s = STRMID(s,0,16)
ENDWHILE
READF,unit,s
lin_qua = FIX(STRMID(s,3,2))
IF verbose THEN PRINT,'Retrieval type =',lin_qua

;..............................................read frequencies
WHILE s NE '# >>>>>>>>Brightness' DO BEGIN
   READF,unit,s
   s = STRMID(s,0,20)
ENDWHILE
READF,unit,s        ; frequencies
s = STRMID(s,3,STRLEN(s)-3)
words = STRSPLIT(s,' ',/extract)
n_f = N_ELEMENTS(words)
freq   = FLOAT(words)
IF verbose THEN PRINT,'Frequencies=',freq

;..............................................read angles
FOR i = 0,1 DO READF,unit,s
angle = FLOAT(STRMID(s,3,15))
IF verbose THEN PRINT,'Angle=',angle

;..............................................read altitudes
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

;..............................................read offset
offset = FLTARR(n_z)
WHILE str NE '# >>>>>>>>>>>>>>>>>>Offset' DO BEGIN
   READF,unit,s
   str = STRMID(s,0,26)
ENDWHILE
FOR i = 0,n_z-1 DO BEGIN
  READF,unit,s
  ind = STRPOS(s, '#')
  offset(i) = FLOAT(STRMID(s, 3, ind-3))
  IF verbose THEN PRINT,'Offset=',z(i),offset(i)
ENDFOR

;..............................................read linear coefficients
c_l = FLTARR(n_z,n_f)
WHILE str NE 'TL=' DO BEGIN
   READF,unit,s
   str = STRMID(s,0,3)
ENDWHILE
ind = STRPOS(s,'#')
s = STRMID(s,3,ind-3)
c_l(0,*) = STRSPLIT(s,' ',/extract)
IF verbose THEN print,'C_l=',z(0),c_l(0,*),FORMAT='(A4,F9.1,7E14.4)'
FOR i_z = 1,n_z-1 DO BEGIN
  READF,unit,s
  ind = STRPOS(s,'#')
  s = STRMID(s,3,ind-3)
  c_l(i_z,*) = STRSPLIT(s,' ',/extract)
  IF verbose THEN PRINT,'C_l=',z(i_z),c_l(i_z,*),FORMAT='(A4,F9.1,7E14.4)'
ENDFOR

;..............................................read quadratic coefficients
c_q = FLTARR(n_z,n_f)
WHILE str NE 'TQ=' DO BEGIN
   READF,unit,s
   str = STRMID(s,0,3)
ENDWHILE
ind = STRPOS(s,'#')
s = STRMID(s,3,ind-3)
c_q(0,*) = STRSPLIT(s,' ',/extract)
IF verbose THEN print,'C_q=',z(0),c_q(0,*),FORMAT='(A4,F9.1,7E14.4)'
FOR i_z = 1,n_z-1 DO BEGIN
  READF,unit,s
  ind = STRPOS(s,'#')
  s = STRMID(s,3,ind-3)
  c_q(i_z,*) = STRSPLIT(s,' ',/extract)
  IF verbose THEN PRINT,'C_q=',z(i_z),c_q(i_z,*),FORMAT='(A4,F9.1,7E14.4)'
ENDFOR

CLOSE,unit
FREE_LUN,unit


END

