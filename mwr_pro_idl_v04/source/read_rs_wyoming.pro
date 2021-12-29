;+
;********************
PRO READ_RS_WYOMING,$
;********************
;INPUT:
filename,$                      ;name of file to be read
;OUTPUT:
time,$                          ;time of ascent (Julian)
iwv                             ;Integrated Water Vapor
; $Id: read_rs_wyoming.pro,v 1.1 2008/10/20 14:00:05 sredl Exp $
; Abstract:
; * return IWV of Wyoming sondes
; Authors:
; U. Loehnert
; Date:
; 2008-09-25
; Dependencies:
; -
; Changes:
; XXXX-XX-XX:
; changed program ...
;-

OPENR, unit, filename, /GET_LUN

s = ''
time = -999.
iwv = -999.

WHILE NOT EOF(unit) DO BEGIN

 WHILE STRMID(s, 27, 17) NE 'Observation time:'  DO BEGIN
  IF EOF(unit) THEN GOTO, DONE
  READF, unit, s
 ENDWHILE 

 year = 2000 + FIX(STRMID(s, 45, 2))
 month = FIX(STRMID(s, 47, 2))
 day = FIX(STRMID(s, 49, 2))
 hour = FIX(STRMID(s, 52, 2))

 time = JULDAY(month, day, year, hour, 0, 0) 

 WHILE STRMID(s, 0, 12)  NE 'Precipitable'  DO BEGIN
   IF EOF(unit) THEN GOTO, DONE
   READF, unit, s
 ENDWHILE 

 pos = STRPOS(s,'.')
 xx = STRMID(s, pos-2, 2)
 yy = STRMID(s, pos+1, 2)
 iwv = FLOAT(xx) + FLOAT(yy)/100.
 GOTO, DONE

ENDWHILE

DONE:

FREE_LUN, unit

RETURN

END
