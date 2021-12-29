;************
PRO GET_DEC,$
;************
verbose,$
;INPUT:
infile,$
;OUPUT:
time_re,$      ;julian date
tb_re,$        ;TB(3channelsXn_measurements) in K 
el_re,$        ;elevation angle in deg 
freq           ;frequencies in GHz

; $Id: $
; Abstract:
; * read raw data from RESCOM MWR data files: .dec format
; Authors:
; U. Loehnert/S. Crewell
; Date:
; 2011-02-01
; Dependencies:
; -
; Changes:
; -
;-

IF verbose THEN print, 'Entering get_dec.pro'

s = ''
yyyy = ''
mm = ''
dd = ''

OPENR, unit, infile, /GET_LUN
FOR i = 0, 1 DO READF, unit, s

READF, unit, s, f1, s, f2, s, f3, format = '(a7, f7.3, a2, f7.3, a2, f7.3)'
READF, unit, s, dd, s, mm, s, yyyy, format = '(a7, a2, a1, a2, a1, a4)'
READF, unit, s, hh, s, mi, s, ss, format = '(a7, f2.0, a1, f2.0, a1, f2.0)'

freq = [f1, f2, f3]

n_samp = 86400l
TB_re = FLTARR(3, n_samp)
time_re = FLTARR(n_samp)
el_re = FLTARR(n_samp)
dt = FLTARR(n_samp-1)
n_re = 0l

print, infile
WHILE NOT EOF(unit) DO BEGIN
 READF, unit, s
 new = STRSPLIT(s, /Extract)
 IF N_ELEMENTS(new) GT 8 THEN BEGIN
  FOR i = 0, 2 DO TB_re(i, n_re) = new(i+2)
  time_re(n_re) = new(0)
  el_re(n_re) = new(1)
  if n_re GT 0 THEN dt(n_re-1) = time_re(n_re)-time_re(n_re-1)
  n_re = n_re + 1l
 ENDIF
ENDWHILE

FREE_LUN, unit


time_re = time_re(0:n_re-1)
el_re = el_re(0:n_re-1) 
tb_re = REFORM(tb_re(*,0:n_re-1))

start_time = hh + mi/60d + ss/3600d
DATE_TIME_TO_JULDAT, yyyy+mm+dd, (time_re+start_time)/3600d, jd
time_re = jd

RETURN

END
