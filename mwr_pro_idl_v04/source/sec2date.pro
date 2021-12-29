;+
;************************************************************************************************************************************
PRO SEC2DATE, $
;************************************************************************************************************************************
; input:
time, $
; output:
date, hours
;-
; transforms seconds from 01/01/2001 to date vector and hours 

jj = where(time eq time) 
date  = replicate(0., 6, n_elements(time))
hours = replicate(0., n_elements(time))

jd = julday(1, 1, 2001, 0, 0, time(jj))
caldat, jd, m, d, y, h, mn, s

date(0, jj) = y
date(1, jj) = m 
date(2, jj) = d
date(3, jj) = h
date(4, jj) = mn
date(5, jj) = s

s = round(s)  
hours(jj) = h(jj) + mn(jj)/60. + s(jj)/3600. 

END