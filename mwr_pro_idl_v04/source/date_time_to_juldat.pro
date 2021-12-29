PRO DATE_TIME_TO_JULDAT, date, time, julian_date

;INPUT SPECFICATIONS
;1.) date must be a character*8 (yyyymmdd)
;2.) time must be decimal in hours (e.g. 23.5 = 23:30)

;OUTPUT
;julian_date

dd = FIX(STRMID(date, 6, 2))
mm = FIX(STRMID(date, 4, 2))
yyyy = FIX(STRMID(date, 0, 4))

fixhour = FIX(time)

min = (time - fixhour)*60.
fixmin = FIX(min)

sec = (min - fixmin)*60.
fixsec = ROUND(sec)

julian_date = JULDAY(mm, dd, yyyy, fixhour, fixmin, fixsec)

;print, date, time
;print, mm, dd, yyyy, fixhour, fixmin, fixsec
;print, julian_date, format = '(f15.5)'

END