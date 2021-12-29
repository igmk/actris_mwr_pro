PRO GET_PREV_DATE, date, year_prev, month_prev, day_prev

;input : date string, yyyymmdd
;ouput : string yyyy, mm, dd (prev day)

DATE_TIME_TO_JULDAT, date, 0.0, jd
jd = jd - 1l
CALDAT, jd, month_prev, day_prev, year_prev
year_prev = STRING(year_prev, format = '(i4)')
IF month_prev LT 10 THEN month_prev = '0'+STRING(month_prev, format = '(i1)')
IF month_prev GE 10 THEN month_prev = STRING(month_prev, format = '(i2)')
IF day_prev LT 10 THEN day_prev = '0'+STRING(day_prev, format = '(i1)')
IF day_prev GE 10 THEN day_prev = STRING(day_prev, format = '(i2)')

RETURN

END