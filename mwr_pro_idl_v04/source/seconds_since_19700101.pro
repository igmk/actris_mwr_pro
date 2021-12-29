PRO SECONDS_SINCE_19700101, date, time, seconds_return

;seconds_return: returns seconds since 01.01.1970, 00:00 UTC
;date: STRING in format 'yyyymmdd'
;time: array of decimal hours since 0 UTC

date_start = JULDAY(1,1,1970,0,0,0)

DATE_TIME_TO_JULDAT, date, time, julian_date

date_s = julian_date - date_start

seconds_return = date_s*86400d

END
