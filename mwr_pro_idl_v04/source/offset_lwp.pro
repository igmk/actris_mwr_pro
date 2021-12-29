;+
;***************
PRO OFFSET_LWP,$
;***************
;input
mwr_dir,$                       ; directory
date,$                          ; yymmdd to be plotted
time,$                          ; time 
ele,$                           ; elevation angle array
lwp,$                           ; liquid water path
lwp_std_thres,$                 ; threshold to be exceed within 2-min intervals to decide whether clear or cloudy
flag_2a,$                       ; 2a threshold - offset correction only carried out if flag_2a is zero
;output
lwp_cor,$                       ; offset corrected lwp
;keywords
verbose=verbose
; $Id:$
; Abstract: 
; * automatic zenith lwp offset correction - every
;   processed day is analysed for clear sky situations. The criterion "clear sky" is set to yes/no
;   for every 20 minute interval based solely on the two minute zenith lwp-standard-deviation [in g/m2]
;   of the original lwp values. If, within the 20-minute-interval, each two-minute-interval shows a
;   lwp-standard-deviation of lower than lwp_std_thres, then the 20-minute interval is considered as
;   "clear-sky". Note, lwp_std_thres is instrument, retrieval and climate-zone dependent and thus must
;   be chosen carefully! If set_lwp_off is set, a file lwp_offset_yyyy.sav is modified containing the
;   dates and lwp offsets during clear sky intervals. Off-zenith measurements are NOT corrected 

; Authors:
; U. Loehnert
; Date:
; 2012-02
; Dependencies:
; - 
; Changes
; 20130910 (UL) : removed faulty off-zenith lwp correction; off-zenith lwp are NOT offset corrected
;-


IF verbose THEN print, 'OFFSET_LWP: calculating lwp offset correction'

;***check for clear-sky times on this day

;*make 20-minute time array for this day
n_twe = 73
twe_min = FINDGEN(n_twe)/3.

;*make 60-minute time array for this day
;n_twe = 25
;twe_min = FINDGEN(n_twe)

two_min = 2./60.
n_off = 0
time_off = -999.
lwp_off = -999.

;**loop over 20-minute intervals and find out if they are clear or cloudy

FOR i = 0, n_twe-2 DO BEGIN
;print, i, twe_min(i)
;*loop over 2-minute intervals within 20-minute interval and calculate lwp std
 n_std = 0
 clear = 0
 cloudy = 0
 lwp_off_x = -999.
 FOR j = 0, 9 DO BEGIN
  i_two = WHERE(time GE twe_min(i)+FLOAT(j)*two_min AND time LT twe_min(i)+FLOAT(j+1)*two_min AND ele GT 89.4 AND ele LT 90.6 AND lwp NE -999. AND flag_2a EQ 0., nn)
;print, twe_min(i), twe_min(i)+FLOAT(j)*two_min, twe_min(i)+FLOAT(j+1)*two_min
  IF nn GE 30. THEN BEGIN
   lwp_std = STDDEV(lwp(i_two))
;print, '**',  lwp_std
   IF lwp_std GE lwp_std_thres/1000. THEN BEGIN
    cloudy = 1 
   ENDIF ELSE BEGIN
    clear = 1
    IF lwp_off_x(0) EQ -999. THEN BEGIN
     lwp_off_x = lwp(i_two)
    ENDIF ELSE BEGIN
     lwp_off_x = [lwp_off_x, lwp(i_two)]
    ENDELSE 
   ENDELSE
   n_std = n_std + 1
  ENDIF
 ENDFOR ; loop over 2-minute intervals

 IF verbose THEN print, 'OFFSET_LWP: CLEAR/CLOUDY ', clear, cloudy

 IF clear EQ 1 AND cloudy EQ 0 THEN BEGIN
  IF lwp_off(0) EQ -999. THEN BEGIN
   lwp_off = MEAN(lwp_off_x)
   DATE_TIME_TO_JULDAT, '20'+date, 0.5*(twe_min(i+1)+twe_min(i)), jd
   time_off = jd
  ENDIF ELSE BEGIN
   lwp_off = [lwp_off, MEAN(lwp_off_x)]
   DATE_TIME_TO_JULDAT, '20'+date, 0.5*(twe_min(i+1)+twe_min(i)), jd
   time_off = [time_off, jd]
  ENDELSE
 ENDIF
ENDFOR ; loop over 20-minute intervals

time_off_c = time_off
lwp_off_c = lwp_off
n_off = N_ELEMENTS(time_off)

;**add time_off and lwp_off to lwp_off_yyyy.sav
yyyy = '20'+STRMID(date, 0, 2)
check_exist = FILE_TEST(mwr_dir+'/lwp_off_'+yyyy+'.sav')

IF check_exist EQ 0 AND time_off(0) NE -999. THEN BEGIN

 SAVE, filename=mwr_dir+'/lwp_off_'+yyyy+'.sav', time_off, lwp_off

ENDIF

IF check_exist EQ 1 THEN BEGIN
 RESTORE, filename=mwr_dir+'/lwp_off_'+yyyy+'.sav'
 CALDAT, time_off, mon_c, day_c, yea_c, hou_c, min_c, sec_c
 yea_d = LONG(STRMID(date, 0, 2)) + 2000l
 mon_d = LONG(STRMID(date, 2, 2))
 day_d = LONG(STRMID(date, 4, 2))
 
;*exclude existing lwp offsets of this day from offset array - these will replaced by the new calculations

 result = WHERE(yea_c EQ yea_d AND mon_c EQ mon_d AND day_c EQ day_d, COMPLEMENT=i_c)
 IF i_c(0) NE -1 THEN BEGIN
  time_off = time_off(i_c)
  lwp_off = lwp_off(i_c)  
;*add new lwp_off_c to lwp_off
  time_off = [time_off, time_off_c]
  lwp_off = [lwp_off, lwp_off_c]
  isort = SORT(time_off)
  time_off = time_off(isort)
  lwp_off = lwp_off(isort)
 ENDIF 

;*save new offsets to lwp_off_yyyy.sav

 iii = WHERE(time_off GT 0.)
 IF iii(0) NE -1 THEN BEGIN
  time_off = time_off(iii)
  lwp_off = lwp_off(iii)
 ENDIF 

 SAVE, filename=mwr_dir+'/lwp_off_'+yyyy+'.sav', time_off, lwp_off
ENDIF

;** in case of month = January, also check data from December
IF STRMID(date, 2, 2) EQ '01' THEN BEGIN

;*check if last year's offset file exists
 ly = FIX(STRMID(date, 0, 2))-1
 IF ly LT 10 THEN yy_ly = '0'+STRING(ly, format='(i1)')
 IF ly GE 10 THEN yy_ly = STRING(ly, format='(i2)')
 yyyy_ly = '20'+yy_ly
 check_exist = FILE_TEST(mwr_dir+'/lwp_off_'+yyyy_ly+'.sav')

 IF check_exist EQ 1 THEN BEGIN
  time_off_ty = time_off
  lwp_off_ty = lwp_off
  RESTORE, filename=mwr_dir+'/lwp_off_'+yyyy_ly+'.sav'
  time_off = [time_off, time_off_ty]
  lwp_off = [lwp_off, lwp_off_ty]
 ENDIF

ENDIF; if month = January


DATE_TIME_TO_JULDAT, '20'+date, time, jd_lwp

;***loop over lwp values
n_lwp = N_ELEMENTS(lwp)

;*search for offsets around time
FOR i = 0l, n_lwp-1l DO BEGIN
 ii_high = WHERE(time_off GE jd_lwp(i), n_high)
 ii_low = WHERE(time_off LT jd_lwp(i), n_low)

;*interpolate offset if offsets are available before and after time
 IF ii_high(0) NE -1 AND ii_low(0) NE -1 THEN BEGIN
  IF lwp_off(ii_high(0)) NE -999. AND lwp_off(ii_low(n_low-1)) NE -999. THEN BEGIN
   m = (lwp_off(ii_high(0))-lwp_off(ii_low(n_low-1)))/(time_off(ii_high(0))-time_off(ii_low(n_low-1)))
   offset = lwp_off(ii_low(n_low-1)) + m*(jd_lwp(i)-time_off(ii_low(n_low-1)))
;*off-zenith lwp are NOT corrected
   IF ele(i) GT 89.4 AND ele(i) LT 90.6 THEN lwp_cor(i) = lwp(i) - offset
  ENDIF

;*use nearest neighbor if no values available after time
 ENDIF ELSE BEGIN
  CLOSEST, jd_lwp(i), time_off, dmin, ind
  offset = lwp_off(ind)
;*off-zenith lwp are NOT corrected
   IF ele(i) GT 89.4 AND ele(i) LT 90.6 THEN lwp_cor(i) = lwp(i) - offset
 ENDELSE

ENDFOR ; loop over all lwp


RETURN

END
