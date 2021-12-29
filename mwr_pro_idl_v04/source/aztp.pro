;+
;*********
PRO AZTP,$
;*********
;keywords
par=par,$                       ; name of parameter to be plotted 
filename=filename,$             ; ps filename to be created
Date=date,$                     ; YYMMDD to be plotted
mwr_meas=mwr_meas,$             ; three letter code
time=time,$                     ; time 
data=data,$                     ; parameter values
ele=ele,$                       ; elevation angle array
azi=azi,$                       ; azimuth angle array
ytit=ytit,$                     ; ytitle bar
el_aztp=el_aztp,$               ; elevation values for AZTP (array for maximum two)
az_scan_thres=az_scan_thres,$   ; delta az for definen start/end of azimuth scan   
color_tab=color_tab,$           ; color table to load (5 for IWV, 1 for LWP optimal)
flag_2a=flag_2a,$               ; 
verbose=verbose           
; $Id:$
; Abstract: 
; * plots IWV or LWP azimuth-time contours (AZTPs)
; Authors:
; U. Loehnert
; Date:
; 2012-01
; Dependencies:
; - 
; Changes
; 2013-10-13 (UL): changed valid plotting range to +-0.6 deg from +-0.5deg 
; 2014-10-27 (UL): plotting range adjusted so that solar flagged values are NOT considered
;-

!P.charsize = 0.8
set_plot, 'ps'
!p.multi = [0, 2, 1]
LOADCT, color_tab
IF N_ELEMENTS(filename) EQ 0 OR filename EQ '' THEN BEGIN
 print, 'AZTP: no filename specified, aborting!'
 GOTO, SKIP_AZTP
ENDIF 

device, /color, file = filename, BITS_PER_PIXEL = 8, /landscape
!P.font = 3

n_el = N_ELEMENTS(el_aztp)
IF n_el EQ 0 THEN BEGIN
 print, 'AZTP: no elevation angles specified, aborting!'
 GOTO, SKIP_AZTP
ENDIF

plot_off_xx = 0.5
n_azi = N_ELEMENTS(azi)

;**identify start and end times of azimuth scans
in_scan = 0
scan_start = -999.
scan_stop = -999.
FOR j = 0l, n_azi-4l DO BEGIN

;**check for scan start times
; da = ABS(azi(j)-azi(j+1))
; IF da GT az_scan_thres AND in_scan EQ 0 THEN BEGIN

 IF ABS(azi(j)-azi(j+1)) NE 0 AND ABS(azi(j+1)-azi(j+2)) NE 0 AND ABS(azi(j+2)-azi(j+3)) NE 0 AND in_scan EQ 0 THEN BEGIN  
  IF scan_start(0) EQ -999. THEN BEGIN
   scan_start = time(j) 
  ENDIF ELSE BEGIN
   scan_start = [scan_start, time(j)]    
  ENDELSE    
  in_scan = 1
 ENDIF 

;**check for scan end times
; da1 = ABS(azi(j+1)-azi(j))
; da2 = ABS(azi(j+2)-azi(j))   
; IF da2 LE az_scan_thres AND da1 LE az_scan_thres AND in_scan EQ 1 THEN BEGIN
 IF in_scan EQ 1 AND ABS(azi(j)-azi(j+1)) LE az_scan_thres AND ABS(azi(j)-azi(j+2)) LE az_scan_thres THEN BEGIN 
  IF scan_stop(0) EQ -999. THEN BEGIN
   scan_stop = time(j)
  ENDIF ELSE BEGIN
   scan_stop = [scan_stop, time(j)]
  ENDELSE
  in_scan = 0   
 ENDIF
ENDFOR

;**validate if number of scan start and scan end times are the same 
n_scan = N_ELEMENTS(scan_start)
n_scan_stop = N_ELEMENTS(scan_stop)
IF n_scan NE n_scan_stop THEN BEGIN
;*check if scan end was not detected due to end of measurement
 IF n_scan-1 EQ n_scan_stop THEN BEGIN
  IF scan_start(n_scan-2) LT scan_stop(n_scan-2) AND scan_start(n_scan-1) GT scan_stop(n_scan-2) THEN BEGIN
   scan_stop = [scan_stop, time(n_azi-1)]
  ENDIF
 ENDIF ELSE BEGIN
  print, 'AZTP: error in scan number counting!'
 ENDELSE
ENDIF

;**determine max/min of all values
r_max = -999.
r_min = 999.
par_plot = -999.
el_plot = -999.

FOR i_el = 0, n_el-1 DO BEGIN
 i_plot = WHERE(ele GT el_aztp(i_el)-0.6 AND ele LT el_aztp(i_el)+0.6 AND ((flag_2a AND 64) EQ 0))

 IF i_plot(0) NE -1 THEN BEGIN
;*perform air-mass correction to zenith
  el_plot = ele(i_plot)
  par_plot = data(i_plot)*SIN(!dpi*el_plot/180.)
 ENDIF

 r_max = MAX([MAX(par_plot), r_max])  
 r_min = MIN([MIN(par_plot), r_min])
 IF r_min LT 0. THEN r_min = 0.
 IF par EQ 'IWV' AND r_max GT 40. THEN r_max = 40.

ENDFOR

;**create array of delta between two scans
IF n_scan GT 1 THEN BEGIN
 dt = REPLICATE(-999., n_scan-1)
 FOR i = 0, n_scan-2 DO dt(i) = scan_start(i+1) - scan_start(i) 

 FOR i_el = 0, n_el-1 DO BEGIN

  plot_off_x = FLOAT(i_el)*plot_off_xx
  el_c = STRING(el_aztp(i_el), format = '(i2)')

  tit = par+' (az_vs_time)@'+el_c+' deg, ' + mwr_meas +', ' + date

  PLOT, [0, 360], [0, 24], /NODATA, Title=tit,$
  xrange = [0, 360], xstyle=1, xcharsize = 1.2, ycharsize = 1.2,$
  xticklen=0.02, xticks=4, xtickname = ['N', 'E', 'S', 'W', 'N'],$
  yrange = [0, 24], ystyle=1, ytitle = 'Time of day [UTC]',$
  yticklen = -0.02, yticks=8, position=[0.1+plot_off_x,0.15,0.45+plot_off_x,0.9]

  i_plot = WHERE(ele GT el_aztp(i_el)-0.6 AND ele LT el_aztp(i_el)+0.6)
  IF i_plot(0) EQ -1 THEN BEGIN
   print, 'AZTP: no elevation angles found @', el_aztp(i_el)
   GOTO, LOOP_END
  ENDIF

  IF i_plot(0) NE -1 THEN BEGIN
   time_plot = time(i_plot)
   flag_plot = flag_2a(i_plot)
   az_plot = azi(i_plot)
   el_plot = ele(i_plot)
;*perform air-mass correction to zenith
   par_plot = data(i_plot)*SIN(!dpi*el_plot/180.)
  ENDIF

;**plot scans in countour
  FOR j = 0, n_scan-2 DO BEGIN
   tim = (scan_start(j)+scan_stop(j))/2.  
   jj = WHERE(time_plot GE scan_start(j) AND time_plot LE scan_stop(j), n_jj)
   FOR k = 0, n_jj-2 DO BEGIN   
    da = ABS(az_plot(jj(k+1)) - az_plot(jj(k)))
    IF da GT 300. THEN da = 360.-da
    x = [az_plot(jj(k)), az_plot(jj(k)), az_plot(jj(k))+da, az_plot(jj(k))+da]
    y = [tim, tim+dt(j), tim+dt(j), tim]
    col = par_plot(jj(k))

;*extra check
;    IF ABS(col-MEAN(par_plot(jj))) GT 2. AND STDDEV(iwv_plot(jj)) GT 1.0 THEN col = r_max
    IF col GT r_max THEN col = r_max
    IF col LT r_min THEN col = r_min

;*flag solar interference black
    IF (flag_plot(jj(k)) AND 64) EQ 64 THEN col = r_max
    POLYFILL, x, y, color=(r_max-col)*255./(r_max-r_min) 

   ENDFOR
  ENDFOR; loop over n_scan 
 LOOP_END:
 ENDFOR ; loop over n_el
ENDIF

bar = 255-INDGEN(255)
n_ticks = 5

FOR i=0,n_ticks DO BEGIN
 z = r_min + (r_max-r_min)*i/n_ticks
 x = 0.2 + i*0.6/n_ticks
 XYOUTS, x, 0.08, '|', alignment=0.5, size=0.8, /NORMAL
 XYOUTS, x, 0.05, STRING(z, FORMAT='(F6.1)'), charsize=0.8, alignment=0.5, /NORMAL
ENDFOR

XYOUTS, 0.5,  0.015, ytit + ' (airmass-corrected)', alignment=0.5, charsize=1.2, /NORMAL
XYOUTS, 0.5, -0.005, 'solar disturbance blocked black', alignment=0.5, charsize=0.7, /NORMAL
TV, bar, 0.2, 0.07, xsize=0.6, ysize=0.02, /NORMAL

DEVICE, /CLOSE

SKIP_AZTP:
END
