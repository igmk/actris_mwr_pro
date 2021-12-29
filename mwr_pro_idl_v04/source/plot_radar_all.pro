;############################
PRO PLOT_RADAR_ALL,$
;this routine reads all radar data for a specific date
;checks also manual filter file and exclude data if needed
;plots HATPRO LWP and radar Doppler moments, if data are available
;****INPUT
path_cr, $
filter_file_radar, $
filename1_plot_radar,$
filename2_plot_radar,$
date,$ ;'YYMMDD'
station, $; 'jue'
time_w, $;
lwp,$
lwp_min,$
lwp_max
;############################

;path_cr = '/data/data_hatpro/jue/data/mmclx/incoming/'
;filter_file_radar = '/data/data_hatpro/jue/data/mmclx/radar_filter.dat'
date_m= STRMID(date, 0, 4)
date_d = STRMID(date, 4, 2)

;****read radar filter file to check if measurements need to be filtered out manually
a1 = ''
n_filter = 0
i_filter = 0
date_filter = REPLICATE('', 1000)
time_filter_b = REPLICATE(-999., 1000, 10)
time_filter_e = REPLICATE(-999., 1000, 10)

OPENR, unit_filter, filter_file_radar, /GET_LUN
FOR i = 0, 13 DO READF, unit_filter, a1
WHILE NOT EOF(unit_filter) DO BEGIN
 READF, unit_filter, a1, n_filter, x1, x2, format = '(a6, 1x, i3, 2f6.2)'
 date_filter(i_filter) = a1
 time_filter_b(i_filter, 0) = x1
 time_filter_e(i_filter, 0) = x2  
 IF n_filter GT 1 THEN BEGIN
  FOR i = 1, n_filter-1 DO BEGIN
   READF, unit_filter, x1, x2, format = '(10x, 2f6.2)'
   time_filter_b(i_filter, i) = x1
   time_filter_e(i_filter, i) = x2
  ENDFOR
 ENDIF
 i_filter = i_filter + 1
ENDWHILE
CLOSE,unit_filter
date_filter=date_filter(0:i_filter-1)
time_filter_b=time_filter_b(0:i_filter-1,*)
time_filter_e=time_filter_e(0:i_filter-1,*)

;****unzip, read and zip radar data: the reading procedure takes into account that several files for one day may exist
SPAWN, 'date +%y%m%d', date_sys
IF date_m+date_d NE date_sys THEN BEGIN
print, path_cr+'/20'+STRMID(date_m, 0, 2)+'/'+date_m+'/20'+date_m+date_d+'*.mmclx*gz'
res = FILE_SEARCH(path_cr+'/20'+STRMID(date_m, 0, 2)+'/'+date_m+'/20'+date_m+date_d+'*.mmclx*gz', count=count_rf)
ENDIF

IF date_m+date_d EQ date_sys THEN BEGIN
print, path_cr+'/20'+STRMID(date_m, 0, 2)+'/'+date_m+'/20'+date_m+date_d+'*.mmclx.??'
res = FILE_SEARCH(path_cr+'/20'+STRMID(date_m, 0, 2)+'/'+date_m+'/20'+date_m+date_d+'*.mmclx.??', count=count_rf)
ENDIF

flag_first_file=0

IF res(0) NE '' THEN BEGIN

 FOR icr=0, count_rf-1 DO BEGIN


IF date_m+date_d NE date_sys THEN BEGIN
  print, 'unzipping radar file: ', res(icr)
  SPAWN, 'gunzip -f '+res(icr)

  ipos = STRPOS(res(icr), '.gz')
  res2 = STRMID(res(icr), 0, ipos)
 
  print, 'PLOT_RADAR_ALL - reading radar file: ', res2 
  READ_SCANCLDRAD_NC, res2(0), time_cr_help, range_help, ze_help, ldr_help, vd_help, sigma_help, elv_help, elvv_help, azi_help, aziv_help
ENDIF

IF date_m+date_d EQ date_sys THEN BEGIN
  print, 'PLOT_RADAR_ALL - reading radar file: ', res(icr) 
  READ_SCANCLDRAD_NC, res(icr), time_cr_help, range_help, ze_help, ldr_help, vd_help, sigma_help, elv_help, elvv_help, azi_help, aziv_help
ENDIF

  time_cr_help = time_cr_help/86400d + JULDAY(1,1,1970,0,0,0)
  CALDAT, time_cr_help, month_cr_help, day_cr_help, year_cr_help, hour_cr_help, minute_cr_help, second_cr_help
  time_cr_help = hour_cr_help + minute_cr_help/60. + second_cr_help/3600.
  IF flag_first_file EQ 0 THEN BEGIN
   range=range_help
  ENDIF ELSE BEGIN
   IF N_ELEMENTS(range_help) GT N_ELEMENTS(range) THEN range=range_help
  ENDELSE

  index_day=WHERE(day_cr_help EQ date_d)
  time_cr_help=time_cr_help(index_day)
  ze_help=ze_help(*,index_day)
  ldr_help=ldr_help(*,index_day)
  vd_help=vd_help(*,index_day)
  sigma_help=sigma_help(*,index_day)
  elv_help=elv_help(index_day)
  elvv_help=elvv_help(index_day)
  azi_help=azi_help(index_day)
  aziv_help=aziv_help(index_day)

      IF flag_first_file EQ 0 THEN BEGIN
      nheightlev=N_ELEMENTS(range_help)
      time_cr=time_cr_help
      ze=ze_help
      ldr=ldr_help
      vd=vd_help
      sigma=sigma_help
      elv=elv_help
      elvv=elvv_help
      azi=azi_help
      aziv=aziv_help
      ENDIF ELSE BEGIN
       nheightlev_help=N_ELEMENTS(range_help)
       IF nheightlev EQ nheightlev_help THEN BEGIN
       time_cr=[time_cr,time_cr_help]
       ze=[[ze],[ze_help]]
       ldr=[[ldr],[ldr_help]]
       vd=[[vd],[vd_help]]
       sigma=[[sigma],[sigma_help]]
       elv=[elv,elv_help]
       elvv=[elvv, elvv_help]
       azi=[azi, azi_help]
       aziv=[aziv,aziv_help]
       ENDIF ELSE BEGIN
        IF nheightlev GT nheightlev_help THEN BEGIN
        ntimeold=N_ELEMENTS(time_cr)
        time_cr=[time_cr,time_cr_help]
        ntimenew=N_ELEMENTS(time_cr)
        fillarray=REPLICATE(-999., nheightlev, N_ELEMENTS(time_cr_help))
        ze=[[ze],[fillarray]]
        ze(0:nheightlev_help-1, ntimeold:ntimenew-1)=ze_help
        ldr=[[ldr],[fillarray]]
        ldr(0:nheightlev_help-1, ntimeold:ntimenew-1)=ldr_help
        vd=[[vd],[fillarray]]
        vd(0:nheightlev_help-1, ntimeold:ntimenew-1)=vd_help
        sigma=[[sigma],[fillarray]]
        sigma(0:nheightlev_help-1, ntimeold:ntimenew-1)=sigma_help
        elv=[elv,elv_help]
        elvv=[elvv, elvv_help]
        azi=[azi, azi_help]
        aziv=[aziv,aziv_help]

        ENDIF
       IF nheightlev LT nheightlev_help THEN BEGIN
        ntimeold=N_ELEMENTS(time_cr)
        time_cr=[time_cr,time_cr_help]
        ntimenew=N_ELEMENTS(time_cr)
        fillarray=REPLICATE(-999., nheightlev_help, ntimenew)
        ze_temp=fillarray
        ze_temp(0:nheightlev-1, 0:ntimeold-1)=ze
        ze_temp(0:nheightlev_help-1, ntimeold:ntimenew-1)=ze_help
        ze=ze_temp
        ldr_temp=fillarray
        ldr_temp(0:nheightlev-1, 0:ntimeold-1)=ldr
        ldr_temp(0:nheightlev_help-1, ntimeold:ntimenew-1)=ldr_help
        ldr=ldr_temp
        vd_temp=fillarray
        vd_temp(0:nheightlev-1, 0:ntimeold-1)=vd
        vd_temp(0:nheightlev_help-1, ntimeold:ntimenew-1)=vd_help
        vd=vd_temp
        sigma_temp=fillarray
        sigma_temp(0:nheightlev-1, 0:ntimeold-1)=sigma
        sigma_temp(0:nheightlev_help-1, ntimeold:ntimenew-1)=sigma_help
        sigma=sigma_temp
        elv=[elv,elv_help]
        elvv=[elvv, elvv_help]
        azi=[azi, azi_help]
        aziv=[aziv,aziv_help]
        nheightlev=nheightlev_help
        ENDIF


       ENDELSE
      ENDELSE



flag_first_file=flag_first_file+1

IF date_m+date_d NE date_sys THEN BEGIN
print, 'zipping radar file: ', res2(0)
SPAWN, 'gzip -f '+res2(0)
ENDIF

 ENDFOR ;end loop over radar files on that day
 

;****create new time grid with temporal resolution of 10 s
time_cr_new=(INDGEN(8640)*10.)/3600.
n_time_cr_new=N_ELEMENTS(time_cr_new)
nd_rad=N_ELEMENTS(range)

;****create new variables on regular grid
ze_new = REPLICATE(-999., nd_rad, n_time_cr_new)
ldr_new =REPLICATE(-999., nd_rad, n_time_cr_new)
vd_new=REPLICATE(-999., nd_rad, n_time_cr_new)
sigma_new=REPLICATE(-999., nd_rad, n_time_cr_new)
elv_new=REPLICATE(-999., n_time_cr_new)
elvv_new=REPLICATE(-999., n_time_cr_new)
azi_new=REPLICATE(-999., n_time_cr_new)
aziv_new=REPLICATE(-999., n_time_cr_new)

;****find nearest value (in time)
tres = 5./3600. ; (5 sec)
count_nomeas=0
FOR i = 0, n_time_cr_new-1 DO BEGIN

ii_f = WHERE(date_filter EQ date)
ii_f = ii_f(0)
IF ii_f NE -1 THEN BEGIN
iii_f = WHERE(time_cr_new(i) GE time_filter_b(ii_f, *) AND time_cr_new(i) LE time_filter_e(ii_f, *))
iii_f = iii_f(0)
ENDIF ELSE iii_f=-1

    IF iii_f NE -1 THEN BEGIN;...check radar manual filter: if time is filtered out then no need to fill in the new array

	  IF count_nomeas EQ 0 THEN BEGIN ;need to remember index of missing measurement so that the right color (grey) is used when plotting
		index_nomeas=[i]
	  ENDIF ELSE BEGIN
		index_nomeas=[index_nomeas,i]
	  ENDELSE
	      count_nomeas=count_nomeas+1

    ENDIF ELSE BEGIN ;else find nearest value 
	minx = MIN(ABS(time_cr_new(i)-time_cr), index)
	IF minx LT tres THEN BEGIN
	  ze_new(*,i)=ze(*,index)
	  ldr_new(*,i)=ldr(*,index)
	  vd_new(*,i)=vd(*,index)
	  sigma_new(*,i)=sigma(*,index)
	  elv_new(i)=elv(index)
	  elvv_new(i)=elvv(index)
	  azi_new(i)=azi(index)
	  aziv_new(i)=aziv(index)
	ENDIF ELSE BEGIN; no nearest value found
	      IF count_nomeas EQ 0 THEN BEGIN ;need to remember index of missing measurement so that the right color (grey) is used when plotting
		index_nomeas=[i]
	      ENDIF ELSE BEGIN
		index_nomeas=[index_nomeas,i]
	      ENDELSE
	      count_nomeas=count_nomeas+1
	ENDELSE
    
    ENDELSE
ENDFOR

;****rename variables
;time_cr=time_cr_new
;ze=ze_new
;ldr=ldr_new
;vd=vd_new
;sigma=sigma_new
;elv=elv_new
;elvv=elvv_new
;azi=azi_new
;aziv=aziv_new

;****plot ze, ldr, vd and sigma time series

pos_lwp = [0.1, 0.85, 0.9, 0.95]
pos_ze = [0.1, 0.65, 0.9, 0.85]
pos_ldr = [0.1, 0.45, 0.9, 0.65]
pos_vd = [0.1, 0.25, 0.9, 0.45]
pos_si = [0.1, 0.05, 0.9, 0.25]
;pos_ze = [0.1, 0.7, 0.9, 0.9]
;pos_ldr = [0.1, 0.5, 0.9, 0.7]
;pos_vd = [0.1, 0.3, 0.9, 0.5]
;pos_si = [0.1, 0.09, 0.9, 0.3]

;***convert Z to dBZ and dB
ze_new = 10d*ALOG10(ze_new)
ldr_new = 10d*ALOG10(ldr_new)

;***parameter ranges
ze_min  = -60.
ze_max  =  20.
ldr_min = -35.
ldr_max =  5.
vd_max =   10.
vd_min =  -10.
sigma_min =  0.
sigma_max =  2.

range=range/1000.

range_max =  14.999
range_min = MIN(range)
range_ind = WHERE(range GE range_min and range LE range_max)
time_ind=WHERE(elv_new LE 89.9 OR elv_new GT 90.1)

ydim = N_ELEMENTS(time_cr_new)
xdim = N_ELEMENTS(range_ind)
arr_ze  = BYTARR(xdim, ydim)
arr_ldr = BYTARR(xdim, ydim)
arr_vd = BYTARR(xdim, ydim)
arr_sigma = BYTARR(xdim, ydim)

;***convert ze to 0...255
arr_ze = (ABS(ze_min) + ze_new(range_ind,*)) / ((ze_max-ze_min)/253)
arr_ldr = (ABS(ldr_min) + ldr_new(range_ind,*)) / ((ldr_max-ldr_min)/253)
arr_vd = (ABS(vd_min) + vd_new(range_ind,*)) / ((vd_max-vd_min)/253);*(-1)+255
arr_sigma = (ABS(sigma_min) + sigma_new(range_ind,*)) / ((sigma_max-sigma_min)/253)


FOR i = 0, ydim-1 DO BEGIN
 ii = WHERE(~FINITE(arr_ze(*, i)) OR arr_ze(*, i) LT 0.)
 IF ii(0) NE -1 THEN arr_ze(ii, i) = 255b
;  arr_ze(*,time_ind) = 255b
 ii = WHERE(~FINITE(arr_vd(*, i)) OR arr_vd(*, i) LT 0.)
 IF ii(0) NE -1 THEN arr_vd(ii, i) = 255b
;  arr_vd(*, time_ind) = 255b
 ii = WHERE(~FINITE(arr_sigma(*, i)) OR arr_sigma(*, i) LT 0.)
 IF ii(0) NE -1 THEN arr_sigma(ii, i) = 255b
;  arr_sigma(*, time_ind) = 255b
  ii = WHERE(~FINITE(arr_ldr(*, i)) OR arr_ldr(*, i) LT 0.)
 IF ii(0) NE -1 THEN arr_ldr(ii, i) = 255b
;  arr_ldr(*, time_ind) = 255b
ENDFOR

;***set missing measurements to fill value
IF count_nomeas NE 0 THEN BEGIN
arr_ze(*,index_nomeas) = 254b
arr_vd(*, index_nomeas) = 254b
arr_sigma(*, index_nomeas) = 254b
arr_ldr(*, index_nomeas) = 254b
ENDIF

;***set off-zenith measurements to fill value since I only want to plot zenith measurements
IF time_ind(0) NE -1 THEN BEGIN
arr_ze(*,time_ind) = 254b
arr_vd(*, time_ind) = 254b
arr_sigma(*, time_ind) = 254b
arr_ldr(*, time_ind) = 254b 
ENDIF

 ;****plot radar&LWP time series
!P.charsize = 1.8

psfile = filename1_plot_radar
set_plot, 'ps'
!p.multi = [0, 1, 4]
LOADCT, 1
device, /color, file = psfile, BITS_PER_PIXEL = 8, /landscape;, xsize=29.6, ysize=26.
!P.font = 3

;**plot ze
LOADCT, 34
TVLCT, red, green, blue, /GET
red(255)=255
green(255)=255
blue(255)=255
red(254)=211
green(254)=211
blue(254)=211
red(0)=0
green(0)=0
blue(0)=0
TVLCT, red, green, blue
TV, TRANSPOSE(arr_ze), pos_ze(0), pos_ze(1), xsize=pos_ze(2)-pos_ze(0), ysize=pos_ze(3)-pos_ze(1), /NORMAL
LOADCT, 0
PLOT, time_cr_new, range,/NODATA,$
      yrange=[range(0), range_max], ytitle='Altitude [km AGL]',$
      xrange=[time_cr_new(0), time_cr_new(ydim-1)], xstyle=1, ystyle=1,  $
      xcharsize=1e-5, /NOERASE, position = pos_ze, yticklen = -0.009, xticklen = 0.04

n_ticks = 4
z = FLTARR(n_ticks+1)
FOR i = 0, n_ticks DO BEGIN
 z(i) = ze_min + (ze_max-ze_min)*i/n_ticks
ENDFOR

LOADCT, 34
TVLCT, red, green, blue, /GET
red(255)=255
green(255)=255
blue(255)=255
red(254)=211
green(254)=211
blue(254)=211
red(0)=0
green(0)=0
blue(0)=0
TVLCT, red, green, blue
COLORBAR1, position = [pos_ze[2]+0.02, pos_ze[1], pos_ze[2]+0.04, pos_ze[3]-(pos_ze[3]-pos_ze[1])/10.],$
          ncolors=256, divisions=n_ticks, ticknames=STRING(z, FORMAT='(I4)'),$
          title = 'Ze [dBZe]', charsize=1.5, color=0 , font=0, right=1, /vertical
XYOUTS, pos_ze(0)+0.02, pos_ze(3)-0.03, 'Radar reflectivity factor', /normal, charsize=1.

;*** plot Doppler velocity
LOADCT, 18

MAKE_BGYR_CT,vd_min,vd_max,red,green,blue
TVLCT, red, green, blue


TV, TRANSPOSE(arr_vd), pos_vd(0), pos_vd(1), xsize=pos_vd(2)-pos_vd(0), ysize=pos_vd(3)-pos_vd(1), /NORMAL
LOADCT, 0
PLOT, time_cr_new, range,/NODATA,$
      yrange=[range(0), range_max], ytitle='Altitude [km AGL]',$
      xrange=[time_cr_new(0), time_cr_new(ydim-1)], xstyle=1, ystyle=1,  $
      xcharsize=1e-5, /NOERASE, position = pos_vd, yticklen = -0.009, xticklen = 0.04

n_ticks = 4
z = FLTARR(n_ticks+1)
FOR i = 0, n_ticks DO BEGIN
 z(i) = vd_min + (vd_max-vd_min)*i/n_ticks
ENDFOR


MAKE_BGYR_CT,vd_min,vd_max,red,green,blue
TVLCT, red, green, blue
COLORBAR1, position = [pos_vd[2]+0.02, pos_vd[1], pos_vd[2]+0.04, pos_vd[3]-(pos_vd[3]-pos_vd[1])/10.],$
          ncolors=256, divisions=n_ticks, ticknames=STRING(z, FORMAT='(f5.1)'),$
          title = 'Dop. vel. [m/s]', charsize=1.5, color=0 , font=0, right=1, /vertical
XYOUTS, pos_vd(0)+0.02, pos_vd(3)-0.03, 'Doppler velocity', /normal, charsize=1.

;*** plot sigma
LOADCT, 34
TVLCT, red, green, blue, /GET
red(255)=255
green(255)=255
blue(255)=255
red(254)=211
green(254)=211
blue(254)=211
red(0)=0
green(0)=0
blue(0)=0
TVLCT, red, green, blue
TV, TRANSPOSE(arr_sigma), pos_si(0), pos_si(1), xsize=pos_si(2)-pos_si(0), ysize=pos_si(3)-pos_si(1), /NORMAL
LOADCT, 0
PLOT, time_cr_new, range,/NODATA,$
      yrange=[range(0), range_max], ytitle='Altitude [km AGL]',$
      xrange=[time_cr_new(0), time_cr_new(ydim-1)], xstyle=1, ystyle=1,  $
      /NOERASE, position = pos_si, yticklen = -0.009, xticklen = 0.04,$
      xtitle = 'Time [UTC] on ' +date
      
n_ticks = 4
z = FLTARR(n_ticks+1)
FOR i = 0, n_ticks DO BEGIN
 z(i) = sigma_min + (sigma_max-sigma_min)*i/n_ticks
ENDFOR

LOADCT, 34
TVLCT, red, green, blue, /GET
red(255)=255
green(255)=255
blue(255)=255
red(254)=211
green(254)=211
blue(254)=211
red(0)=0
green(0)=0
blue(0)=0
TVLCT, red, green, blue
COLORBAR1, position = [pos_si[2]+0.02, pos_si[1], pos_si[2]+0.04, pos_si[3]-(pos_si[3]-pos_si[1])/10.],$
          ncolors=256, divisions=n_ticks, ticknames=STRING(z, FORMAT='(F4.1)'),$
          title = 'Spectral width [m/s]', charsize=1.5, color=0 , font=0, right=1, /vertical
XYOUTS, pos_si(0)+0.02, pos_si(3)-0.03, 'Spectral width', /normal, charsize=1.

;*** plot Linear depolarization ratio
LOADCT, 34
TVLCT, red, green, blue, /GET
red(255)=255
green(255)=255
blue(255)=255
red(254)=211
green(254)=211
blue(254)=211
red(0)=0
green(0)=0
blue(0)=0
TVLCT, red, green, blue
TV, TRANSPOSE(arr_ldr), pos_ldr(0), pos_ldr(1), xsize=pos_ldr(2)-pos_ldr(0), ysize=pos_ldr(3)-pos_ldr(1), /NORMAL
LOADCT, 0

PLOT, time_cr_new, range,/NODATA,$
      yrange=[range(0), range_max], ytitle='Altitude [km AGL]',$
      xrange=[time_cr_new(0), time_cr_new(ydim-1)], xstyle=1, ystyle=1,  $
      xcharsize=1e-5, /NOERASE, position = pos_ldr, yticklen = -0.009, xticklen = 0.04

n_ticks = 4
z = FLTARR(n_ticks+1)
FOR i = 0, n_ticks DO BEGIN
 z(i) = ldr_min + (ldr_max-ldr_min)*i/n_ticks
ENDFOR

LOADCT, 34
TVLCT, red, green, blue, /GET
red(255)=255
green(255)=255
blue(255)=255
red(254)=211
green(254)=211
blue(254)=211
red(0)=0
green(0)=0
blue(0)=0
TVLCT, red, green, blue
COLORBAR1, position = [pos_ldr[2]+0.02, pos_ldr[1], pos_ldr[2]+0.04, pos_ldr[3]-(pos_ldr[3]-pos_ldr[1])/10.],$
          ncolors=256, divisions=n_ticks, ticknames=STRING(z, FORMAT='(f5.1)'),$
          title = 'LDR [dB]', charsize=1.5, color=0 , font=0, right=1, /vertical
XYOUTS, pos_ldr(0)+0.02, pos_ldr(3)-0.03, 'Linear depolarization ratio', /normal, charsize=1.

;***plot LWP time series
 

 
  ii = WHERE(lwp GT -50.)
  tit = 'TOPHAT LWP & cloud radar zenith measurements, ' + station +', ' + date
 
 
  PLOT, [0, 0], [24, 1], /NODATA, yrange=[lwp_min, lwp_max], ystyle = 1, $
        xrange = [time_cr_new(0), time_cr_new(ydim-1)], xstyle = 1, xcharsize = 0.0000001, color=0, title=tit,$
        position = pos_lwp, xticks = 1, ytitle = 'LWP [gm!e-2!N]', /NOERASE, /NORMAL
      
 
  OPLOT, time_w, lwp, psym = 3, color = 50
  ;OPLOT, [time_w(0), time_w(N_ELEMENTS(time_w)-1)], [0., 0.], linestyle = 1
  XYOUTS, pos_lwp(0)+0.02, pos_lwp(3)-0.03, 'Liquid water path', /normal, charsize=1.


DEVICE, /CLOSE


;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;plot all radar data
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;****create new time grid with temporal resolution of 10 s
time_cr_new=(INDGEN(8640)*10.)/3600.
n_time_cr_new=N_ELEMENTS(time_cr_new)
nd_rad=N_ELEMENTS(range)

;****create new variables on regular grid
ze_new = REPLICATE(-999., nd_rad, n_time_cr_new)
ldr_new =REPLICATE(-999., nd_rad, n_time_cr_new)
vd_new=REPLICATE(-999., nd_rad, n_time_cr_new)
sigma_new=REPLICATE(-999., nd_rad, n_time_cr_new)
elv_new=REPLICATE(-999., n_time_cr_new)
elvv_new=REPLICATE(-999., n_time_cr_new)
azi_new=REPLICATE(-999., n_time_cr_new)
aziv_new=REPLICATE(-999., n_time_cr_new)

;****find nearest value (in time)
tres = 5./3600. ; (5 sec)
count_nomeas=0
FOR i = 0, n_time_cr_new-1 DO BEGIN

ii_f = WHERE(date_filter EQ date)
ii_f = ii_f(0)
IF ii_f NE -1 THEN BEGIN
iii_f = WHERE(time_cr_new(i) GE time_filter_b(ii_f, *) AND time_cr_new(i) LE time_filter_e(ii_f, *))
iii_f = iii_f(0)
ENDIF ELSE iii_f=-1

    IF iii_f NE -1 THEN BEGIN;...check radar manual filter: if time is filtered out then no need to fill in the new array

	  IF count_nomeas EQ 0 THEN BEGIN ;need to remember index of missing measurement so that the right color (grey) is used when plotting
		index_nomeas=[i]
	  ENDIF ELSE BEGIN
		index_nomeas=[index_nomeas,i]
	  ENDELSE
	      count_nomeas=count_nomeas+1

    ENDIF ELSE BEGIN ;else find nearest value 
	minx = MIN(ABS(time_cr_new(i)-time_cr), index)
	IF minx LT tres THEN BEGIN
	  ze_new(*,i)=ze(*,index)
	  ldr_new(*,i)=ldr(*,index)
	  vd_new(*,i)=vd(*,index)
	  sigma_new(*,i)=sigma(*,index)
	  elv_new(i)=elv(index)
	  elvv_new(i)=elvv(index)
	  azi_new(i)=azi(index)
	  aziv_new(i)=aziv(index)
	ENDIF ELSE BEGIN; no nearest value found
	      IF count_nomeas EQ 0 THEN BEGIN ;need to remember index of missing measurement so that the right color (grey) is used when plotting
		index_nomeas=[i]
	      ENDIF ELSE BEGIN
		index_nomeas=[index_nomeas,i]
	      ENDELSE
	      count_nomeas=count_nomeas+1
	ENDELSE
    
    ENDELSE
ENDFOR

;****rename variables
;time_cr=time_cr_new
;ze=ze_new
;ldr=ldr_new
;vd=vd_new
;sigma=sigma_new
;elv=elv_new
;elvv=elvv_new
;azi=azi_new
;aziv=aziv_new

;****plot ze, ldr, vd and sigma time series

pos_lwp = [0.1, 0.85, 0.9, 0.95]
pos_ze = [0.1, 0.65, 0.9, 0.85]
pos_ldr = [0.1, 0.45, 0.9, 0.65]
pos_vd = [0.1, 0.25, 0.9, 0.45]
pos_si = [0.1, 0.05, 0.9, 0.25]
;pos_ze = [0.1, 0.7, 0.9, 0.9]
;pos_ldr = [0.1, 0.5, 0.9, 0.7]
;pos_vd = [0.1, 0.3, 0.9, 0.5]
;pos_si = [0.1, 0.09, 0.9, 0.3]

;***convert Z to dBZ and dB
ze_new = 10d*ALOG10(ze_new)
ldr_new = 10d*ALOG10(ldr_new)

;***parameter ranges
ze_min  = -60.
ze_max  =  20.
ldr_min = -35.
ldr_max =  5.
vd_max =   10.
vd_min =  -10.
sigma_min =  0.
sigma_max =  2.

range=range/1000.

range_max =  14.999
range_min = MIN(range)
range_ind = WHERE(range GE range_min and range LE range_max)
time_ind=WHERE(elv NE 90.)

ydim = N_ELEMENTS(time_cr_new)
xdim = N_ELEMENTS(range_ind)
arr_ze  = BYTARR(xdim, ydim)
arr_ldr = BYTARR(xdim, ydim)
arr_vd = BYTARR(xdim, ydim)
arr_sigma = BYTARR(xdim, ydim)

;***convert ze to 0...255
arr_ze = (ABS(ze_min) + ze_new(range_ind,*)) / ((ze_max-ze_min)/253)
arr_ldr = (ABS(ldr_min) + ldr_new(range_ind,*)) / ((ldr_max-ldr_min)/253)
arr_vd = (ABS(vd_min) + vd_new(range_ind,*)) / ((vd_max-vd_min)/253);*(-1)+255
arr_sigma = (ABS(sigma_min) + sigma_new(range_ind,*)) / ((sigma_max-sigma_min)/253)


FOR i = 0, ydim-1 DO BEGIN
 ii = WHERE(~FINITE(arr_ze(*, i)) OR arr_ze(*, i) LT 0.)
 IF ii(0) NE -1 THEN arr_ze(ii, i) = 255b
;  arr_ze(*,time_ind) = 255b
 ii = WHERE(~FINITE(arr_vd(*, i)) OR arr_vd(*, i) LT 0.)
 IF ii(0) NE -1 THEN arr_vd(ii, i) = 255b
;  arr_vd(*, time_ind) = 255b
 ii = WHERE(~FINITE(arr_sigma(*, i)) OR arr_sigma(*, i) LT 0.)
 IF ii(0) NE -1 THEN arr_sigma(ii, i) = 255b
;  arr_sigma(*, time_ind) = 255b
  ii = WHERE(~FINITE(arr_ldr(*, i)) OR arr_ldr(*, i) LT 0.)
 IF ii(0) NE -1 THEN arr_ldr(ii, i) = 255b
;  arr_ldr(*, time_ind) = 255b
ENDFOR

;***set missing measurements to fill value
IF count_nomeas NE 0 THEN BEGIN
arr_ze(*,index_nomeas) = 254b
arr_vd(*, index_nomeas) = 254b
arr_sigma(*, index_nomeas) = 254b
arr_ldr(*, index_nomeas) = 254b
ENDIF


arr_ze_all=arr_ze
arr_vd_all=arr_vd
arr_sigma_all=arr_sigma
arr_ldr_all=arr_ldr


;*****plot all radar data

set_plot, 'ps'
!P.charsize = 1.8
!p.multi = [0, 1, 4]
TEK_COLOR
device, /color, file = filename2_plot_radar, BITS_PER_PIXEL = 8, /landscape;, xsize=29.6, ysize=26.
!P.font = 3
tit = 'Cloud radar measurements, ' +station+ ', ' + date

;****plot ze, ldr, vd and sigma time series

pos_azi = [0.1, 0.85, 0.9, 0.95]
pos_elv = [0.1, 0.75, 0.9, 0.85]
pos_ze = [0.1, 0.575, 0.9, 0.75]
pos_ldr = [0.1, 0.4, 0.9, 0.575]
pos_vd = [0.1, 0.225, 0.9, 0.4]
pos_si = [0.1, 0.05, 0.9, 0.225]
;pos_ze = [0.1, 0.7, 0.9, 0.9]
;pos_ldr = [0.1, 0.5, 0.9, 0.7]
;pos_vd = [0.1, 0.3, 0.9, 0.5]
;pos_si = [0.1, 0.09, 0.9, 0.3]


;****plot azimuth angle
  PLOT, [0, 0], [24, 1], /NODATA, yrange=[0, 360], ystyle = 1,$
        xrange = [time_cr_new(0), time_cr_new(ydim-1)], xstyle = 1, xcharsize = 0.0000001, color=0, title=tit,$
        position = pos_azi, xticks = 1, ytitle = 'Azi.[deg]', /NOERASE, /NORMAL


  OPLOT, time_cr_new, azi_new, color = 0, psym=2, symsize=0.4
  OPLOT, [time_cr_new(0), time_cr_new(ydim-1)], [136.5, 136.5], linestyle=1, color=2 ;North
  OPLOT, [time_cr_new(0), time_cr_new(ydim-1)], [226.5, 226.5], linestyle=1, color=2;East
  OPLOT, [time_cr_new(0), time_cr_new(ydim-1)], [316.5, 316.5], linestyle=1, color=2 ;South
  OPLOT, [time_cr_new(0), time_cr_new(ydim-1)], [46.5, 46.5], linestyle=1, color=2 ;West

  ;OPLOT, [time_w(0), time_w(N_ELEMENTS(time_w)-1)], [0., 0.], linestyle = 1
  XYOUTS, 24.1, 136.5, 'North', charsize=0.7;, alignment=0.5
  XYOUTS, 24.1, 226.5, 'East', charsize=0.7;, alignment=0.5
  XYOUTS, 24.1, 316.5, 'South', charsize=0.7;, alignment=0.5
  XYOUTS, 24.1, 46.5, 'West', charsize=0.7;, alignment=0.5

;****plot elevation angle
  PLOT, [0, 0], [24, 1], /NODATA, yrange=[0, 180], ystyle = 1,$
        xrange = [time_cr_new(0), time_cr_new(ydim-1)], xstyle = 1, xcharsize = 0.0000001, color=0,$
        position = pos_elv, xticks = 1, ytitle = 'Elev.[deg]', /NOERASE, /NORMAL


  OPLOT, time_cr_new, elv_new, color = 0, psym=2, symsize=0.4
  OPLOT, [time_cr_new(0), time_cr_new(ydim-1)], [90., 90.], linestyle=1, color=2 ;Zenith
  ;OPLOT, [time_w(0), time_w(N_ELEMENTS(time_w)-1)], [0., 0.], linestyle = 1
  ;XYOUTS, pos_elv(0)+0.02, pos_elv(3)-0.03, 'Elevation angle', /normal, charsize=1.


;**plot ze
LOADCT, 34
TVLCT, red, green, blue, /GET
red(255)=255
green(255)=255
blue(255)=255
red(254)=211
green(254)=211
blue(254)=211
red(0)=0
green(0)=0
blue(0)=0
TVLCT, red, green, blue
TV, TRANSPOSE(arr_ze_all), pos_ze(0), pos_ze(1), xsize=pos_ze(2)-pos_ze(0), ysize=pos_ze(3)-pos_ze(1), /NORMAL
LOADCT, 0
PLOT, time_cr_new, range,/NODATA,$
      yrange=[range(0), range_max], ytitle='Altitude [km AGL]',$
      xrange=[time_cr_new(0), time_cr_new(ydim-1)], xstyle=1, ystyle=1,$
      xcharsize=1e-5, /NOERASE, position = pos_ze, yticklen = -0.009, xticklen = 0.04

n_ticks = 4
z = FLTARR(n_ticks+1)
FOR i = 0, n_ticks DO BEGIN
 z(i) = ze_min + (ze_max-ze_min)*i/n_ticks
ENDFOR

LOADCT, 34
TVLCT, red, green, blue, /GET
red(255)=255
green(255)=255
blue(255)=255
red(254)=211
green(254)=211
blue(254)=211
red(0)=0
green(0)=0
blue(0)=0
TVLCT, red, green, blue
COLORBAR1, position = [pos_ze[2]+0.02, pos_ze[1], pos_ze[2]+0.04, pos_ze[3]-(pos_ze[3]-pos_ze[1])/10.],$
          ncolors=256, divisions=n_ticks, ticknames=STRING(z, FORMAT='(I4)'),$
          title = 'Ze [dBZe]', charsize=1.5, color=0 , font=0, right=1, /vertical
XYOUTS, pos_ze(0)+0.02, pos_ze(3)-0.03, 'Radar reflectivity factor', /normal, charsize=1.

;*** plot Doppler velocity
LOADCT, 18

MAKE_BGYR_CT,vd_min,vd_max,red,green,blue
TVLCT, red, green, blue


TV, TRANSPOSE(arr_vd_all), pos_vd(0), pos_vd(1), xsize=pos_vd(2)-pos_vd(0), ysize=pos_vd(3)-pos_vd(1), /NORMAL
LOADCT, 0
PLOT, time_cr_new, range,/NODATA,$
      yrange=[range(0), range_max], ytitle='Altitude [km AGL]',$
      xrange=[time_cr_new(0), time_cr_new(ydim-1)], xstyle=1, ystyle=1, $
      xcharsize=1e-5, /NOERASE, position = pos_vd, yticklen = -0.009, xticklen = 0.04

n_ticks = 4
z = FLTARR(n_ticks+1)
FOR i = 0, n_ticks DO BEGIN
 z(i) = vd_min + (vd_max-vd_min)*i/n_ticks
ENDFOR


MAKE_BGYR_CT,vd_min,vd_max,red,green,blue
TVLCT, red, green, blue
COLORBAR1, position = [pos_vd[2]+0.02, pos_vd[1], pos_vd[2]+0.04, pos_vd[3]-(pos_vd[3]-pos_vd[1])/10.],$
          ncolors=256, divisions=n_ticks, ticknames=STRING(z, FORMAT='(f5.1)'),$
          title = 'Dop. vel. [m/s]', charsize=1.5, color=0 , font=0, right=1, /vertical
XYOUTS, pos_vd(0)+0.02, pos_vd(3)-0.03, 'Doppler velocity', /normal, charsize=1.


;*** plot sigma
LOADCT, 34
TVLCT, red, green, blue, /GET
red(255)=255
green(255)=255
blue(255)=255
red(254)=211
green(254)=211
blue(254)=211
red(0)=0
green(0)=0
blue(0)=0
TVLCT, red, green, blue
TV, TRANSPOSE(arr_sigma_all), pos_si(0), pos_si(1), xsize=pos_si(2)-pos_si(0), ysize=pos_si(3)-pos_si(1), /NORMAL
LOADCT, 0
PLOT, time_cr_new, range,/NODATA,$
      yrange=[range(0), range_max], ytitle='Altitude [km AGL]',$
      xrange=[time_cr_new(0), time_cr_new(ydim-1)], xstyle=1, ystyle=1,$
      /NOERASE, position = pos_si, yticklen = -0.009, xticklen = 0.04,$
      xtitle = 'Time [UTC] on ' +date 

n_ticks = 4
z = FLTARR(n_ticks+1)
FOR i = 0, n_ticks DO BEGIN
 z(i) = sigma_min + (sigma_max-sigma_min)*i/n_ticks
ENDFOR

LOADCT, 34
TVLCT, red, green, blue, /GET
red(255)=255
green(255)=255
blue(255)=255
red(254)=211
green(254)=211
blue(254)=211
red(0)=0
green(0)=0
blue(0)=0
TVLCT, red, green, blue
COLORBAR1, position = [pos_si[2]+0.02, pos_si[1], pos_si[2]+0.04, pos_si[3]-(pos_si[3]-pos_si[1])/10.],$
          ncolors=256, divisions=n_ticks, ticknames=STRING(z, FORMAT='(F4.1)'),$
          title = 'Spectral width [m/s]', charsize=1.5, color=0 , font=0, right=1, /vertical
XYOUTS, pos_si(0)+0.02, pos_si(3)-0.03, 'Spectral width', /normal, charsize=1.

;*** plot Linear depolarization ratio
LOADCT, 34
TVLCT, red, green, blue, /GET
red(255)=255
green(255)=255
blue(255)=255
red(254)=211
green(254)=211
blue(254)=211
red(0)=0
green(0)=0
blue(0)=0
TVLCT, red, green, blue
TV, TRANSPOSE(arr_ldr_all), pos_ldr(0), pos_ldr(1), xsize=pos_ldr(2)-pos_ldr(0), ysize=pos_ldr(3)-pos_ldr(1), /NORMAL
LOADCT, 0

PLOT, time_cr_new, range,/NODATA,$
      yrange=[range(0), range_max], ytitle='Altitude [km AGL]',$
      xrange=[time_cr_new(0), time_cr_new(ydim-1)], xstyle=1, ystyle=1,$
      xcharsize=1e-5, /NOERASE, position = pos_ldr, yticklen = -0.009, xticklen = 0.04

n_ticks = 4
z = FLTARR(n_ticks+1)
FOR i = 0, n_ticks DO BEGIN
 z(i) = ldr_min + (ldr_max-ldr_min)*i/n_ticks
ENDFOR

LOADCT, 34
TVLCT, red, green, blue, /GET
red(255)=255
green(255)=255
blue(255)=255
red(254)=211
green(254)=211
blue(254)=211
red(0)=0
green(0)=0
blue(0)=0
TVLCT, red, green, blue
COLORBAR1, position = [pos_ldr[2]+0.02, pos_ldr[1], pos_ldr[2]+0.04, pos_ldr[3]-(pos_ldr[3]-pos_ldr[1])/10.],$
          ncolors=256, divisions=n_ticks, ticknames=STRING(z, FORMAT='(f5.1)'),$
          title = 'LDR [dB]', charsize=1.5, color=0 , font=0, right=1, /vertical
XYOUTS, pos_ldr(0)+0.02, pos_ldr(3)-0.03, 'Linear depolarization ratio', /normal, charsize=1.


DEVICE, /CLOSE


ENDIF ELSE BEGIN

print, 'no radar files found for ', date

ENDELSE

END
