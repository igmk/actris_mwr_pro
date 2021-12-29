;+
;*****************
PRO PLOT_CONT_TS,$
;*****************
;*INPUT
par1,$                 ; parameter to be plotted with contours
par2,$                 ; parameter to be plotted with time series
max_p,$                ; maximum par1 plotting range
min_p,$                ; minimum par1 plotting range
ymax,$                 ; maximum par2 plotting range
ymin,$                 ; minimum par2 plotting range
time1,$                ; time axis of plot par1
time2,$                ; time axis of plot par2    
flag,$                 ; data quality flags par2
z_final,$              ; height axis of plot par1
h_min,$                ; minimum height to be plotted
h_max,$                ; maximum height to be plotted
time_start,$           ; minimum time to be plotted
time_end,$             ; maximum time to be plotted
dx,$                   ; contour lines distance (max_p-min_p)/dx
bar_tit,$              ; title color bar
ts_tit,$               ; yaxis title par2
date,$                 ; data yymmdd
pos1,$                 ; position of contour plot
pos2,$                 ; position of color bar
pos3,$                 ; position of flag plot
pos4,$                 ; position of ts plot
tit,$                  ; title of plot
dt2,$                  ; temporal resolution for plotting
home_path              ; path where color table is located

;*OUTPUT
;*KEYWORDS
; Abstract:
;* plot time series of 2D variable (time-height) + time series of a 1D variable
; Author:
; U. Loehnert
; Date:
; 2011-12-13
; Dependencies:
; -
; Changes:
; XXXX-XX-XX: ???
;-

!P.charsize = 1.0
!p.multi = [0, 1, 1]
!P.thick = 3

n_Tz = WHERE(z_final LE h_max)
nz = N_ELEMENTS(n_Tz)

;**contour plot
PLOT, [0, 1], [0, 1], yrange = [h_min/1000., h_max/1000.], xrange = [time_start, time_end],$
      color = 0, ycharsize = 0.7, xtitle = 'hours on '+date+' [UTC]',$
      charsize = 1.2, xstyle = 1, ystyle = 1, ytitle = 'Height [km]',xtickinterval=6,$
      Position = pos1, /Nodata, /Normal, font = 3, xticklen = -0.023, yticklen = -0.013, xcharsize = 0.8

i_valid = WHERE(par1(0, *) GT 0.)
dxx = (max_p - min_p)/dx
IF dxx LT 1 THEN dxx = 1
ll = FINDGEN(FIX(dxx))*dx + min_p
lll = N_ELEMENTS(ll)
clab = REPLICATE(1, lll)
ccol = FLTARR(lll)
FOR ic = 0, lll-1 DO BEGIN
 ccol(ic) = DBZ_2_BYTE_T(FLOAT(ll(ic)), max_p, min_p) 
ENDFOR

IF i_valid(0) NE -1 THEN BEGIN
 nn = N_ELEMENTS(i_valid)

;***create new array on dt2 resolution
 time_new = time_start + FINDGEN(((time_end-time_start)*3600./dt2)+0)*dt2/(60.*60.)
 n_new = N_ELEMENTS(time_new)

;***initialize plot arrays
 time_cont = [time_start]
 p_cont = [par1(n_Tz, i_valid(0))]
 flag_cont = [flag(i_valid(0))]

 j = 1
 FOR k = 0l, n_new-1l DO BEGIN
  
;***deltat to look for nearest neighbour
;  CLOSEST, time_new(k), time1, dmin, itx
;  IF dmin LE 24.*dt2/(86400.*2.) THEN BEGIN    

;   time_cont = [time_cont, time1(itx)]   
;   p_cont = [[p_cont], [par1(n_Tz, itx)]]
;   flag_cont = [flag_cont, flag(itx)]
;  ENDIF 

  xx = 24.*dt2/(86400.*2.)
  ii = WHERE(time1 LT time_new(k)+xx AND time1 GT time_new(k)-xx, nii)

  IF nii EQ 1 THEN BEGIN
   time_cont = [time_cont, time_new(k)]
   p_cont = [[p_cont], [par1(n_Tz, ii(0))]]
   flag_cont = [flag_cont, flag(ii(0))]
  ENDIF ELSE IF nii GT 1 THEN BEGIN
   time_cont = [time_cont, time_new(k)]
   flag_cont = [flag_cont, MAX(flag(ii))]
   x = AVERAGE(par1(0:nz-1, ii), 2)
   p_cont = [[p_cont], [x]]
  ENDIF

 ENDFOR 

 time_cont = [time_cont, time1(nn-1)]
 p_cont = [[p_cont], [par1(n_Tz, nn-1)]]
 flag_cont = [flag_cont, flag(nn-1)]

 n_p_cont_x = N_ELEMENTS(p_cont(*, 0))
 n_p_cont_y = N_ELEMENTS(p_cont(0, *)) 

 FOR i = 0, n_p_cont_x-1 DO BEGIN
  FOR j = 0, n_p_cont_y-1 DO BEGIN
   IF p_cont(i, j) LT min_p THEN p_cont(i, j) = min_p + 1e-4
  ENDFOR
 ENDFOR

 CONTOUR, TRANSPOSE(p_cont), time_cont, z_final(n_Tz)/1000.,$                                                                                                      
          Position = pos1, xrange = [time_start, time_end],$
          Levels = ll, xstyle = 5, ystyle = 5,$
          /Noerase, /FOLLOW, color = 0, /FILL, C_color = ccol

 CONTOUR, TRANSPOSE(p_cont), time_cont, z_final(n_Tz)/1000.,$
          Position = pos1, xrange = [time_start, time_end],$
          Levels = ll, C_labels = clab, xstyle = 5, ystyle = 5,$
          /Noerase, /FOLLOW, color = 0

ENDIF
ncolor = 254
	BarColor = BINDGEN(ncolor)
BarLabel = [min_p, max_p]
!P.Color = 0

IF max_p GE 100. THEN BEGIN

 COLLABEL2, BarColor, BarLabel, $
            CHARSIZE = 0.5, $
            DIRECTION = 1, FONT = 3, $
            LabFormat = '(I3)',LABOFFSET = 35, $
            POSITION = pos2

ENDIF ELSE BEGIN
 COLLABEL2, BarColor, BarLabel, $
            CHARSIZE = 0.5, $
            DIRECTION = 1, FONT = 3, $
            LabFormat = '(f4.1)',LABOFFSET = 35, $
            POSITION = pos2
ENDELSE

XYOUTS, pos2(0)+(pos2(2)-pos2(0))/2., pos2(1)-0.05, bar_tit, color = 0, font = 3, charsize = 0.9, alignment=0.5, /Normal

;**1D time series
PLOT, [0, 1], [0, 1], yrange = [ymin, ymax], xrange = [time_start, time_end],$
      color = 0, ycharsize = 0.7, xcharsize = 1e-5,xtickinterval=6,$
      charsize = 1.2, xstyle = 1, ystyle = 1, ytitle = ts_tit,$
      Position = pos4, /NODATA, /NOERASE, /NORMAL, font = 3, ticklen = -0.013

OPLOT, time2, par2, psym = 3, color = 0 

;**data flag

X = [-1, 1, 1, -1]
Y = [-1, -1, 1, 1]
USERSYM, X, Y, /FILL

LOADCOL, home_path+'/mwr_pro/source/col1'
PLOT, [0, 0], [24, 1], /NODATA, yrange=[0.9,1.1], ystyle = 1, title = tit,xtickinterval=6,$
      xrange = [time_start, time_end], xstyle = 1, xcharsize = 0.0000001, ycharsize = 0.00001,$
      position = pos3, yticks = 1, xticks = 1, /NOERASE, /NORMAL

IF N_ELEMENTS(flag_cont) GT 0 THEN BEGIN

 index = WHERE(flag_cont GT 0, count)
 IF count GT 0 THEN BEGIN
  IF count EQ 1 THEN BEGIN
   count = 2
   index = [index(0), index(0)]
  ENDIF
  h = REPLICATE(1, count)
  OPLOT, time_cont(index), h, psym=8, thick=5, color=7
 ENDIF

 index = WHERE(flag_cont EQ 0, count)
 IF count GT 0 THEN BEGIN
  IF count EQ 1 THEN BEGIN
   count = 2
   index = [index(0), index(0)]
  ENDIF
  h = REPLICATE(1, count)
  OPLOT, time_cont(index), h, psym=8, thick=5, color=4
 ENDIF

 XYOUTS, 0.05, pos3(1)+(pos3(3)-pos3(1))/2., 'DATA FLAG', alignment=0.5, charsize=0.6, /Normal

ENDIF

RETURN

END