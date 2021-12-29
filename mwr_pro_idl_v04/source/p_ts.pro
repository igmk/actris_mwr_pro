PRO P_TS, par, max_p, min_p, time, z_final,$
          h_min, h_max, time_start, time_end, station,$
          base, top, dx, bar_tit, date, var_c

;plot
;time series of 2D variable (time-height)
;UL 02/07

!P.charsize = 1.0
!p.multi = [0, 1, 1]
!P.thick = 3

dummy = LABEL_DATE(DATE_FORMAT=['%H'])

n_Tz = WHERE(z_final/1000. LE h_max)
nz = N_ELEMENTS(n_Tz)

PLOT, time, base, YRANGE = [h_min, h_max], xrange = [time_start, time_end],$
      color = 0, ycharsize = 0.7, title = station + ', IPT '+var_c, xtitle = 'hours on '+date+' [UTC]',$
      charsize = 1.2, xstyle = 1, ystyle = 1, ytitle = 'Height [km]',$
      Position = pos1, /Nodata, /Normal, font = 3, ticklen = -0.013, xcharsize = 0.8,$
      xtickunits= 'hours', xtickformat = 'LABEL_DATE', xminor = 2, xtickinterval = 2, xticks = 12


i_valid = WHERE(par(0, *) GT 0.)

dxx = (max_p - min_p)/dx
ll = FINDGEN(FIX(dxx))*dx + min_p
lll = N_ELEMENTS(ll)
clab = REPLICATE(1, lll)
ccol = FLTARR(lll)
FOR ic = 0, lll-1 DO BEGIN
 ccol(ic) = DBZ_2_BYTE_T(FLOAT(ll(ic)), max_p, min_p) 
ENDFOR

IF i_valid(0) NE -1 THEN BEGIN

 nn = N_ELEMENTS(i_valid)
 p_cont = FLTARR(nz, nn)

 FOR k = 0, nn-1 DO p_cont(*, k) = par(0:nz-1, i_valid(k))

help, p_cont
help, time(i_valid)
help, z_final(n_Tz)/1000.
 CONTOUR, TRANSPOSE(p_cont), time(i_valid), z_final(n_Tz)/1000.,$                                                                                                      
          Position = pos1, xrange = [time_start, time_end],$
          Levels = ll, xstyle = 5, ystyle = 5,$
          /Noerase, /FOLLOW, color = 0, /FILL, C_color = ccol

 CONTOUR, TRANSPOSE(p_cont), time(i_valid), z_final(n_Tz)/1000.,$
          Position = pos1, xrange = [time_start, time_end],$
          Levels = ll, C_labels = clab, xstyle = 5, ystyle = 5,$
          /Noerase, /FOLLOW, color = 0

ENDIF

OPLOT, time, base/1000., psym = 4, color = 0, symsize = 0.4
OPLOT, time, top/1000., psym = 5, color = 0, symsize = 0.4

ncolor = 254
BarColor = BINDGEN(ncolor)
BarLabel = [min_p, max_p]
!P.Color = 0

COLLABEL2, BarColor, BarLabel, $
           CHARSIZE = 0.5, $
           DIRECTION = 1, FONT = 3, $
           LabFormat = '(I3)',LABOFFSET = 35, $
           POSITION = [0.90, 0.5, 0.95, 0.85]

XYOUTS, 0.90,0.46, bar_tit, color = 0, font = 3, charsize = 0.9, /Normal

RETURN

END