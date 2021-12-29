;+
;*************
PRO COSMO_BIAS;, $
;*************

station = 'jue';, $
abscal  = '100819'

;-

outfile   = '/home/gmasch/idl/savs/tb_bc/'
plot_path = '/home/gmasch/idl/plots/'

GET_TB_BC, station, date, start_d, stop_d, caldate  
jd = julday(strmid(date,2,2), strmid(date,4,2), '20' + strmid(date,0,2))

n_freq = 14
n_freq_wv = 7
n_freq_o2 = 7 
n_ang  = 6

cs = where(start_d ge 0. and caldate ge abscal)
n_cs = n_elements(cs)

dummy = -999.
TB_meas_mean = REPLICATE(dummy, n_freq, n_ang, n_cs)
TB_meas_std  = REPLICATE(dummy, n_freq, n_ang, n_cs)
TB_calc      = REPLICATE(dummy, n_freq, n_ang, n_cs)
TB_diff      = REPLICATE(dummy, n_freq, n_ang, n_cs)
ii           = REPLICATE(dummy, n_cs)

for i = 0, n_cs-1 do begin
  
;   ***COSMO profiles
    READ_COSMO_COLUMNS, station, date(cs(i)), time, z, T, p, q
    if z(0) eq dummy then goto, ABORT

    if station eq 'jue' then begin
      f = 0                ;time since model forecast
      c = [19, 20, 25, 26] ;averad columns
      
      z = average(z(*, c), 2)
      T = reform(average(T(*, f, *, c), 4))
      p = reform(average(p(*, f, *, c), 4))
      q = reform(average(q(*, f, *, c), 4))
      
    endif

     SEC2DATE1970, time, x, h_cosmo
         
     jj = where(h_cosmo ge start_d(cs(i)) and h_cosmo le stop_d(cs(i)))
     
     if jj(0) ne -1 then begin
	z  = reverse(z)
	T  = reverse(average(T(jj, *), 1))
	p  = reverse(average(p(jj, *), 1))
	q  = reverse(average(q(jj, *), 1))
     endif else goto, ABORT 
     
    data_dir = '/data/data_hatpro/'+station+'/data/level0c/'

   ;***find mwr data

   ;**actual date
    xx = data_dir + strmid(date(cs(i)), 0, 4) + '/' + date(cs(i)) + '*l0c.nc' 

   ;**next date
    GET_NEXT_DATE,  date(cs(i)), jahr_next, mtp_next, tag_next
    yy = data_dir + jahr_next+mtp_next+tag_next+'*l0c.nc'

    IF yy NE '' THEN xx = [yy, xx]
    ress = FILE_SEARCH(xx)
    n_ress = N_ELEMENTS(ress)

    IF ress(0) NE '' THEN BEGIN
        
    FOR i_ress = 0, n_ress-1 DO BEGIN
      READ_HATPRO_LEVEL0C_NC, ress(0), commment, lat, lon, alt, angles, freq, time, tb_x, flag_x, temp, pres, relh
      print,date(i)
      print, angles
      IF i_ress EQ 0 THEN BEGIN
      h_hatpro  = time
     
      tb = tb_x
      flag = flag_x
      ENDIF ELSE BEGIN
      h_hatpro = [h_hatpro + 24., time]
      tb = [tb, tb_x]
      flag = [flag, flag_x]
      ENDELSE
    ENDFOR

    i_comp = WHERE(h_hatpro ge start_d(cs(i)) and h_hatpro le stop_d(cs(i)) and flag ne 1)
   
   ENDIF ELSE BEGIN
      print, 'no comparable radiometer data found'
      GOTO, ABORT     
    ENDELSE
  
    IF N_ELEMENTS(i_comp) GE 2 THEN BEGIN
    FOR k = 0, n_freq-1 DO BEGIN
      FOR kk = 0, n_ang-1 DO BEGIN
      TB_meas_mean(k, kk, i) = MEAN(tb(k, kk, i_comp))
      TB_meas_std(k, kk, i)  = STDDEV(tb(k, kk, i_comp))
    
      ENDFOR
    ENDFOR
    ENDIF ELSE BEGIN 
    
       print, 'no comparable radiometer data found'
       GOTO, ABORT       
    ENDELSE

    jk = WHERE(FINITE(z) AND FINITE(T) AND FINITE(p) AND FINITE(q))
    IF N_ELEMENTS(jk) EQ -1 THEN GOTO, ABORT  

    IF N_ELEMENTS(angles) EQ 0 THEN GOTO, ABORT 

   ;***calculate RT 
    tc = 1
    FOR kk = 0, n_ang-1 DO BEGIN    
	IF kk GT 0 THEN tc = 0
	
	ABSHUM_TO_MIXR, T, p, q/1000., m
	
	STP, z(jk), T(jk), p(jk), q(jk)/1000., m(jk), -999., -999., -999., -999., $
		  angles(kk)-90., freq, 'r98', 1, TB_calcx, tau, tau_wv, tau_o2, tau_liq, T_mr, linew_22='lil05', cont_corr='tur09' 
      
	icheck = WHERE(TB_meas_mean(*, kk, i) LT 0. OR TB_calcx LT 0.)
	IF icheck(0) EQ -1 THEN BEGIN
	  TB_calc(*, kk, i) = TB_calcx  
	  TB_diff(*, kk, i) = TB_meas_mean(*, kk, i) - TB_calcx
	  
	  ii(i) = i
	ENDIF ELSE GOTO, ABORT
    ENDFOR
   ABORT:
endfor

n_ii = total(n_elements(where(ii ne dummy)))
TB_meas_minus_calc = REPLICATE(dummy, n_ii, n_freq_o2, n_ang) 

;***plot results as time_series

device, /color, filename=plot_path + 'plot_cosmo_bias_corr'

;**K-Band
dummyy = LABEL_DATE(DATE_FORMAT=['%D-%M','%Y'])

FOR aa = 0, n_ang-1 DO BEGIN

 miny = 999.
 maxy = -999.
 IF aa EQ 0 THEN mean_tb_diff = FLTARR(7)

 FOR k = 0, 6 DO BEGIN
  IF ii(0) NE -1 THEN BEGIN
   IF MIN(TB_diff(k, aa, ii)) LT miny THEN miny = MIN(TB_diff(k, aa, ii))
   IF MAX(TB_diff(k, aa, ii)) GT maxy THEN maxy = MAX(TB_diff(k, aa, ii)) 
  ENDIF
  IF aa EQ 0 THEN mean_tb_diff(k) = MEAN(TB_diff(k, aa, ii))
 ENDFOR

 delta = (maxy-miny)/10.

;  jj = WHERE(TB_diff(0, 0, *) NE dummy AND TB_diff(7, 0, *) NE dummy AND ABS(TB_diff(0, 0, *)-mean_tb_diff(0)) LT 4. AND ABS(TB_diff(6, 0, *)-mean_tb_diff(6)) LT 3.) ;original
   jj = WHERE(TB_diff(0, 0, *) NE dummy AND TB_diff(7, 0, *) NE dummy) ; test
 PLOT, jd(jj), TB_diff(0, aa, jj), xtitle = 'Date', ytitle = 'TB_meas-TB_calc [K]',$ 
       title = 'MWR-RS comparison, '+station+', elev. '+STRING(angles(aa), format = '(f6.1)'),$
       yrange = [miny-delta, maxy+delta], ystyle = 1, xstyle = 1, /NODATA,$
       XTICKFORMAT='LABEL_DATE', XTICKUNITS = ['Time', 'Time'], XTICKS=6
 OPLOT, [0., 400.], [0., 0.], linestyle = 1

 FOR j = 0, 6 DO BEGIN
  OPLOT, jd(jj), TB_diff(j, aa, jj), color = j+1, psym = -4
 ENDFOR

;**V-Band
 miny = 999.
 maxy = -999.

 FOR k = 7, 13 DO BEGIN
  IF ii(0) NE -1 THEN BEGIN
   IF MIN(TB_diff(k, aa, ii)) LT miny THEN miny = MIN(TB_diff(k, aa, ii))
   IF MAX(TB_diff(k, aa, ii)) GT maxy THEN maxy = MAX(TB_diff(k, aa, ii)) 
  ENDIF
 ENDFOR

 delta = (maxy-miny)/10.
 PLOT, jd(jj), TB_diff(0, aa, jj), xtitle = 'Date', ytitle = 'TB_meas-TB_calc [K]',$
       title = 'MWR-RS comparison,'+station+', elev. '+STRING(angles(aa), format = '(f6.1)'),$
       yrange = [miny-delta, maxy+delta], ystyle = 1, xstyle = 1, /NODATA,$
       XTICKFORMAT='LABEL_DATE', XTICKS=6, XTICKUNITS = ['Time', 'Time']
 OPLOT, [0., 400.], [0., 0.], linestyle = 1

 FOR j = 7, 13 DO BEGIN
  OPLOT, jd(jj), TB_diff(j, aa, jj), color = j+1-6, psym = -4
 ENDFOR
ENDFOR ; angles loop aa

;***create & plot offsets specs 

maxy_tot = 10.
miny_tot = -10.
 
 FOR i = 0, n_ii-1 DO BEGIN
    delta_tot = (maxy_tot-miny_tot)/10. 
    FOR aa = 0, n_ang-1 DO BEGIN
      FOR j = n_freq_wv, n_freq-1 DO BEGIN
      TB_meas_minus_calc(*, j-n_freq_wv, aa) = MEAN(TB_diff(j, aa, ii))
      ENDFOR
      IF aa EQ 0 THEN BEGIN

      PLOT, freq(7:n_freq-1), TB_meas_minus_calc(i, *, aa), xtitle = 'Frequency [GHz]', ytitle = 'TB_meas-TB_calc',$
; 	      title = 'MWR-RS comparison, '+time_span+', '+station, yrange = [miny_tot-delta_tot, maxy_tot+delta_tot], ystyle = 1,$
              title = 'MWR-RS comparison, '+station, yrange = [miny_tot-delta_tot, maxy_tot+delta_tot], ystyle = 1,$
	      xrange = [freq(7)-1, freq(13)+1], xstyle = 1 
      ENDIF
      OPLOT, freq(7:n_freq-1), TB_meas_minus_calc(i, *, aa), color = aa+1, psym = -4   
      XYOUTS, freq(7)-0.9, maxy_tot-(delta_tot*FLOAT(aa)), STRING(angles(aa), format = '(f6.2)'), color = aa+1
      OPLOT, [0., 100.], [0., 0.], linestyle = 1
    ENDFOR 

ENDFOR

;***save bias-offsets (V-Band)
comment = 'TB_meas_minus_calc = FLTARR(n_crit_ind-1, n_freq, n_ang)'
comment2 = 'freq = 51.26, 52.28, 53.86, 54.94, 56.66, 57.3, 58.0'
comment3 = 'ang = 90.0 42.0 30.0 19.2 10.2 5.4'

SAVE, filename = outfile+'.sav', comment, comment2, comment3, TB_meas_minus_calc, freq, angles, TB_calc
device, /close

STOP
END