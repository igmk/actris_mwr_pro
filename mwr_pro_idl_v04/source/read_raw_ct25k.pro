; --------------------------------------------------------------------------------
PRO READ_RAW_CT25K, $
; --------------------------------------------------------------------------------
    filename, $ ; name of input file
    time,     $ ; time as julian day, DBLARR( numprof )
    height,   $ ; height in m above ground, FLTARR( nz )
    detection_status, $ ; number of detected cloud levels (0..3) or flag for vertical visibility (4) 
              $ ; or something else (5), or data missing or suspect (99 - the original ascii is '/') INTARR(numprof)  
    status_flag,$ ; OK ('0'), Warning ('W') or Alarm ('A') , BYTARR(numprof)
    first_cbh, second_cbh, third_cbh,$ ; height above gnd (m) of detected cloud layers, FLTARR(numprof)
    vertical_vis, max_bscat,$ ; height above gnd (m) of derived vertical visibility, FLTARR(numprof), NEW Nov 2007
    laser_pulse_energy,$ ; laser pulse energy (% of optimum ?), FLTARR(numprof)
    laser_temperature, $ ; laser temperature (°C), FLTARR(numprof)
    receiver_sensitivity, $ ; receiver sensitivity (% of optimum ?), FLTARR(numprof)
    window_contamination, $ ; window contamination (0..2500mV), FLTARR(numprof)
    tilt_angle,           $ ; tilt angle of the instrument (° against vertical) , FLTARR(numprof)
    background_light,     $ ; background light (0..2500mV), FLTARR(numprof)
    sum_backscatter,      $ ; integrated backscatter from instrument (units ?), FLTARR(numprof)
    bscat,                $ ; backscatter profiles in (1e-4/km/sr), FLTARR(nz,numprof)
    measurement_parameters,$; information about actual device setup (string[stat_len]), BYTARR(stat_len,numprof)
    status_string,        $ ; current alarms and warnings Coded as a hex string (string[stat_len]), BYTARR(stat_len,numprof)
    ceilometername=ceilometername, $ ; every data record starts with 'CT<unit_number><software_level><message_number><message_subclass>' the procedure returns the string of the first record, NEW oct 2007
    cloud_data_string=cloud_data_string, $ ; line 2 of the last data message = cloud data in the form <N><W|A> <cld_hgt[0]> <cld_hgt[1]> <cld_hgt[2]> <statusstring>
    ceilo_para_string=ceilo_para_string, $ ; line 3 of the last data message = instrument paramters in the form <scale_factor> <bksc_prof_res> <length_prof_in_samples> <laser_pulse_energy> <laser_temp.> <window_transm_est> <tilt_angle> <bkgnd_light> <meas_parameters> <sum_of_bksct>
    scale_factor=scale_factor, $ ; first parameter in the status line of the ceilometer (intarr(numprof))
    measurement_mode=measurement_mode, $ ; second parameter in the status line of the ceilometer (bytarr(numprof))
    verbose=verbose ;-)
; --------------------------------------------------------------------------------
;+
; $Id: read_raw_ct25k.pro,v 1.7 2012/01/19 10:01:04 jschween Exp $
; Abstract:
;   read vaisala ct25k ceilometer raw ascii data
; Authors:
;   Susanne Crewell, Jan Schween
; Date:
;   2006-02-04
; Dependencies:
;   
; Changes:
; 2014-10-27 - Ulrich Loehnert: /COMPRESS keyword added to OPENR statements in order to cope with gzipped files now unter /data/tr32
; 2012-01-19 - Jan Schween
;  data format changed ????
;  tail dat does not start with date-time
;  try to catch it (introduce iLine and read until iLine gt 1 ... )
;  does not work ... :-(
; 2010-06-08 - Jan Schween
;  * if no data are in file return 'nodata' i.e. set variables to novalue 
;    => if numprof eq 0 ... do not try to make arrays with zero elements
; 2009-06-18 - Jan Schween
;  * added optional output of scale_factor and measurement_mode 
;  * replaced the ultimate STOP command by a warning message if scale_factor <> 100
; 2009-06-05 - Jan Schween
;   * added (optional) output of cloud_data_string and ceilo_para_string
; 2009-04-05 Jan Schween / Bernhard Prospichal
;  comments on the accuracy of dz and tiltangle ...
; 2008-12-12 Jan Schween:
;  + variable dummy replaced by novalue
;  + instrument parameters with security check for non values like '///' (see eg. file 081202 10:32:29)
;  + function f_hex2int replaced by reads... because it is faster
; 2007-10-XX Jan Schween:
;   + interpretation of detection_status = 4 for cloud levels
;     vertical visibility  and  altitude of maximum backscatter 
;     are now stored seperately
;     -> additional output variables vertical_vis and max_bscat
;   + new function f_hex2int to convert ascii hex numbers
;     (more efficient)
;   + more efficient counting of datarecords
;   + safer identification of datarecords (-> start_char)
;   + check for '/' sign in detection status
;     -> 'Raw data input to algorithm missing or suspect'
;     -> detection_status set to 99 
;   + added optional variable  verbose 
;       verbose =  0 => no messages
;       verbose >= 1 => some messages
;       verbose >= 5 => more messages and backscatter image for every 15min
;       verbose > 10 => many messages, every profile as plot
;     (the image is really fancy - try number 5 !)
;   + comments...
;-
; --------------------------------------------------------------------------------


;----------------------------------------------------------------------
; init vars
;----------------------------------------------------------------------

; external functions should be declared
; FORWARD_FUNCTION f_hex2int

; if not present set verbose level
IF NOT KEYWORD_SET(verbose) THEN verbose = 0

; non valid value
novalue = -9999.

; counter for scale_factors NE 100
cnt_scale_fac_warning = 0

; define str as string
str   = ' '

; conversion factor feet to meter
m2ft = 0.30479

; length of the status strings 
stat_len = 8

dz  = 30.          ; 30 m range bins 
; see Manual page 50: 
; return signal is sampled every 100ns with c=299 792 458m/s => dz = 29.9792458m
; or page 9: 
; dz = 50ft => 2*15.2395m = 30.479m
; or page 37:
; dz = 100ft = 30m ...
; we use here the latter !
nz  = 16*16        ; sixteen lines with sixteen altitudes = 256 values
n   = 24*60*60/15  ; maximum number of time steps in file

n_bksct = 252 ; the upper four values 252..255 of the backscatter signal are allways zero

; --------------------------------------------------------------------------------
; start 
; --------------------------------------------------------------------------------
if verbose ge 1 then begin
  print,'READ_RAW_CT25K:'
  print, 'Working on ',filename
endif

; --------------------------------------------------------------------------------
; count data records
; --------------------------------------------------------------------------------
; every data record starts with <x>CTA1210<y>
; where <x> is ascii 2 (CTRL+B) and <y> is ascii 3 (CTRL+C)
; we count the lines where the first char is <x>
; the string CTA1210 contains information about software level - 
; see gathering of variable ceilometername  below 
start_char = 1

; reset counter
numprof = 0
; open file and count ...
OPENR, unit, filename, /GET_LUN, /COMPRESS
WHILE NOT EOF(unit) DO BEGIN
  readf, unit, str
  if byte(STRMID(str,0,1)) eq start_char then numprof = numprof + 1
ENDWHILE
CLOSE,unit
FREE_LUN,unit

; message 
IF verbose GE 1 THEN begin 
  print, 'Number of backscatter profiles to read:', numprof
  ; every 30-data-minutes some information about progress
  n_progress = 4*30
  if verbose ge 5 then begin
    SET_PLOT,'x' ; open x-window for plots
    DEVICE, RETAIN=2, /decomposed ; set the retain flags to keep window content
  endif
endif

; if there is data in the file ...
if numprof gt 0 then begin

  ; allocate arrays 
  detection_status      = INTARR(numprof)
  status_flag           = BYTARR(numprof)
  scale_factor          = INTARR(numprof)
  measurement_mode      = BYTARR(numprof)
  laser_pulse_energy    = FLTARR(numprof)
  laser_temperature     = FLTARR(numprof)
  receiver_sensitivity  = FLTARR(numprof)
  window_contamination  = FLTARR(numprof)
  tilt_angle            = FLTARR(numprof)
  background_light      = FLTARR(numprof)
  measurement_parameters= BYTARR(stat_len,numprof)
  status_string         = BYTARR(stat_len,numprof)
  sum_backscatter       = FLTARR(numprof)

  time    = DBLARR(numprof)
  height  = FLTARR(nz)
  cloud_b = FLTARR(3,numprof)
  vert_vis= FLTARR(2,numprof)
  bscat   = FLTARR(nz,numprof)

  ; --------------------------------------------------------------------------------
  ; read data
  ; --------------------------------------------------------------------------------

  ; open file
  OPENR,unit,filename,/GET_LUN, /COMPRESS
  ; line number
  iLine = 0L

  ; define str2 - necessary to keep the line before the record start identifier
  str2 = '' 

  IF verbose GE 1 THEN t0 = systime(/julian)

  ; go through the file
  FOR nn = 0 , numprof-1 DO BEGIN

    ; search start of record = search start_char
    ; str must have an element s(0)
    str = 'this string is empty ;-)'
    WHILE (NOT EOF(unit)) AND (byte(STRMID(str,0,1)) NE start_char) DO BEGIN
      str2 = str
      READF,unit,str
      iLine += 1
      ENDWHILE

    ; --------------------------------------------------------------------------------
    ; date and time of record
    ; --------------------------------------------------------------------------------
    ; one line before the record start marker is the time 
    ;  - that is in str2 
    ; it looks like:
    ; -2007-05-24 14:40:12
    ;  123456789 123456789  
    year  = STRMID( str2, 3, 2)
    month = STRMID( str2, 6, 2)
    day   = STRMID( str2, 9, 2)
    date  = year + month + day

    hour  = STRMID( str2, 12, 2)
    minut = STRMID( str2, 15, 2)
    sec   = STRMID( str2, 18, 2)

    if verbose ge 5 then print, 'try to convert timestr. "',str2,'"'

    ; convert date-time to floating point number
    time(nn) = JULDAY(FIX(month),FIX(day),2000+FIX(year),FIX(hour),FIX(minut),FIX(sec))

    ; extract here the name of the ceilometer - once !
    if nn eq 0 then begin
      ceilometername = strmid(str,1,strlen(str)-2)
      ; actually this string is build as 
      ; CT<unit_number><software_level><message_number><message_subclass>
      ; where 
      ; <unit_number>      can be 0..9,A..Z
      ; <software_level>   can be 00..99
      ; <message_number>   we are analysing her message number 2 
      ; <message_subclass> we are analysing her subclass 3
      ; message number and subclass could be used as a check whether 
      ; we are dealing with the right dataformat for this record 
      ; ... just in case for future problems ?
      msg_format = strmid(ceilometername,5,1)+'.'+strmid(ceilometername,6,1)
      if verbose ge 1 then print, 'ceilometername="',ceilometername,'" => message format=',msg_format
      if msg_format ne '2.3' then begin 
        print, 'WARNING: message format is',msg_format,' ... it should be 2.3 !'
      endif
      
    end


    ; message to screen 
    IF verbose ge 10 THEN print, 'i=',nn,' time=',day,'.',month,'.',year,' ',hour,':',minut,':',sec

    ; --------------------------------------------------------------------------------
    ; analysed cloud altitudes or vertical visibilities
    ; --------------------------------------------------------------------------------
    ; the second line of the data record looks like
    ; '30 01230 12340 23450 FEDCBA98<CR>'
    ; or
    ; '50 ///// ///// ///// 00000800<CR>'
    ;   123456789 123456789 123456789'

    ; read text line 
    READF,unit,str
    iLine += 1

    ; save data string
    cloud_data_string = str

    ; detection status
    ; information about the number of detected 
    ; cloud layers or detected fog - see manual
    ; a '/' indicates 'Raw data input to algorithm missing or suspect'
    if STRMID(str,0,1) eq '/' then detstat = 99 $
    else detstat = FIX(STRMID(str,0,1))

    detection_status(nn) = detstat

    ; the status flag indicates any Warning ('W') or any alarm ('A') or all ok ('0')
    ; converted to byte for easier handling 
    ; bytearrays are similar to character-arrays the content is interpreted
    ; only sometimes differently (print, '+', ...)
    status_flag(NN) = BYTE(STRMID(str,1,1))

    ; initialize cloud_b with 'novalue' 
    FOR ib = 0 , 2 DO cloud_b(ib,nn) = novalue
    ; init vert_vis with 'novalue' 
    FOR ib = 0 , 1 DO vert_vis(ib,nn) = novalue
    ; evaluate cloud layers
    IF 0 LT detstat AND detstat LT 4 THEN BEGIN 
      FOR ib = 0, detstat-1 DO BEGIN
        str2  = STRMID(str,3+ib*6,5)
        IF str2 NE '/////' THEN BEGIN
          cloud_b(ib,nn) = FLOAT(str2)*m2ft
        ENDIF ELSE BEGIN ; this should not happen ... 
          cloud_b(ib,nn) = novalue
        ENDELSE
      ENDFOR ; FOR ib
    ENDIF ELSE IF detstat EQ 4 THEN BEGIN
        ; in this case the first field in the string is vertical visibility
        str2  = STRMID(str,3,5)
        IF str2 NE '/////' THEN vert_vis(0,nn) = FLOAT(str2) * m2ft
        ; the second field is (according to the manual) maximum backscatter (?)
        ; according to a comment in the AMF-netcdf files "Altitude of highest signal"
        str2  = STRMID(str,9,5)
        IF str2 NE '/////' THEN vert_vis(1,nn) = FLOAT(str2) * m2ft
      ENDIF ; ELSE IF destat

    ; decode the status string into bytes
    ; the status string contains information about current alarms and warnings
    ; Coded as a hex string. Each character code 4 bits which are flags for the 
    ; status of the subsystem. See manual page 35.
    ; The conversion to byte does not change the binary content of the variable
    ; this must be done for handling in the netcdf routines
    FOR i = 0,7 DO status_string(i,nn) = BYTE(STRMID(str,21+i,1))

    ; if verbose is set print first info about first datarecord
    IF nn EQ 0 AND verbose GE 5 THEN BEGIN
      print, 'Detection status    =', detection_status(nn)
      print, 'Status flag         =', STRING(status_flag(nn))
      print, 'Status string       =', STRING(status_string(0:7,nn))
      ; measurement_parameters are here not yet read !
      ; print, 'meas.params         =', STRING(measurement_parameters(0:7,nn))
      ;print,nn,' ',day,hour,minut,sec,'',time(nn),cloud_b(*,nn)
    ENDIF

    ; --------------------------------------------------------------------------------
    ; status line
    ; --------------------------------------------------------------------------------
    ; the status line contains info about laser temperature, window contamination etc. 
    ; it looks like
    ; '100 N  53 +34 204  146  +2  621 LF7HN1 139<CR>'
    ;   123456789 123456789 123456789 123456789 12
    ; coding for:
    ; <scale_factor> <measurement_mode> <laser_pulse_energy> <laser_temp.> ...
    ;   ...<receiver_sens.> <window_cont.> <tilt_angle> <bkgnd_light> <meas_parameters> <sum_of_bksct>

    ; read status line
    READF,unit,str
    iLine += 1

    ; save status string
    ceilo_para_string = str

    IF nn EQ 0 AND verbose GE 5 THEN print,'read status line:"',str,'"'

    ; decode elements
    scale_factor[nn] = FIX(STRMID(str,0,3))
    IF nn EQ 0 AND verbose GE 5 THEN print, 'scale_factor        =', scale_factor[nn]

    measurement_mode(nn)  = BYTE(STRMID(str,4,1))
    IF nn EQ 0 AND verbose GE 5 THEN print, 'measurement mode    =', measurement_mode

    ; laser pulse energy
    x = STRMID(str,5,4)
    IF STRPOS( x , '/' ) NE -1 THEN x = novalue
    IF x GE 10. AND x LE 110. THEN BEGIN
      laser_pulse_energy(nn) = x
    ENDIF ELSE laser_pulse_energy(nn) = novalue
    IF nn EQ 0 AND verbose GE 5 THEN print, 'laser pules energy  =', laser_pulse_energy(nn)

    ; laser temperature
    x = STRMID(str,9,4)
    IF STRPOS( x , '/' ) NE -1 THEN x = novalue
    IF x GE 10. AND x LE 110. THEN BEGIN
      laser_temperature(nn) = x
    ENDIF ELSE laser_temperature(nn) = novalue
    IF nn EQ 0 AND verbose GE 5 THEN print, 'laser temperature   =', laser_temperature(nn)

    ; receiver sensitivity
    x = STRMID(str,13,4)
    IF STRPOS( x , '/' ) NE -1 THEN x = novalue
    IF x GE 0. AND x LE 999. THEN BEGIN
      receiver_sensitivity(nn) = x
    ENDIF ELSE receiver_sensitivity(nn) = novalue
    IF nn EQ 0 AND verbose GE 5 THEN print, 'receiver sensitivity=', receiver_sensitivity(nn)

    ; window contamination
    x = STRMID(str,17,5)
    IF STRPOS( x , '/' ) NE -1 THEN x = novalue
    IF x GT 0. AND x LT 2500. THEN BEGIN
      window_contamination(nn) = x
    ENDIF ELSE window_contamination(nn) = novalue
    IF nn EQ 0 AND verbose GE 5 THEN print, 'window contamination=', window_contamination(nn)

    ; tilt angle
    x = STRMID(str,22,4)
    IF STRPOS( x , '/' ) NE -1 THEN x = novalue
    IF x GE -15. AND x LE 90. THEN BEGIN
      tilt_angle(nn) = x
    ENDIF ELSE tilt_angle(nn) = novalue
    IF nn EQ 0 AND verbose GE 5 THEN print, 'tilt angle          =', tilt_angle(nn)

    ; background light
    x = STRMID(str,26,5)
    IF STRPOS( x , '/' ) NE -1 THEN x = novalue
    IF x GE 0. AND x LE 2500. THEN BEGIN
      background_light(nn) = x
    ENDIF ELSE background_light(nn) = novalue
    IF nn EQ 0 AND verbose GE 5 THEN print, 'background light    =', background_light(nn)

    ; measurement parameters
    ; information about actual device setup
    FOR i=0,5 DO measurement_parameters(i,nn)= BYTE(STRMID(str,32+i,1))
    IF nn EQ 0 AND verbose GE 5 THEN print, 'meas.parameters     =', string(measurement_parameters(*,nn))

    ; sum of backscatter signal
    x = STRMID(str,38,4)
    IF STRPOS( x , '/' ) NE -1 THEN x = novalue
    IF x GE 0. AND x LE 999. THEN BEGIN
      sum_backscatter(nn) = x
    ENDIF ELSE sum_backscatter(nn) = novalue
    IF nn EQ 0 AND verbose GE 5 THEN print, 'sum backscatter     =', sum_backscatter(nn)

    ; Error message and program STOP if scale_factor is not 100
    ; according to the manual the scale_factor can be set to any value between 0 and 9999
    ; but this is usually not done by us ...
    IF scale_factor[nn] NE 100 THEN BEGIN
      ; STOP,' Scale factor is not set to 100 as should be; programm interrupted'
      if cnt_scale_fac_warning LE 3 then message, ' Scale factor is NOT 100 as usually !', /informational
      cnt_scale_fac_warning += 1
    ENDIF ELSE IF (scale_factor[nn] EQ 100) AND (cnt_scale_fac_warning ne 0) THEN BEGIN
      message, string('There were ',cnt_scale_fac_warning,' scale_factors NE 100 !'), /informational
      cnt_scale_fac_warning = 0
    ENDIF

    ; --------------------------------------------------------------------------------
    ; profile data
    ; --------------------------------------------------------------------------------
    ; the profile data comes in 16 lines each with values for 16 layers.
    ; the first three digits are the height index, it follows an ascii string 
    ; with 16*4characters with the values as Hexstring for 16bit signed integers.
    ; The range of the instrument is 25,000 ft with a vertical resolution 
    ; of 100ft => 251 values => the last four values in the datastring are 
    ; allways set to zero 
    ; the datastring looks like:
    ; '00000100010000E000E000E000D000E000F000E000C000D000B000D000F000D000E'
    ;   123456789 123456789 123456789 123456789 123456789 123456789 123456

    ; reset height index
    iz = 0
    ; define bxxx as signed 2 byte integer
    bxxx = 0S
    ; go through lines 
    FOR j = 0,15 DO BEGIN
      READF,unit,str
      iLine += 1
      ; first 3digits are the index of the first element
      alt = FIX(STRMID(str,0,3))
      ; go through columns 
      FOR jj=0,15 DO BEGIN
        height(iz) = dz/2. + (alt + jj)*dz
        ; convert 4 character hex number - takes about 4ms/ceilometerrecord - too slow
        ; bscat(iz,nn) = f_hex2int( STRMID(str,3+(jj*4),4) )
        ; alternative with reads - faster by more than a facotr two
        READS, STRMID(str,3+(jj*4),4), bxxx, FORMAT='(Z)'
        bscat(iz,nn) = bxxx
        ; increment height index
        iz = iz + 1
      ENDFOR ; FOR jj
    ENDFOR ; FOR j



    ; show progress graphically
    IF (verbose ge 1) && ((nn mod n_progress) eq 0) then begin
      ; a simple progress bar made of '*' on the terminal
      print,'*',format='(A,$)'
      IF (verbose ge 5) && (verbose lt 10) then begin
        ; this is a nice game: 
        ; display currently porcessed backscater profiles
        ; and scroll them through the display ...
        ; select only the part of the data that fits into the screen 
        img = transpose( bscat( * , max([0,nn-!D.x_size]):nn ) )
        ; erase values smaller than 1.0 (we are using log and dont want to have negative values)
        idx = where( img lt 1 , cnt )
        if cnt gt 0 then img( idx ) = 1.
        ; log scaling
        img = 255 * alog(img) / alog(600)
        ; display image 
        tv, img
      endif ; verbose ge 5
    endif ;  verbose ge 1
    ; for very high verbose values every profile is shown
    IF verbose GE 10 THEN PLOT, bscat(*,nn), height, xrange=[-5,200], xstyle=1

  ENDFOR ; FOR nn ...

  IF verbose GE 1 THEN BEGIN
    t1 = systime(/julian)
    print,'--------------------------------------------------'
    print,(t1-t0)*24*60*60/nn,format='("this took ",F,"sec per record")'
    print,'--------------------------------------------------'
  ENDIF


  ; close file
  CLOSE, unit
  FREE_LUN, unit

  ; new line after simple progress bar
  if verbose ge 1 then print

  ; correct height for tilting angle of the instrument
  ; this is a bit ambigous:
  ; 1st: the variable in the AMF netcdf files is 'range'
  ; - i.e. the distance of signal to the instrument
  ; but not the height above ground
  ; 2nd: the tilt angle might vary in one file 
  ; we consider here only the last value ...
  ; 3rd the tilt angle might change from time to time by 1degree (inaccuracy, themral effects ...)
  ; this results for tilt angles around 20deg in errors around 0.6% 
  if tilt_angle(numprof-1) GT novalue THEN BEGIN
    height = height * COS(tilt_angle(numprof-1) * !DTOR ) ; !DTOR = 'degrees to radians' = pi/180
    ENDIF

  ; reformat the arrays to the actual number of values
  ; or copy them in the appropriate vars (cloud_b and vert_vis)
  detection_status      = detection_status(0:nn-1)
  status_flag           = status_flag(0:nn-1)
  scale_factor          = scale_factor(0:nn-1)
  measurement_mode      = measurement_mode(0:nn-1)
  laser_pulse_energy    = laser_pulse_energy(0:nn-1)
  laser_temperature     = laser_temperature(0:nn-1)
  receiver_sensitivity  = receiver_sensitivity(0:nn-1)
  window_contamination  = window_contamination(0:nn-1)
  tilt_angle            = tilt_angle(0:nn-1)
  background_light      = background_light(0:nn-1)
  measurement_parameters= measurement_parameters(*,0:nn-1)
  status_string         = status_string(*,0:nn-1)
  sum_backscatter       = sum_backscatter(0:nn-1)
  time                  = time(0:nn-1)
  height                = height(0:n_bksct-1) ; the upper four values of the backscatter signal are allways zero.
  first_cbh             = REFORM(cloud_b(0,0:nn-1))
  second_cbh            = REFORM(cloud_b(1,0:nn-1))
  third_cbh             = REFORM(cloud_b(2,0:nn-1))
  vertical_vis          = REFORM(vert_vis(0,0:nn-1))
  max_bscat             = REFORM(vert_vis(1,0:nn-1))
  bscat                 = bscat(0:n_bksct-1,0:nn-1)


endif else begin ; no data in file

  ; delete all variables
  detection_status = novalue
  status_flag = novalue
  scale_factor          = novalue
  laser_pulse_energy = novalue
  laser_temperature = novalue
  receiver_sensitivity = novalue
  window_contamination  = novalue
  window_transm_est = novalue
  tilt_angle = novalue
  background_light = novalue
  measurement_parameters = novalue
  status_string = novalue
  sum_backscatter = novalue

  time = novalue
  height = novalue
  first_cbh = novalue
  second_cbh = novalue
  third_cbh = novalue
  vertical_vis = novalue
  max_bscat = novalue
  bscat = novalue

endelse ; no data in file



END

