;------------------------------------------------------------
PRO GET_LOG_MOD, $
   record,  $ ; INPUT: flag for talkative mode: the higher the value the more information
   filename,$ ; INPUT: name of file to be read
   n_f,     $ ; OUTPUT: number of frequencies
   freq,    $ ; OUTPUT: frequencies of instrument; FLTARR(n_f)
   cal_type,$ ; OUTPUT: type of calibration: 0=gain, 1=noise, 2,3=tip curve results, 2=.. with short and 3=.. with long information; DBLARR(n)
   cal_time,$ ; OUTPUT: time of calibration in seconds since 1.1.2000 00:00:00; DBLARR(n)
   tip_curve_stat,$ ; OUTPUT: status of tip curve calibration: 2=success, 3=failed (only set if cal_type=2 or 3, otherwise 9d-33) DBLARR(n) 
   gain,    $ ; OUTPUT: slope of the linear Response per channel; FLTARR(n,n_f)
   tsys,    $ ; OUTPUT: system noise temperature, offset of the linear Response per channel; FLTARR(n,n_f)
   chisqu,  $ ; OUTPUT: mean normalized square of deviation from the sky tipping curve per channel; FLTARR(n,n_f) 
   lincorr, $ ; OUTPUT: correlation coefficient of the sky tipping curve per channel ; FLTARR(n,n_f)
   tnoise,  $ ; OUTPUT: noise source temperature per channel; FLTARR(n,n_f)
   n_ang,   $ ; OUTPUT: number of angles per tip curve FLTARR(n)
   airmass, $ ; OUTPUT: airmass per tipping curve and angle; FLTARR(n, max_ang)
   volts,   $ ; OUTPUT: voltages of tipping curve,last entry is voltage measured on hot load FLTARR(n, max_ang+1)
   success, $ ; OUTPUT: success flag of tipping curve cal. per channel: 2=success, 3=failed FLTARR(n, n_r1)
   tau,     $ ; OUTPUT: optical thicknesses per tip curve cal. and channel and angle FLTARR(n, n_r1) 
   constant,$ ; OUTPUT: parameter A of linear fit to tip curve: TB = A+B*tau ; FLTARR(n, n_r1)
   coeff,   $ ; OUTPUT: parameter B of linear fit to tip curve: TB = A+B*tau ; FLTARR(n, n_r1)
   success_test=success_test, $ ; INPUT, OPTIONAL: flag to check for TauSuccess(0). According to the manual ....
   no_lincorr_test=no_lincorr_test, $ ; INPUT, OPTIONAL: flag to check NOT for lincorr[1]=-1 or lincorr[2]=-1
   n_warnings=n_warnings ; OUTPUT, OPTIONAL: number of warnings occured
;..............................................................
;+
; $Id: get_log_mod.pro,v 1.1 2010/03/09 13:45:18 jschween Exp $
; Abstract:
;   SC 07/02/20
;     programm to read CAL.LOG files produced by HATPRO or DPR
;
; INPUT:
; record   flag for talkative mode
; filename name of log file (program works for linux and windows)
;
; OUTPUT:
; for software versions starting around 5.30 (mid 2006)  an
; additional calibration type (3) is possible which includes
; details on tipping curve performance
; n_ang     number of angles used in tipping curve (max 10!)
; airmass   airmass factor for the different angles
; volts     voltage for each frequency and angle + hot load
; success   indicates success for different frequencies
; tau       opacity as function of frequency and angle
; constant, coeff   tau = constant + coeff*airmass
;
; CHANGES
; JS 9.3.2010
;  * check for number of calibrations in header and calibration found in the file...
;  => warning message and reformat of files if inconsistent
; JS 2.6.2008
; * inconsistency with log file description found:
;   - the old version (SC 07/02/20) reads the detailed tip curve calibration 
;     information in every case
;   - the manual (Version 7.40 April 2008 p.130) says that this information is
;       only given if TauSuccess(0)=1
;   - according to Uli there has been at some time a change in the format ...
;   + included an additional if block which tests for TauSuccess(0)=1
;     this test must be switched on by setting the /success_test flag:
;       GET_LOG_MOD, ... ,/success_test
;    if not given the program works in the old way,
;    but a warning message is given if NOT TauSuccess(0) = 1 and the data is read.
; * There is a PROBLEM with the GOTO,SKIP_160 command:
;   if lincorr[1]=-1 or lincorr[2]=-1 variable tau_flag is set to 1 or 2 and 
;   reading of the sky tip information is interupted by GOTO,SKIP_160
;   if this happens the file pointer is positioned wrong and further reading gives
;   wrong values
;   => optinal parameter /no_lincorr_test to switch off the test for lincorr
;   with
;     GET_LOG_MOD, ... ,/no_lincorr_test 
;   if no_lincorr_test is not set and lincorr=-1 is encountered a WARNING is
;   printed and a jump to skip_160 is performed
; * comments, comments ....
; * Error message and exit if file identifier is wrong
; * warnings if read past end of file (EOF) or if more records are found 
;     in file than inidcated in the header file header
; * warnings of undocumented values of calibration type and cal_stat are found
; * warnings if the number of read calibrations is not identical to the 
;     number given in the file header
; * included the check for i_c in the WHILE condition of the MAIN reading loop
;   => replaces the 'GOTO, LEAVE' command
; * the talkative mode has now four levels: 
;    0 => quite, 
;    1 => messages at start and statisitcs at the end
;    5,10 => more blablabla ...
; * n_warnings :optional variable to count the warnings
;-
;------------------------------------------------------------

; if not present define and 'unset' the success_test flag
if ~keyword_set(success_test) then success_test=0

; if not present define and 'unset' the no_lincorr_test
if ~keyword_set(no_lincorr_test) then no_lincorr_test=0

; talkative messages
if record ge 1 then BEGIN
  print,'GET_LOG_MOD:'
  print,'Reading calib file "',filename,'"'
  print, 'success-test-flag=',success_test
ENDIF

n_warnings = 0

; open file
OPENR,unit,filename,/GET_LUN

; read Header information

; read file identifier - must be 657643
code=1l
READU,unit,code
; if the file idnetifier is wrong drop a message and break
IF code NE 657643 THEN MESSAGE, 'ERROR: "'+filename+'" is not a hatpro LOG file !'

n_gain=0l
READU,unit,n_gain
IF record ge 1 THEN print,'Number of recorded gain samples =',n_gain

n_noise=0l
READU,unit,n_noise
IF record ge 1 THEN print,'Number of recorded noise calibration samples =',n_noise

n_skydip=0l
READU,unit,n_skydip
IF record ge 1 THEN print,'Number of recorded tip curve samples =',n_skydip

n_r1=0l
READU,unit,n_r1
IF record ge 1 THEN print,'Number of receiver 1 channels =',n_r1

n_r2=0l
READU,unit,n_r2
IF record ge 1 THEN print,'Number of receiver 2 channels =',n_r2

n_f = n_r1 + n_r2
freq = FLTARR(n_f)
READU,unit,freq
IF record ge 1 THEN BEGIN
  print,'Number of frequencies =', n_f
  print,'Frequencies=',freq
ENDIF

; total number of calibrations
n = n_gain + n_noise + n_skydip

IF record ge 1 THEN BEGIN
  print,'n = n_gain + n_noise + n_skydip = ',n_gain,' + ',n_noise,' + ',n_skydip,' = ',n
ENDIF

IF record ge 5 THEN BEGIN
  s = ''
  Print,'press <return> to continue'
  READ,s
ENDIF


; init data fields as double, filled with a 'novalue' of type double
; question:
;   cal_type and tip_curve_stat in the CAL.LOG files are just integers 
;   with maxima four possible values - why are they predefined here as double precision fields ?
;   JS 2.6.2008
cal_type       = REPLICATE(-9d+33, n)
cal_time       = Replicate(-9d+33, n)
tip_curve_stat = Replicate(-9d+33, n)

; init other fields as float array 
gain    = FLTARR(n,n_f)
Tsys    = FLTARR(n,n_f)
chisqu  = FLTARR(n,n_f)
lincorr = FLTARR(n,n_f)
Tnoise  = FLTARR(n,n_f)
; preset them with a single precision 'novalue'
gain[*,*]    = -9e33
Tsys[*,*]    = -9e33
chisqu[*,*]  = -9e33
lincorr[*,*] = -9e33
Tnoise[*,*]  = -9e33

;........................following parameters only if cal type = 3
; maximum number of angles for tip curve 
max_ang = 20

; init fields for (succesful) sky tip calibrations ...
n_ang   = FLTARR(n)                   ; number of angles per tip curve
airmass = FLTARR(n, max_ang)          ; maximum number of 10 angles is assumed
volts   = FLTARR(n, n_r1, max_ang+1)  ; last entry is voltage measured on hot load
success = FLTARR(n, n_r1)             ; success flag per channel: 2=success, 3=failed
tau     = FLTARR(n, n_r1, max_ang)    ; thicknesses per per channel and tip curve
constant= FLTARR(n, n_r1)             ; parameter A of linear fit to tip curve: TB = A+B*tau
coeff   = FLTARR(n, n_r1)             ; parameter B of linear fit to tip curve: TB = A+B*tau

; an unknown flag (JS 2.6.2008)
; it has something to do with a check of the sky tip calibs...
tau_flag = 0

; counter for calibrations
i_c = 0l
n_gain_c   = 0l
n_noise_c  = 0l
n_skydip_c = 0l
n_skydipX_c= 0l
n_skytip_succes_c = 0l

; last time 
time = 0.d

; --------------------------------------------------------------------------------
; MAIN reading loop ...
; --------------------------------------------------------------------------------
WHILE NOT EOF(unit) AND (i_c LT n) DO BEGIN

  ; read calibration type
  y = 0l
  READU,unit,y
  cal_type(i_c) = y

  ; save last time
  t_last = time
  ; read time of calibration - as allways with hatpro: date/time in seconds since 1.1.2001
  READU,unit,y
  cal_time(i_c) = y
  ; convert to julian day
  time = cal_time(i_c)/86400.d + JULDAY(1,1,2001,0,0,0)
  ; convert time to string ...
  time_str = string( time , format='(C(CYI04,"/",CMOI02,"/",CDI02,"/",CHI02,":",CMI02,":",CSI02))' )

  if record ge 5 then print,'t[',i_c,']=',time,' = ',time_str 

; old conversion
; replaced by JS June 2008
;   CALDAT, h, M1, D1, Y1, h1, mi1, s1
;   time_str = STRING(y1,FORMAT='(I4)')+'/'+ STRING(m1,FORMAT='(I2)')+'/'+$
;              STRING(d1,FORMAT='(I2)')+' '+ STRING(h1,FORMAT='(I2)')+':'+$
;              STRING(mi1,FORMAT='(I2)')+':'+STRING(s1,FORMAT='(I2)')

  IF t_last GT time THEN BEGIN
    MESSAGE,'WARNING: inconsistency in time at i='+string(i_c,format='(I0)')+' last_time='+string(t_last,format='(F0)')+' > time='+string(time,format='(F0)'),/Informational
    ; MESSAGE,time_str+' > '+string( t_last , format='(C(CYI04,"/",CMOI02,"/",CDI02,"/",CHI02,":",CMI02,":",CSI02))' ),/Informational
    n_warnings += 1
  ENDIF

  ; verbose message
  IF record ge 5 THEN print,'Cal_type=',cal_type(i_c),' at ',time_str

  ; count calibration types
  CASE cal_type(i_c) OF
    0: N_gain_c += 1
    1: N_noise_c += 1
    2: N_skydip_c += 1
    3: N_skydipX_c += 1
    ELSE: BEGIN 
      MESSAGE,'WARNING: unknown calibration type found !',/Informational
      print,'time=',time_str,' calibration type =',cal_type(i_c)
      n_warnings += 1
      ENDELSE
  ENDCASE

  ; if calibration is tip curve (with short (cal_type=2) or long information (3) ...
  IF cal_type(i_c) EQ 2 OR cal_type(i_c) EQ 3 THEN BEGIN
    ; .... read status of tip curve calibration
    y = 0l
    READU,unit,y
    tip_curve_stat(i_c) = y
    ; evaluate this flag:
    CASE tip_curve_stat(i_c) OF
      3: BEGIN ; stat=3 will be more frequent - check for it first
         ; create Message
         s = ' -> Failed'
         END
      2: BEGIN
         ; count succesful calibrations
         n_skytip_succes_c += 1
         ; create message
         s = ' -> Success'
         END
      ELSE: BEGIN
        MESSAGE,'WARNING: tip_curve_stat has invalid value'+string(tip_curve_stat(i_c))+' at '+time_str,/INFORMATIONAL
        n_warnings += 1
        END
    ENDCASE
    ; print message
    IF record ge 5 THEN print,'Tip curve result at ',time_str, ' :', tip_curve_stat(i_c), s
  ENDIF ; cal_type EQ 2 ...

  ; read calibration gains
  x = FLTARR(n_f)
  READU,unit,x
  gain(i_c,*) = x
  IF record ge 5 THEN print,  'Gain               =',gain(i_c,0),i_c,n,cal_type(i_c),n_gain,n_noise,n_skydip

  ; if calibration is not a gain calibration ... 
  IF cal_type(i_c) NE 0 THEN BEGIN
    ; ... read system noise temperaure of all channels
    x = FLTARR(n_f)
    READU,unit,x
    tsys(i_c,*) = x
    IF record ge 5 THEN print,  'Tsys               =',tsys(i_c,0), ' K'
  ENDIF

  ; if calibration is a sky tip curve ... 
  IF cal_type(i_c) EQ 2 OR cal_type(i_c) EQ 3 THEN BEGIN
    ; ... read linear correlation (coefficients ?) 
    x = FLTARR(n_f)
    READU,unit,x
    lincorr(i_c,*) = x
    IF record ge 5 THEN print,'Linear Correlation =',lincorr(i_c,0)
    ; set tau_flag if correlation of channels 1 or 2 has flagged value -1 i.e. no fit was possible
    IF NOT no_lincorr_test THEN BEGIN
      IF lincorr(i_c,2) EQ -1.0 THEN tau_flag = 2
      IF lincorr(i_c,1) EQ -1.0 THEN tau_flag = 1
      IF tau_flag NE 0 THEN BEGIN
        MESSAGE,'WARNING: encountered a lincorr['+string(tau_flag,format='(I1)')+']=-1 : will exit reading sky tip information - this may cause problems!',/INFORMATIONAL
        n_warnings += 1
      ENDIF
    ENDIF ELSE BEGIN ; no_lincorr_test
      IF (lincorr(i_c,1) EQ -1.0) OR (lincorr(i_c,2) EQ -1.0) THEN BEGIN
        MESSAGE,'WARNING: encountered a lincorr=-1 : will ignore this',/INFORMATIONAL 
      ENDIF 
    ENDELSE 
    ; read chi-square values of all channels
    x = FLTARR(n_f)
    READU,unit,x
    chisqu(i_c,*) = x
    IF record ge 5 THEN print,'Chi Square         =',chisqu(i_c,0)
    ; read noise temperatures of all channels
    x = FLTARR(n_f)
    READU,unit,x
    Tnoise(i_c,*) = x
    ; message
    IF record ge 5 THEN print,'Tnoise             =',tnoise(i_c,0)
  ENDIF

  ; if calibration contains long information about the sky tip ...
  IF cal_type(i_c) EQ 3 THEN BEGIN
    ; ... read number of angles
    y = 0l
    READU,unit,y
    n_ang(i_c) = FIX(y)
    IF record ge 5 THEN print,'Number of angles =',n_ang(i_c)
    IF n_ang(i_c) GT max_ang THEN Print,'Too many angles!',n_ang(i_c)
    ; read airmasses corresponding to angles 
    yy = FLTARR(n_ang(i_c))
    READU,unit,yy
    airmass(i_c,0:n_ang(i_c)-1) = yy
    IF record ge 5 THEN print,'Airmasses            =',airmass(i_c,*)
    ; read detector voltages correspondng to angles/airmasses + voltage of hot target
    xx = FLTARR(n_r1,n_ang(i_c)+1)
    READU,unit,xx
    volts(i_c,*,0:n_ang(i_c)) = xx
    IF record ge 5 THEN print,'Volts@f1             =',volts(i_c,*,0)
    ; read succes flag - name in documetaion is 'TauSuccess'
    h = REPLICATE(0l, n_r1)
    READU,unit,h
    Success(i_c,*) = h
    IF record ge 5 THEN print,'Success (0=no, 1=yes for tau calc., 2=yes for tau calc. and thresholds met)=',success(i_c,*)


  IF (Success(i_c,0) GT 0) OR (NOT success_test) THEN BEGIN
; !!!!!!!!!!!!!!!!!!!!!!!!!! 
; the manual says that the following information comes 
; only if Success(i_c,0) = TauSuccess(first channel) = 1
; the version of Susanne from Feb 2006 ignores this
; may be that is the reason for the "if NOT EOF..." security things.
; According to Uli there has been at some time a change in the format ...
; and some files have errors (which i dont believe)
; the if-then block here tests on success
; must be switched on with the /success_test flag:
;   GET_LOG_MOD, ... , /success_test
; if not switched on the programm behaves as before...
; JS 2.6.2008
; !!!!!!!!!!!!!!!!!!!!!!!!!!

    ; Warning if we are doing things not in accordance with the manual 
    IF Success(i_c,0) EQ 0 THEN BEGIN
      MESSAGE,'WARNING: reading sky tip information although calibration was not succesful',/INFORMATIONAL
      MESSAGE,'According to the manual (v.7.40 April 2008) this information should not appear in the file',/INFORMATIONAL
      n_warnings += 1
    ENDIF 

    ; read information of sky tip for all channels
    FOR i_f = 0, n_r1-1 DO BEGIN
     ; ?????????????????????????????
     ; ... why check here for EOF ?
     ; JS June 2008
     ; ?????????????????????????????
     IF NOT EOF(unit) THEN BEGIN
       ; ... read tau array for channel i_f
       yy = FLTARR(n_ang(i_c))
       READU,unit,yy
       tau(i_c,i_f,0:n_ang(i_c)-1) = yy
       IF record ge 5 THEN print,'Opacity=',tau(i_c,i_f,0:n_ang(i_c)-1)
       ; read parameter A (offset) for all channels of sky tip fit (Volt = A + B*tau )
       xx = 0.
       READU,unit,xx
       constant(i_c,i_f) = xx
       ; read paramter B (slope) for all channels 
       xx = 0.
       READU,unit,xx
       coeff(i_c, i_f) = xx
       IF record ge 5 THEN print,'Fit=',constant(i_c,i_f),' + ',coeff(i_c,i_f),'*airmass'

      ; in case no fit to the opacity values was possible (Lincor =-1) no coefficients are stored
       If tau_flag EQ 1 AND i_f EQ 0 THEN BEGIN
         MESSAGE, 'WARNING: no fit to the tau values was possible !',/INFORMATIONAL
         n_warnings += 1
         ; actually i dont like this - but i dont touch it.
         ; JS June 2008
         GOTO,SKIP_160 
       ENDIF ; if tau_flag
       If tau_flag EQ 2 AND i_f EQ 1 THEN BEGIN
         MESSAGE, 'WARNING: no fit to the tau values was possible !',/INFORMATIONAL
         n_warnings += 1
         ; actually i dont like this - but i dont touch it.
         ; JS June 2008
         GOTO,SKIP_160
       ENDIF ; if tau_flag

      ENDIF $ ; not EOF
     ELSE BEGIN ; EOF
      MESSAGE, 'WARNING : read past EOF while reading sky tip data !',/INFORMATIONAL
      n_warnings += 1
     ENDELSE ; EOF

    ENDFOR ; i_f

  ENDIF ; success gt 0

  ENDIF ; cal_type = 3 => log sky tip info

  ; jump target for erroneous paramter reading...
  SKIP_160:
  ; reset flag ...
  tau_flag = 0
  ; index for next calibration
  i_c = i_c + 1l

ENDWHILE ; NOT(EOF) - Main file reading loop



; reformat fields
if i_c ne N then begin
  MESSAGE,'inconsistency in number of calibrations:', /informational
  print,'        ','n_gain',' + ','n_noise',' + ','n_sky',' + ','n_skyx',' = ',N,format='(5(A,A6))'
  print,' Header ',n_gain  ,' + ',n_noise  ,' + ',n_skydip  ,' + ',  0       ,' = ',N,format='(5(A,I6))'
  print,' found  ',n_gain_c,' + ',n_noise_c,' + ',n_skydip_c,' + ',n_skydipX_c,' = ',i_c,format='(5(A,I6))'
  print,' field lengths will be corrected !'

  N = i_c

  cal_time = cal_time[0:N-1]
  cal_type = cal_type[0:N-1]

  tip_curve_stat = tip_curve_stat[0:N-1]

  gain    = gain[   0:N-1,*]
  Tsys    = Tsys[   0:N-1,*]
  chisqu  = chisqu[ 0:N-1,*]
  lincorr = lincorr[0:N-1,*]
  Tnoise  = Tnoise[ 0:N-1,*]

endif ; i_c lt n


; check times
time = cal_time/86400.d + JULDAY(1,1,2001,0,0,0)
N_time = N_elements(time)
t0 = time[0]
tN = time[N_time-1]
dt = tN - t0

; some messages about statisitcs
IF record ge 1 THEN BEGIN
  print,'Calibration Times:'
  print,'N = ',N_time
  IF (t0 LT -1095) OR (t0 GT 1827933925) OR (tN LT -1095) OR (tN GT 1827933925) THEN BEGIN
    print, 'INVALID Jullian date(s) ! t[0]=',t0  ,' t[N-1]=',tN
  ENDIF ELSE BEGIN
    print, time[0]  ,time[N_time-1], format='("  0, N-1:",2(C(CDI02,".",CMOI02,".",CYI04," ",CHI02,":",CMI02," ")) )'
    print, min(time),max(time),where(time eq min(time)),where(time eq max(time)), format='("min, max:",2(C(CDI02,".",CMOI02,".",CYI04," ",CHI02,":",CMI02," ")),I,I )'
    print, dt/N_time*24*60, dt/N_time-0.5, format='("i.e. in average every ",F8.3," min =",C(CHI02,":",CMI02,":",CSI02," hh:mm:ss"))'
  ENDELSE

  print, 'Found following calibrations:'

  print, n_gain_c,   ' gain ',format='(I5,A,$)'
  if n_gain_c gt 0 then print, dt/n_gain_c-0.5,format='(" i.e. every ",C(CHI02,":",CMI02,":",CSI02)," hh:mm:ss")' else print,''

  print, n_noise_c,   ' noise ',format='(I5,A,$)'
  if n_noise_c gt 0 then print, dt/n_noise_c-0.5,format='(" i.e. every ",C(CHI02,":",CMI02,":",CSI02)," hh:mm:ss")' else print,''

  print, n_skydip_c,   ' sky tip with short info ',format='(I5,A,$)'
  if n_skydip_c gt 0 then print, dt/n_skydip_c-0.5,format='(" i.e. every ",C(CHI02,":",CMI02,":",CSI02)," hh:mm:ss")' else print,''

  print, n_skydipX_c,   ' sky tip with long info ',format='(I5,A,$)'
  if n_skydipX_c gt 0 then print, dt/n_skydipX_c-0.5,format='(" i.e. every ",C(CHI02,":",CMI02,":",CSI02)," hh:mm:ss")' else print,''

  print, n_skytip_succes_c,   ' succesful sky tip calibrations ',format='(I5,A,$)'
  if n_skytip_succes_c gt 0 then print, dt/n_skytip_succes_c-0.5,format='(" i.e. every ",C(CHI02,":",CMI02,":",CSI02)," hh:mm:ss")' else print,''

ENDIF

; warnings if something went wrong ...
; with time
IF (t0 LT -1095) OR (t0 GT 1827933925) OR (tN LT -1095) OR (tN GT 1827933925) THEN BEGIN
  MESSAGE,'WARNING : INVALID JULLIAN DATE(S) !',/INFORMATIONAL
  print,'t[0]=',t0,' t[N-1]=',tN
  n_warnings += 1
ENDIF
IF (min(time) LT t0) OR (max(time) GT tN) THEN BEGIN
  MESSAGE, 'WARNING : inconsistency in time: min < t[0] or t[N-1] < max !',/INFORMATIONAL
  print,'min,max(t)=',min(time),',',max(time),'t[0]=',t0,' t[N-1]=',tN
  n_warnings += 1
ENDIF

IF (i_c GT n-1) AND (not EOF(unit)) then BEGIN 
  MESSAGE,'WARNING : more calibration records in file than indicated in the header !',/INFORMATIONAL
  n_warnings += 1
ENDIF

IF n_gain NE n_gain_c then BEGIN 
  MESSAGE,'WARNING : found'+string(n_gain_c)+' gain calibrations instead expected '+string(n_gain),/INFORMATIONAL
  n_warnings += 1
ENDIF

IF n_noise NE n_noise_c then BEGIN
  MESSAGE,'WARNING : found'+string(n_noise_c)+' noise calibrations instead expected '+string(n_noise),/INFORMATIONAL
  n_warnings += 1
ENDIF

IF n_skydip_c+n_skydipX_c NE n_skydip then BEGIN
  MESSAGE,'WARNING : found'+string(n_skydip_c+n_skydipX_c)+' sky tip calibrations instead expected '+string(n_skydip),/INFORMATIONAL
  n_warnings += 1
ENDIF



CLOSE,unit
FREE_LUN,unit

END
