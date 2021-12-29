;+
;*************************
PRO GET_COEFF_LEVEL2C_NC,$
;*************************
;INPUT
filename,$             ;retrieval file
;OUTPUT
ret_par,$              ;retrieval parameters of MLR 
angles,$               ;elevation angles used for BL scanning
freq,$                 ;frequency array in GHz
height_grid,$          ;height grid in m
offset,$               ;MLR offset (n_z)
c_l,$                  ;linear retrieval coefficients (n_z x (n_f + n_bl*n_ang-1))
freq_bl,$              ;freqeuncies used for BL scanning
predictand_err,$       ;standard error of predictand
;KEYWORDS
verbose=verbose
; Abstract:
; * read multi-linear retrieval (MLR) coefficients for T, q line-of-sight retrieval
;   from netcdf format
; Author:
; U. Loehnert
; Date:
; 2014-10-30
; Dependencies:
; -
; Changes:
; XXXX-XX-XX (UL): nothing yet
;-

;****open file
infid = NCDF_OPEN(filename, /NOWRITE)

;****identifiy variables
freqid = NCDF_VARID(infid, 'freq')
freqblid = NCDF_VARID(infid, 'freq_bl')
latid = NCDF_VARID(infid, 'lat')
lonid = NCDF_VARID(infid, 'lon') 
hgid = NCDF_VARID(infid, 'height_grid') 
prdmxid = NCDF_VARID(infid, 'prdmx')
prdmnid = NCDF_VARID(infid, 'prdmn')
prrmxid = NCDF_VARID(infid, 'prrmx')
prrmnid = NCDF_VARID(infid, 'prrmn')
aslid = NCDF_VARID(infid, 'asl')
elpredictandid = NCDF_VARID(infid, 'elevation_predictand')
elpredictorid = NCDF_VARID(infid, 'elevation_predictor')
predictorerrid = NCDF_VARID(infid, 'predictor_err')
surferrid = NCDF_VARID(infid, 'surface_err')
predictanderrid = NCDF_VARID(infid, 'predictand_err')
predictanderrsysid = NCDF_VARID(infid, 'predictand_err_sys')
coeffid = NCDF_VARID(infid, 'coefficient_mvr')
offsetid = NCDF_VARID(infid, 'offset_mvr')

;****get variables
NCDF_VARGET, infid, freqid, freq
NCDF_VARGET, infid, freqblid, freq_bl
NCDF_VARGET, infid, latid, lat
NCDF_VARGET, infid, lonid, lon
NCDF_VARGET, infid, hgid, height_grid
NCDF_VARGET, infid, prdmxid, prdmx
NCDF_VARGET, infid, prdmnid, prdmn
NCDF_VARGET, infid, prrmxid, prrmx
NCDF_VARGET, infid, prrmnid, prrmn
NCDF_VARGET, infid, aslid, asl
NCDF_VARGET, infid, elpredictandid, elevation_predictand
NCDF_VARGET, infid, elpredictorid, elevation_predictor
NCDF_VARGET, infid, predictorerrid, predictor_err
NCDF_VARGET, infid, surferrid, surface_err
NCDF_VARGET, infid, predictanderrid, predictand_err
NCDF_VARGET, infid, predictanderrsysid, predictand_err_sys
NCDF_VARGET, infid, coeffid, coefficient_mvr
NCDF_VARGET, infid, offsetid, offset_mvr

;****get global attributes
NCDF_ATTGET, infid, /GLOBAL, 'processing_date', com
processing_date = STRING(BYTE(com))
NCDF_ATTGET, infid, /GLOBAL, 'rt_data_path', com
rt_data_path = STRING(BYTE(com))
NCDF_ATTGET, infid, /GLOBAL, 'site', com
site = STRING(BYTE(com))
NCDF_ATTGET, infid, /GLOBAL, 'number_of_profiles_used', com
number_of_profiles_used = com
NCDF_ATTGET, infid, /GLOBAL, 'date_start', com
date_start = com
NCDF_ATTGET, infid, /GLOBAL, 'date_end', com
date_end = com
NCDF_ATTGET, infid, /GLOBAL, 'rt_calc_cut_off_height_in_m', com
rt_calc_cut_off_height_in_m = com
NCDF_ATTGET, infid, /GLOBAL, 'gas_absorption_model', com
gas_absorption_model = STRING(BYTE(com))
NCDF_ATTGET, infid, /GLOBAL, 'cloud_absorption_model', com
cloud_absorption_model = STRING(BYTE(com))
NCDF_ATTGET, infid, /GLOBAL, 'wv_linewidth', com
wv_linewidth = STRING(BYTE(com))
NCDF_ATTGET, infid, /GLOBAL, 'wv_continuum_correction', com
wv_continuum_correction = STRING(BYTE(com))
NCDF_ATTGET, infid, /GLOBAL, 'air_mass_correction', com
air_mass_correction = STRING(BYTE(com))
NCDF_ATTGET, infid, /GLOBAL, 'predictand', com
predictand = STRING(BYTE(com))
NCDF_ATTGET, infid, /GLOBAL, 'predictand_unit', com
predictand_unit = STRING(BYTE(com))
NCDF_ATTGET, infid, /GLOBAL, 'predictor', com
predictor = STRING(BYTE(com))
NCDF_ATTGET, infid, /GLOBAL, 'predictor_unit', com
predictor_unit = STRING(BYTE(com))
NCDF_ATTGET, infid, /GLOBAL, 'retrieval_version', com
retrieval_version = STRING(BYTE(com))
NCDF_ATTGET, infid, /GLOBAL, 'surface_mode', com
surface_mode = STRING(BYTE(com))
NCDF_ATTGET, infid, /GLOBAL, 'regression_type', com
regression_type = STRING(BYTE(com))
NCDF_ATTGET, infid, /GLOBAL, 'cloudy_clear', com
cloudy_clear = STRING(BYTE(com))
NCDF_ATTGET, infid, /GLOBAL, 'cloud_diagnosis', com
cloud_diagnosis = STRING(BYTE(com))
NCDF_ATTGET, infid, /GLOBAL, 'cloud_diagnosis_rh_threshold', com
cloud_diagnosis_rh_threshold = com
NCDF_ATTGET, infid, /GLOBAL, 'bandwidth_correction', com
bandwidth_correction = STRING(BYTE(com))
NCDF_ATTGET, infid, /GLOBAL, 'beamwidth_correction', com
beamwidth_correction = STRING(BYTE(com))

NCDF_CLOSE, infid

;**for level2c products: elevation_predictand 90deg zenith
;                        elevation_predictor variable
angles = elevation_predictor
n_ang = N_ELEMENTS(elevation_predictor)
n_freq = N_ELEMENTS(freq)
n_freq_bl = N_ELEMENTS(freq_bl)
offset = offset_mvr

IF N_ELEMENTS(regression_type) EQ 0 THEN BEGIN
 print, 'regression_type not specified'
 print, 'ABORT in GET_COEFF_LEVEL_2C_NC'
 stop
ENDIF

n_ret = N_ELEMENTS(height_grid)
c_l = coefficient_mvr(*, 0:(n_freq_bl*(n_ang-1) + n_freq)-1)

IF surface_mode EQ 'in_surface' THEN BEGIN
 print, 'GET_COEFF_LEVEL2C_NC: Warning - surface_mode in_surface not included in tel processing!' 
ENDIF

;**define retrieval parameter structure (ret_par)

ret_par = {process:processing_date, path:rt_data_path, s:site, np:number_of_profiles_used, sd:date_start, ed:date_end,$
           co:rt_calc_cut_off_height_in_m, ga:gas_absorption_model, ca:cloud_absorption_model, wvl:wv_linewidth,$
           wvc:wv_continuum_correction, amc:air_mass_correction, pd:predictand, pdu:predictand_unit, pr:predictor,$
           pru:predictor_unit, vers:retrieval_version, surf:surface_mode, reg:regression_type, clcl:cloudy_clear,$
           cd:cloud_diagnosis, cdthres:cloud_diagnosis_rh_threshold, bdw:bandwidth_correction, bmw:beamwidth_correction}

END

