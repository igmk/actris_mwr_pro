;+
;**************
;PAR_MWR_PRO
;**************
; Abstract:
; * user defined parameters for MWR_PRO processing
; * PLEASE READ THIS FILE CAREFULLY
; Author:
; U. Loehnert
; Date:
; 2011-01-18
; Dependencies:
; -
; Changes:
; XXXX-XX-XX: ???
;-

;****output file naming - hdcp2 compliant; for details and options see hdcp2 observation product standard
;kind of measurement type
kkk = 'sups'
;supersite or owner institute
sss = 'joy'
;sss = 'rao'
;instrument or synergy product, retrieval algorithm
inst = 'mwr'
instBL = 'mwrBL'
;serial number of instrument or synergy product, retrieval algorithm
instnn = '00'
;level of data
lll1 = 'l1'
lll2 = 'l2'
;variable name
var_l1b = 'tb'
var_l1c = 'tb'
var_l2a_iwv = 'prw'
var_l2a_lwp = 'clwvi'
var_l2b_t = 'ta'
var_l2b_q = 'hua'
var_l2c_tel = 'ta'

;hdcp2 version of data level starting with v00
vnn_l1b = 'v01'
vnn_l1c = 'v01'

vnn_l2a_iwv = 'v01'
vnn_l2a_lwp = 'v01'
vnn_l2b_t = 'v01'
vnn_l2b_q = 'v01'
vnn_l2c_t = 'v01'

;OR internal (testing) version
;vnn_l2a_iwv = 'i00'
;vnn_l2a_lwp = 'i00'
;vnn_l2b_t = 'i00'
;vnn_l2b_q = 'i00'
;vnn_l2c_t = 'i00'

;*define file naming structures
file_naming_l1b = {kkk:kkk, sss:sss, inst:inst, instnn:instnn, lll1:lll1, var:var_l1b, vnn:vnn_l1b}
file_naming_l1c = {kkk:kkk, sss:sss, inst:instBL, instnn:instnn, lll1:lll1, var:var_l1c, vnn:vnn_l1c}
file_naming_2a_iwv = {kkk:kkk, sss:sss, inst:inst, instnn:instnn, lll2:lll2, var:var_l2a_iwv, vnn:vnn_l2a_iwv}
file_naming_2a_lwp = {kkk:kkk, sss:sss, inst:inst, instnn:instnn, lll2:lll2, var:var_l2a_lwp, vnn:vnn_l2a_lwp}
file_naming_2b_tze = {kkk:kkk, sss:sss, inst:inst, instnn:instnn, lll2:lll2, var:var_l2b_t, vnn:vnn_l2b_t}
file_naming_2b_hze = {kkk:kkk, sss:sss, inst:inst, instnn:instnn, lll2:lll2, var:var_l2b_q, vnn:vnn_l2b_q}
file_naming_2c_tel = {kkk:kkk, sss:sss, inst:instBL, instnn:instnn, lll2:lll2, var:var_l2c_tel, vnn:vnn_l2c_t}

;*dependencies of level2 data on version number of level1 data
dep_l2a_iwv = 'v01' ; l1b 
dep_l2a_lwp = 'v01' ; l1b
dep_l2b_t = 'v01' ; l1b
dep_l2b_q = 'v01' ; l1b
dep_l2c_t = 'v01' ; l1c

;****threshold of TB sanity flags in K
TB_min = 2.7
TB_max = 330.

;****thresholds for LWP & IWV sanity flags in kgm^-2
LWP_min = -0.2
LWP_max = 3.
IWV_min = 0.
IWV_max = 100.

;****thresholds for retrieved temperature sanity flags in K
T_min = 180.
T_max = 330.

;****RETRIEVAL coefficient file specification for multiple linear regression
;All retrieval files must be present in retrieval directory specified in setup.sh. 
;Within the retrieval files, all frequency and angle coefficients must be specified with monotonically
;increasing value!
;If you are using an RPG-formatted *.RET file set 
;ret_file = 'RET' 
;If you are using a netcdf retrieval file set 
ret_file = 'nc'
;NAMING convention:
;algo_ppp=xxx_yyyy_zzz
;ppp=parameter
;lwp: liquid water path (level2a product) --> uses TB and TB^2 as predictand
;iwv: integrated water vapor (level2a product) --> uses TB and TB^2 as predictand
;wdl: wet delay (level2a product) --> uses TB and TB^2 as predictand
;taa: total atmospheric attenuation (level2a product) --> uses TB and TB^2 as predictand
;tze: temperature profile using only one elevation angle (level2b product) --> uses TB and TB^2 as predictand
;hze: humidity profile using only one elevation angle (level2b product) --> uses TB and TB^2 as predictand
;tel: temperature profile from multiple elevation angles (boundary layer scanning mode --> level2c product) --> uses TB as predictand
;xxx=climatology specification
;yyyy=user defined description
;zzz=elevation angle information (i.e. 090: zenith observation, 005: close to horizon observation, BL1: spec. BL elevations)
;**NOTE for LEVEL2a products: zzz SHOULD NOT be specified here --> pl_mk will calculate the specified level1 parameter 
;using all retrieval files containing the string 'ppp_yyyy'. Each *.RET only contains retrieval file coefficients at one 
;specific elevation angle. Retrievals at arbitrary elevation angles are possible because regression coefficients are 
;linearly interpolated between elevation angles. From past experience, it is recommended that the elevation angles specified 
;in the retrieval files and contained in the measurements should not differ too much (e.g. more that 5 deg), 
;otherwise systematic errors due to interpolation may arise. IF only one *.RET file is found pl_mk will only calculate
;level1 parameters for measured elevation angles +.0.6 around the one contained in the *.RET file. 
;SPECIFICATION EXMAPLE for LEVEL1 product (all but TAA):
;algo_lwp = 'cab_mwr1'
;            ppp_yyyy       
;All retrieval files containing the string ppp_yyyy are used for calculating the level1 product. These files must ONLY
;differ in elevation angle - other parameters must stay the same.
;An exception are retrievals for total atmospheric attenuation (TAA). The are derived for one specified frequency and angle 
;The user must specify the array algo_taa as follows:
;SPECIFICATION EXMAPLE for LEVEL2 TAA product :
;algo_taa=['deb_0222_070', 'deb_0512_070']
;           ppp_yyyy_zzz
;yyyy is specified by the user - this case yyyy could be the frequency multiplied by 10 
;For each element of algo_taa specified, the user must also specify a frequency
;f_taa = [22.24, 51.28]
;**NOTE for level2b products: zzz SHOULD be specified here --> Only data containing exactly the specified elevation angle
;will be converted to level2b product. It is possible to define a string of retrieval parameters. Then level2b product will 
;again only be derived at these multiple defined angles. If measured elevation angle and retrieval elevation angle differ, 
;then the retrieval will be aborted.
;SPECIFICATION EXMAPLE for LEVEL2b product:
;algo_tze = 'jue_0112_090'
;algo_hze = ['jue_0112_090', 'jue_0112_045']
;**NOTE for level2c products: zzz SHOULD be specified here --> Only one retrieveal file can be specified and only data 
;containing exactly the specified elevation angle will be converted to a level2c product. 
;In the level2c (=tel) case zzz specifies a set of elevation angles and frequencies
;that are used to derive the temperature profile assuming horizontal homogeneous conditions. Only zenith profiles can be
;retrieved as a level2c product. If measured elevation angle and retrieval elevation angle differ, then the retrieval
;will be aborted. The 'tel' product consists of temperature, relative humidity, pot. temperature & equiv. potential 
;temperature profiles and is considered the most accurate profile product.
;SPECIFICATION EXMAPLE for LEVEL2c product:
;algo_tel = 'jue_0112_BL1'
;SPECIFY retrieval algorithm parameters staring here. Specify only if you wish the parameters to be retrieved within
;the MWR_PRO run
algo_iwv='deb_rt00'
algo_lwp='deb_rt00'
;algo_taa=['deb_r213_070', 'deb_r238_070', 'deb_r316_070'] 
;f_taa = [21.3, 23.8, 31.6]
algo_tze=['deb_rt00_90', 'deb_rt00_80', 'deb_rt00_75', 'deb_rt00_71', 'deb_rt00_61', 'deb_rt00_60', 'deb_rt00_52', 'deb_rt00_45', 'deb_rt00_42', 'deb_rt00_32', 'deb_rt00_30', 'deb_rt00_23']
algo_hze=['deb_rt00_90', 'deb_rt00_80', 'deb_rt00_75', 'deb_rt00_71', 'deb_rt00_61', 'deb_rt00_60', 'deb_rt00_52', 'deb_rt00_45', 'deb_rt00_42', 'deb_rt00_32', 'deb_rt00_30', 'deb_rt00_23']
algo_tel='deb_rt00'

;****raw data file type
;RPG format
rd_file='brt'
;RPG format
;rd_file='wvl'
;MWP format (not yet implemented)
;rd_file='xxx'
;RESCOM format
;rd_file = 'dec'

;****integration time BLB files
;The integration time for a boundary layer scan (RPG-specific) must be specified here. This is only relevant if BLB files
;for boundary layer temperature profiles are created. You can find this number with the RPG software when you load your
;Measurement Definition File (MDF) and open "Definition of Measurement and Calibration Parameters" --> "Products and Integrat
;--> "Total Integration Time" --> "Brightness Temperatures (BL)"
;Specify integration time in seconds
int_time_1c = 100.

;****azimuth time contours
; plot iwv azimuth-time contour when aztp_iwv=1
aztp_iwv = 1
; plot lwp azimuth-time contour when aztp_lwp=1
aztp_lwp = 1
; specify array of angles (up to two) in el_aztp for which azimuth-time contours shall be plotted
el_aztp = [30.]
; az_scan_thres is a criterion for finding azimuth scans in your data; if two sequential measurements differ by more than
; az_scan_thres, an scan is assumed (see aztp.pro for details)
az_scan_thres = 1.

;****files containing TB uncertainty specification
;Note: the following files MUST exist in the data directory specified in setup_new_meas.sh under data/uncertainty. It should be
;the responsibility of each operator to specifiy these uncertainties as accurately as possible.

;1.) The tb_bias_file must contain the variable tb_bias and contain an array corresponding to the number of frequency channels
;of the radiometer at stake. This variable is an estimate of the one-standard-deviation calibration error to be expected from an
;absolute calibration, i.e. the likely systematic error of brightness temperature. E.g. for the TOPHAT system, the numbers from
;Maschwitz et al. 2012, AMT (Tab. 5) are used. However, these numbers of course depend on each instrument and should be adapted accordingly.
;Optimally the file should contain the variable date_tb_bias='yyyymmdd', a variable 'instrument_tb_bias', and a variable 'comment_tb_bias'.
;The format of this file is IDL 'sav'.
tb_bias_file = 'tb_bias_jue_20130101.sav'
;tb_bias_file = 'tb_bias_xxx_20XXXXXX_mwp.sav'
;2.) The tb_cov_file must contain the variable tb_cov and contain a matrix corresponding to the number of frequency channels (n_freq x n_freq)
;of the radiometer at stake. This variable is calculated from brightness temperature observations of an internal black body whose physical
;temperature is known. The square root of the matrix diagonal gives the brightness temperature random error of each frequency channel.
;However, these numbers depend on each individual instrument and should be adapted accordingly.
;Optimally the file should contain the variable date_tb_cov='yyyymmdd', a variable 'instrument_tb_cov', and a variable 'comment_tb_cov'.
;The format of this file is IDL 'sav'.
tb_cov_file = 'tb_cov_jue_20130613.sav'
;tb_cov_file = 'tb_cov_xxx_20XXXXXX_mwp.sav'

;****offset correction for TBs, i.e. adjustement of observation to nominal frequency
;'**Note: RPG offers a frequency shift within the radiometer software. This shift will modify the TBs calculated with the RPG software. 
;         Original TBs CANNOT be reconstructed, unless you have recorded the voltages & calibration parameters and possess adequate
;         software routines to perform a re-calibration. Hence, the author of mwr_pro DOES NOT recommend to apply the frequency shifts
;         in general. If you have applied these frequency shifts, you may still give these to protocol so that they are contained in the
;         resulting level1 netcdf files. The variable freq_shift specifies the frequency shifts applied [in MHz]. Set all to 0 if no 
;         frequency shifts were applied.
freq_shift = [0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]      

;**Alteranatively, mwr_pro makes possible a TB offset correction in TB space. I.e., you may apply a TB offset correction to the original
;  TBs, whereby the offset is stored in the level1 netcdf files. In this way, the original TBs may be easily reconstructed. The TB offsets
;  are SUBTRACTED from the original level1 products. Level2 products will be calculated with the modified TBs if use_oc_l2 = 1.  
;  If certain measurements do not correspond to specified dates, frequencies and elevation angles for offset correction: 
;  offset correction is not applied for these measurements. Note, offset correction is linearly interpolated to elevation angles of measurements - however not for 1c data!

;**Note: tb_offset_file must be of IDL binary format (.sav) and contain the parameters defined in the following:
; date_offset=DBLARR(n_date_offset): specifies the times AFTER which a certain offset correction is valid (Julian day format!)
; freq_offset=DBLARR(n_freq_offset): frequency channels that offset correction is applied to (GHz)
; ang_offset=DBLARR(n_ang_offset): elevation angles that offset correction is applied to (elevation angle)
; tb_offset=DBLARR(n_date_offset, n_freq_offset, n_ang_offset): tb offset correction values (K) to be SUBTRACTED from 
; original values to ensure correct values

;  The variable offset_index specifies the channels to be corrected; 1: correct, 0: do not correct
;  Un-comment the following two lines in order to apply the offset correction  
tb_offset_file = 'tb_offset_jue.sav'
offset_index = [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1]
;offset_index = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
use_oc_l2 = 1; offset correction applied to calculate level2 data
;use_oc_l2 = 0; offset correction NOT applied to calculate level2 data (default)

;*This variable will be written to the netcdf output file
offset_tb_type = 'COSMO-DE analysis'
;offset_tb_type = '-'

;****lwp offset correction
; An automatic lwp offset correction can be set by setting set_lwp_off=1. In this case every
; processed day is analysed for clear sky situations. The criterion "clear sky" is set to yes/no 
; for every 20 minute interval based solely on the two minute zenith lwp-standard-deviation [in g/m2] 
; of the original lwp values. If, within the 20-minute-interval, each two-minute-interval shows a
; lwp-standard-deviation of lower than lwp_std_thres, then the 20-minute interval is considered as
; "clear-sky". Note, lwp_std_thres is instrument, retrieval and climate-zone dependent and thus must
; be chosen carefully! If set_lwp_off is set, a file lwp_offset_yyyy.sav is modified containing the
; dates and lwp offsets during clear sky intervals. Also, if set_lwp_off = 1 the the offset corrected lwp
; is automatically shown in the level2a quicklooks and an additional variable (lwp_cor) in written to the 
; netcdf files. Further, it is important to note, that the lwp offset correction is only applied to cases when
; flag_2a ist 0, meaning the measurement and LWP retrieval is flagged as "OK". If, i.e. the LWP offset during 
; clear sky is smaller than lwp_min (see above), then these data are NOT used for the lwp offset correction. 
; In this case the user should check the calibration or update the retrieval file.
set_lwp_off = 1
lwp_std_thres = 2.

;****plotting ranges
;all: all derived retrieval values are plotted
;flag: only those derived retrieval values are plotted which are flagged as "OK"
plot_range = 'all'
;plot_range = 'flag'
;user_range_min and user_range_max are used for y-axis ranges if set to values NE -999.
;otherwise minimum and maximum are used
;index   :  0          1          2      3    4         5
;variable: iwv[kgm^-2] lwp[gm^-2] wdl[m] T[K] q[gm^-3]  RH[%]
user_range_min = [-999., -100., -999., -999., -999., 0.]
user_range_max = [-999., 999., -999., -999., -999., 100.] 
;****contour line distances
;index:       0     1   2 
;variable:    T[K]  q[gm^-3]    RH[%]
cont_dist =   [2.,  1., 10.]

;*****IRT plotting ranges
irt_range_max = 10.
irt_range_min = -60.

;****elevation angle plotting range for quicklooks of data time series
;ang_low = 69.
;ang_high = 71.
ang_low = 89.
ang_high = 91.

;****height range in m for plotting T & q profile time series
h_max = 5000.
h_min = 0.

;****temporal resolution (in seconds) for plotting 2b products (nearest neighbor)
deltat2b = 600.

;****temporal resolution (in seconds) for plotting 2c products (nearest neighbor)
deltat2c = 1200.

;****extended screen output
verbose = 1   ; yes
;verbose = 0  ; no

;****general information about the measurements that will be archived as meta data in HD(CP)2-compliant netcdf format
;!Note: if no information is given, the variables must contain a character with at least one blank.
institution = 'Research Center Juelich Institute for Energy and Climate Research IEK8'
contact_person = 'Ulrich Loehnert (loehnert@meteo.uni-koeln.de)'
author = 'Ulrich Loehnert (loehnert@meteo.uni-koeln.de)'
source = 'RPG HATPRO G2 microwave radiometer TOPHAT'
history = 'Data processed with mwrpro-v04 by University of Cologne; '
conventions = 'CF-1.6'
license = 'For non-commercial use only. This data is subject to the HD(CP)2 data policy to be found at hdcp2.zmaw.de and in the HD(CP)2 Observation Data Product Standard.
measurement_site = 'JOYCE Juelich Observatory for Cloud Evolution'
comments = ' '
global_atts = {i:institution, cp:contact_person, s:source, h:history, cv:conventions, l:license, au:author, co:comments, ms:measurement_site}

latitude = 50.909 ; degrees north
longitude = 6.414 ; degrees east
altitude = 111 ; m MSL

;meas_name = 'MWP@RAO'
;institution = 'Richard Assmann Observatory Lindenberg (DWD)'
;contact_person = 'Dr. Juergen Guelnder (juergen.guelndner@dwd.de'
;source = 'Radiometrics MWP'
;history = 'Data processed with mwr_pro_v04 by University of Cologne.'
;conventions = 'based on CF-1.6'
;license = 'For non-commercial use only - this data is subject to the HD(CP)2 data policy to be found at hdcp2.zmaw.de and in the HD(CP)2 Observation Data Product Standard.
;measurement_site = 'Richard Assmann Observatory Lindenberg (DWD)'
;comment = 'Output data complies with HD(CP)2 Observation Data Product Standard.'

;global_atts = {i:institution, cp:contact_person, s:source, h:history, cv:conventions, l:license, ms: measurement_site, co:comment}

;latitude = 52.209 ; degrees north
;longitude = 14.122 ; degrees east
;altitude = 73 ; m MSL

;!!!!!the general information for the following, once already processed MWR data, need to be made HD(CP)2 compliant...
;meas_name = 'RESCOM MWR deployed at Cabauw by A. Martelucci'
;instrument = 'RESCOM7'
;institution = 'ESA/ESTEC'
;www = ' '
;latitude = '51.971' 
;longitude = '4.927' 
;altitude = '1m MSL' 

;meas_name = 'ATPROP MWR deployed at Cabauw by A. Martelucci'
;instrument = 'ATPROP'
;institution = 'ESA/ESTEC'
;www = ' '
;latitude = '51.971' 
;longitude = '4.927' 
;altitude = '1m MSL' 

;meas_name = 'SUNHAT Barbados deployment starting Dec. 2010'
;instrument = 'HATPRO-SUNHAT'
;institution = 'IGMK'
;www = 'http://www.mpimet.mpg.de/en/wissenschaft/atmosphaere-im-erdsystem/initiativen/barbadosstation.html'
;latitude = '13.165' 
;longitude = '59.432'
;altitude = '25m MSL'
;!!!!!

;****fixed azimuth angle
;If the azimuth angle of observation is not specified in the raw data files (i.e. in the RESCOM data), 
;az_fix is used to characterize all measurements.
az_fix = 0.

;****azimuth correction
;Azimuth angle is transformed to geographical coordinates (E=90 and W=270), currently only aplicable for RPG scanners.
;Set az_cor to the angle that RPG software gives when instrument is pointing to the North.
;If you do not want to transform the coordinates set az_cor to -999.
az_cor = 0.

;****plot radar data
;***yes (radar=1), no (radar=0)
radar = 0
;***path mmclx data
data_path_radar = '/data/joyce/data/mmclx'
;***file contains times with radar problems 
filter_file_radar = 'radar_filter.dat'

;****plot ceilo data
;Ceilometer data are plotted in a multi-plot together with IRT TBs and IWV/LWP (each one only if available).
;Raw ceilometer data (corresponding to the manufactorer type) should be stored under the ceilo_path/YYYYMMDD directory.
ceilo_path = '/net/aure/ceilo/data/'
;Choose ceilo_type between 'ct25' (Vaisala CT25K) and 'jen' (Jenoptik) 
ceilo_type = 'ct25'

;****overide filter file condition
;If offc=1 then ncdf-files are written no matter what entries are contained in filer_jue.dat
;If not, netcdf files are only created when quicklook-data has been checked by eye
offc=1

;****solar angle flagging
;PL_MK calculate solar zenith and azimuth angles as a function of latitude, longitude and date/time.
;saf parameter lets the user determine within what angular region (in deg) around the sun the 
;MWR TBs and products are flagged, e.g. saf = 5 means that all values with +-5 deg elevation and +-5 deg
;azimuth will be flagged; note that this flag only works if radiometer position is correctly specified and the 
;radiometer is aligned towards the North 
saf = 7.

;****end of parameter file 
