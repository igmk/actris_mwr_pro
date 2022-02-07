import numpy as np


params = {
    'station_altitude': 108., # m MSL
    'station_longitude': 6.407, # degrees east
    'station_latitude': 50.906, # degrees north

    # integration time BLB files: 
    # The integration time for a boundary layer scan (RPG-specific) must be specified here. 
    # This is only relevant if BLB files for boundary layer temperature profiles are created. 
    # You can find this number with the RPG software when you load your Measurement Definition File (MDF) and open "Definition of Measurement and Calibration Parameters" 
    # --> "Products and Integration" --> "Total Integration Time" --> "Brightness Temperatures (BL). 
    # Specify integration time in seconds.
    'scan_time' : 50.,

    # integration time of measurements in seconds
    'int_time': 1, 

    # bandwidth of the central frequency in GHz (center frequency of single of upper side-band)
    'bandwidth': np.array([230, 230, 230, 230, 230, 230, 230, 230, 230, 230, 230, 600, 1000, 2000]),

    # single, double, or double-double sideband (1, 2, 4 are possible values)
    'sideband_count': 1,

    # 56.xx +/- X +/- Y
    'sideband_IF_separation': np.array([[0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.], [0.]]),

    # offset correction for TBs, i.e. adjustement of observation to nominal frequency:
    # Note: RPG offers a frequency shift within the radiometer software. 
    # This shift will modify the TBs calculated with the RPG software. 
    # Original TBs CANNOT be reconstructed, unless you have recorded the voltages & calibration parameters and possess adequate software routines to perform a re-calibration.
    # Hence, the authors of actris_mwr_pro DO NOT recommend to apply the frequency shifts in general. 
    # If you have applied these frequency shifts, you may still give these to protocol so that they are contained in the resulting level1 netcdf files. 
    # The variable freq_shift specifies the frequency shifts applied [in MHz]. 
    # Set all to 0 if no frequency shifts were applied."""
    'freq_shift' : np.array([0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0.]),

    # [lower, upper] bound threshold of TB quality flags in K
    'TB_threshold': np.array([2.7, 330.]), 

    # solar angle flagging:
    # saf parameter lets the user determine within what angular region (in deg) around the sun the MWR TBs and products are flagged, 
    # e.g. saf = 5 means that all values with +-5 deg elevation and +-5 deg azimuth will be flagged. 
    # Note that this flag only works if radiometer position is correctly specified and the radiometer is aligned towards the North.
    'saf': 7.,

    # azimuth correction:
    # Azimuth angle is transformed to geographical coordinates (E=90 and W=270), currently only for RPG scanners.
    # Set az_cor to the angle that RPG software gives when instrument is pointing to the North.
    # If you do not want to transform the coordinates set azi_cor to -999.
    'azi_cor': 0.,

    # active quality flags for level 1 data; 1: flag active
    # Bit 1: Missing TB-value\n'
    # Bit 2: TB threshold (lower range)\n'
    # Bit 3: TB threshold (upper range)\n'
    # Bit 4: Spectral consistency threshold\n'
    # Bit 5: Receiver sanity\n'
    # Bit 6: Rain flag\n'
    # Bit 7: Solar flag\n'
    # Bit 8: TB offset threshold'
    'flags_active': [1, 1, 1, 1, 1, 1, 1, 0],

    # retrieval coefficient files
    # specify retrieval coefficient path and files to be used in retrieving Level 2 products here
    'path_coeff': '/home/hatpro/mwr_pro_jue/mwr_pro/retrievals/',
    'algo_lwp': ['deb_rt00_90'],
    'algo_iwv': ['deb_rt00_90'],
    'algo_tel': ['deb_rt00'],
    'algo_hze': ['deb_rt00_90', 'deb_rt00_80', 'deb_rt00_75', 'deb_rt00_71', 'deb_rt00_61', 'deb_rt00_60', 'deb_rt00_52', 'deb_rt00_45', 'deb_rt00_42', 'deb_rt00_32', 'deb_rt00_30', 'deb_rt00_23'],
    
    # spectral consistency coefficient files
    # specify spectral consistency coefficient path and files to be used in Level 1 quality control here
    'path_spec': '/home/hatpro/idl/tbx_ret/',
    'algo_spec': ['tbx_deb_rt00_90'],
    'threshold_spec': np.array([2., 1., 1., 1., 1., 1., 2., 3., 3.5, 4., 2.5, 1., 1.5, 1.]),
    'factor_spec': np.array([.5, 1., 1., 1., 1., 1., 1., 1.5, 1.5, 1., 1., 1., 1., 1.]),
    
}

global_specs = {
    # Name of the conventions followed by the dataset  
    'conventions' : 'CF-1.8',
    
    # A succinct description of what is in the dataset, composed of instrument type and site name
    'title' : 'Juelich RPG HATPRO G5',
    
    # Versioning of the datasets (containing date and software version)
    'history' : '',
    
    # Where the original data was produced
    'institution' : 'University of Cologne',
    
    # The method of production of the original data
    'source' : 'Ground Based Remote Sensing',
    
    # Miscellaneous Information about the dataset or methods used to produce it
    'comment' : '',
    
    # References that describe the data or methods used to produce it
    'references' : '',
    
    # Name of measurement station
    'site_location' : 'Juelich, Germany',
    
    # E-PROFILE instrument identifier
    'instrument_id' : '',
    
    # WIGOS Station identifier acording to WIGOS convention
    'wigos_station_id' : '0-276-0-10508',
    
    # Department responsible for the instrument
    'principal_investigator' : 'Institute for Geophysics and Meteorology (IGMK)',
    
    # Manufacturer of the instrument
    'instrument_manufacturer' : 'Radiometer Physics (RPG)',
    
    # Instrument model
    'instrument_model' : 'RPG-HATPRO',
    
    # Instrument generation
    'instrument_generation' : 'G5',
    
    # Specific to mainboard
    'instrument_hw_id' : '',
    
    # Name of network(s) that instrument may be part of
    'network_name' : 'ACTRIS',
    
    # Name of campaign instrument may collect data for
    'campaign_name' : '',
    
    # List of files the data set is depending on
    'dependencies' : '',
    
    # Data license
    'license' : '',
    
    # Instrument calibration status
    'instrument_calibration_status' : '',
    
    # Time of last (automatic or manual) absolute calibration; LN2 or sky tipping as YYYYMMDD
    'date_of_last_absolute_calibration' : '',
    
    # Time of last covariance update as YYYYMMDD
    'date_of_last_covariance_matrix' : '',
    
    # Description of the type of automatic calibrations performed including information at calibration interval and respective integration time
    'type_of_automatic_calibrations_performed' : '',
    
    # Logbook repair/replacement work performed
    'factory_history' : '',
    
    # Manufacturer of the infrared radiometer
    'ir_instrument_manufacturer' : '',
    
    # Infrared radiometer model
    'ir_instrument_model' : 'RPG-HATPRO',   
    
    # Fabrication year of the infrared radiometer'
    'ir_instrument_fabrication_year' : '',  
    
    # Manufacturer of the weather station
    'met_instrument_manufacturer' : '',
    
    # Weather station model
    'met_instrument_model' : '',    
    
    # Fabrication year of the weather station
    'met_instrument_fabrication_year' : '',     
    
    # Air temperature accuracy. Unit: K.
    'air_temperature_accuracy' : '',
    
    # Relative humidity accuracy. Unit: 1.
    'relative_humidity_accuracy' : '',
    
    # Air pressure accuracy. Unit: Pa.
    'air_pressure_accuracy' : '',
    
    # Rain rate accuracy. Unit: mm/h.
    'rain_rate_accuracy' : '',
    
    # Wind direction accuracy. Unit: degrees.
    'wind_direction_accuracy' : '',
    
    # Wind speed accuracy. Unit: m/s.
    'wind_speed_accuracy' : '',
}