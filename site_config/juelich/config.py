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
    'bandwidth': np.array([230., 230., 230., 230., 230., 230., 230., 
                           230., 230., 230., 230., 600., 1000., 2000.]),

    # single, double, or double-double sideband (1, 2, 4 are possible values)
    'sideband_count': 1,

    # 56.xx +/- X +/- Y
    'sideband_IF_separation': np.array([[0.], [0.], [0.], [0.], [0.], [0.], [0.], 
                                        [0.], [0.], [0.], [0.], [0.], [0.], [0.]]),
    
    # Beam width (3 dB) of the microwave radiometer
    'beam_width': -999.,

    # offset correction for TBs, i.e. adjustement of observation to nominal frequency:
    # Note: RPG offers a frequency shift within the radiometer software. 
    # This shift will modify the TBs calculated with the RPG software. 
    # Original TBs CANNOT be reconstructed, unless you have recorded the voltages & calibration parameters and possess adequate software routines to perform a re-calibration.
    # Hence, the authors of actris_mwr_pro DO NOT recommend to apply the frequency shifts in general. 
    # If you have applied these frequency shifts, you may still give these to protocol so that they are contained in the resulting level1 netcdf files. 
    # The variable freq_shift specifies the frequency shifts applied [in MHz]. 
    # Set all to 0 if no frequency shifts were applied.
    'freq_shift' : np.array([0., 0., 0., 0., 0., 0., 0., 
                             0., 0., 0., 0., 0., 0., 0.]),

    # [lower, upper] bound threshold of TB quality flags in K
    'TB_threshold': np.array([2.7, 330.]), 
    
    # IRT parameters
    'ir_bandwidth' : -999.,
    'ir_beamwidth': -999.,

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

    # quality flag status for level 1 data; 0: flag active
    # Bit 1: missing_tb
    # Bit 2: tb_below_threshold
    # Bit 3: tb_above_threshold
    # Bit 4: spectral_consistency_above_threshold
    # Bit 5: receiver_sanity_failed
    # Bit 6: rain_detected
    # Bit 7: sun_in_beam
    # Bit 8: tb_offset_above_threshold
    'flag_status': [0, 0, 0, 0, 0, 0, 0, 1],
    
    # factor for spectral consistency quality flag 
    # multiplied with channel TB retrieval error
    # compared to absolute difference of retrieved and observed brightness temperatures
    # (subtracting a 5min mean for noise reduction)
    'tbx_f': np.array([2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 
                        3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5]),     
    
    # thresholds for met quality flags
    'met_thresholds': np.array([[213.15, 333.15], # air_temperature [K]
                                [0., 1.],         # relative_humidity [1]
                                [900., 1100.],    # air_pressure [hPa]
                                [0., 300.],       # rain_rate [mm/h]
                                [0., 360.],       # wind_direction [°]
                                [0., 100.]]),     # wind_speed [m/s]
}


global_specs = {
    # Name of the conventions followed by the dataset  
    'conventions' : 'CF-1.8',
    
    # A succinct description of what is in the dataset, composed of instrument type and site name
    'title' : 'HATPRO G5 MWR at Juelich, Germany',
    
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
    
    # E-PROFILE instrument identifier. “A” if there is only one instrument on the station. Additional instruments are identified with the letters B, C, etc.
    'instrument_id' : '',
    
    # WIGOS Station identifier acording to WIGOS convention
    'wigos_station_id' : '0-276-0-10508',
    
    # Department responsible for the instrument
    'principal_investigator' : 'Institute for Geophysics and Meteorology (IGMK)',
    
    # Manufacturer of the instrument
    'instrument_manufacturer' : 'Radiometer Physics (RPG)',
    
    # Instrument model
    'instrument_model' : 'HATPRO',
    
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
    
    # Status of instrument absolute calibration
    'instrument_calibration_status' : '',
    
    # Time of last (automatic or manual) absolute calibration; LN2 or sky tipping as YYYYMMDD
    'date_of_last_absolute_calibration' : '',
    
    # Time of last covariance update as YYYYMMDD
    'date_of_last_covariance_matrix' : '',
    
    # Type of automatic calibrations including information on calibration interval and respective integration time
    'type_of_automatic_calibrations' : '',
    
    # Logbook repair/replacement work performed
    'instrument_history' : '',
    
    # Manufacturer of the infrared radiometer
    'ir_instrument_manufacturer' : '',
    
    # Infrared radiometer model
    'ir_instrument_model' : 'RPG-HATPRO',   
    
    # Fabrication year of the infrared radiometer
    'ir_instrument_fabrication_year' : '',  
    
    # Total absolute calibration uncertainty of infrared brightness temperature, one standard deviation. Unit: K
    'ir_accuracy' : '',    
    
    # Logbook repair/replacement work performed
    'ir_instrument_history' : '',    
    
    # Manufacturer of the weather station
    'met_instrument_manufacturer' : '',
    
    # Weather station model
    'met_instrument_model' : '',    
    
    # Fabrication year of the weather station
    'met_instrument_fabrication_year' : '',     
    
    # Logbook repair/replacement work performed
    'met_instrument_history' : '',      
    
    # Air temperature accuracy. Unit: K.
    'air_temperature_accuracy' : '',
    
    # Relative humidity accuracy. Unit: 1.
    'relative_humidity_accuracy' : '',
    
    # Air pressure accuracy. Unit: hPa.
    'air_pressure_accuracy' : '',
    
    # Rain rate accuracy. Unit: mm/h.
    'rain_rate_accuracy' : '',
    
    # Wind direction accuracy. Unit: degrees.
    'wind_direction_accuracy' : '',
    
    # Wind speed accuracy. Unit: m/s.
    'wind_speed_accuracy' : '',
}
