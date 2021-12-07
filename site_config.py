from global_nc import GLOBAL_ALL
import numpy as np

def get_site_specs(site: str, 
                   data_type: str) -> dict: 
    """ This function initializes site specific global attributes and parameters of RPG MWR Level 1 data for NetCDF file writing.
    Args:
        site: Name of site.
        data_type: Data type of the netCDF file.
        
    Returns:
        Dictionary
      
    Raises:
        RuntimeError: Specified site or data type is not supported.
    
    Example:
        from level1.site_config import get_global_attributes
        att, param = get_site_specs('site_name','data_type')
       
    """

    specs, params = site_specs(site, data_type)
    add_global_description(specs, GLOBAL_ALL)
    
    return specs, params
            

def add_global_description(site_specs: dict, 
                           global_description: dict) -> None:
    """Adds global attribute description.
    Args:
        site_specs: Site specific global attributes.
        attributes: Global attribute description.
    """
    
    for key in site_specs:
        if key in global_description:
            site_specs[key] = site_specs[key]+'    ['+''.join(global_description[key])+']'
            

def site_specs(site: str, 
               data_type: str) -> dict:
    """Site specific global attributes and parameters"""

    if site == 'juelich':
        
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
            
            # Bandwidth of the central frequency in GHz (center frequency of single of upper side-band)
            'bandwidth': np.array([230, 230, 230, 230, 230, 230, 230, 230, 230, 230, 230, 600, 1000, 2000]),
            
            # Single, double, or double-double sideband (1, 2, 4 are possible values)
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
            'flags_active': [1,1,1,0,1,1,1,0],
            
            # Retrieval coefficient files
            # Specify retrieval algorithm parameter path and files to be used in retrieving Level 2 products here.
            'path_coeff': '/home/hatpro/mwr_pro_jue/mwr_pro/retrievals/',
            'algo_lwp': ['deb_rt00_90'],
            'algo_iwv': ['deb_rt00_90'],
            'algo_tze': ['deb_rt00_90', 'deb_rt00_80', 'deb_rt00_75', 'deb_rt00_71', 'deb_rt00_61', 'deb_rt00_60', 'deb_rt00_52', 'deb_rt00_45', 'deb_rt00_42', 'deb_rt00_32', 'deb_rt00_30', 'deb_rt00_23'],
            'algo_hze': ['deb_rt00_90', 'deb_rt00_80', 'deb_rt00_75', 'deb_rt00_71', 'deb_rt00_61', 'deb_rt00_60', 'deb_rt00_52', 'deb_rt00_45', 'deb_rt00_42', 'deb_rt00_32', 'deb_rt00_30', 'deb_rt00_23'],
        }

        global_specs = {
        'conventions' : 'CF-1.8',
        'title' : 'Juelich RPG HATPRO G5',
        'history' : '',
        'institution' : 'University of Cologne',
        'source' : 'Ground Based Remote Sensing',
        'comment' : '',
        'references' : '',
        'site_location' : 'Juelich, Germany',
        'instrument_id' : '',
        'wigos_station_id' : '0-276-0-10508',
        'principal_investigator' : 'Institute for Geophysics and Meteorology (IGMK)',
        'instrument_manufacturer' : 'Radiometer Physics (RPG)',
        'instrument_model' : 'RPG-HATPRO',
        'instrument_generation' : 'G5',
        'instrument_hw_id' : '',
        'network_name' : 'ACTRIS',
        'campaign_name' : '',
        'dependencies' : '',
        'license' : '',
        'instrument_calibration_status' : '',
        'date_of_last_absolute_calibration' : '',
        'date_of_last_covariance_matrix' : '',
        'type_of_automatic_calibrations_performed' : '',
        'factory_history' : '',
        'ir_instrument_manufacturer' : '',
        'ir_instrument_model' : 'RPG-HATPRO',     
        'ir_instrument_fabrication_year' : '',   
        'met_instrument_manufacturer' : '',
        'met_instrument_model' : '',    
        'met_instrument_fabrication_year' : '',        
        'air_temperature_accuracy' : '',
        'relative_humidity_accuracy' : '',
        'air_pressure_accuracy' : '',
        'rain_rate_accuracy' : '',
        'wind_direction_accuracy' : '',
        'wind_speed_accuracy' : '',                  
        # Level 2 specific:
        'retrieval_type' : '',
        'retrieval_elevation_angles' : '',
        'retrieval_frequencies' : '',
        'retrieval_auxiliary_input' : '',
        'retrieval_description': '',
        } 
        

    else:
        raise RuntimeError(['Site '+ site +' not supported'])        
        
    if data_type == '1B01':
        
        values = list(global_specs.values())[0:24]
        keys = list(global_specs.keys())[0:24]
        specs = {k: v for k, v in zip(keys, values)}

    elif data_type == '1B11':
        
        vv = list(global_specs.values())[:]
        kk = list(global_specs.keys())[:]
        keys = kk[0:11] + kk[24:26] + kk[15:17] + kk[18:19] + kk[23:24] + kk[26:27]
        values = vv[0:11] + vv[24:26] + vv[15:17] + vv[18:19] + vv[23:24] + vv[26:27]
        specs = {k: v for k, v in zip(keys, values)}

    elif data_type == '1B21':
        
        vv = list(global_specs.values())[:]
        kk = list(global_specs.keys())[:]
        keys = kk[0:11] + kk[27:29] + kk[15:17] + kk[18:19] + kk[23:24] + kk[30:36] + kk[29:30]
        values = vv[0:11] + vv[27:29] + vv[15:17] + vv[18:19] + vv[23:24] + vv[30:36] + vv[29:30]
        specs = {k: v for k, v in zip(keys, values)}

    elif data_type == '1C01':
        
        vv = list(global_specs.values())[:]
        kk = list(global_specs.keys())[:]
        keys = kk[0:36]
        values = vv[0:36]
        specs = {k: v for k, v in zip(keys, values)}
        
    elif data_type in ('2P00','2I00','2S00'):
        
        vv = list(global_specs.values())[:]
        kk = list(global_specs.keys())[:]
        keys = kk[0:11] + kk[37:41]
        values = vv[0:11] + vv[37:41]
        specs = {k: v for k, v in zip(keys, values)}        

    else:
        raise RuntimeError(['Data type '+ data_type +' not supported for file writing.'])

    return specs, params