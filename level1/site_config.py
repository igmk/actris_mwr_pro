from level1.global_nc import GLOBAL_ALL

def get_global_attributes(site: str, 
                          data_type: str) -> dict: 
    """ This function initializes site specific global attributes of RPG MWR Level 1 data for NetCDF file writing.
    Args:
        site: Name of site.
        data_type: Data type of the netCDF file.
        
    Returns:
        Dictionary
      
    Raises:
        RuntimeError: Specified site or data type is not supported.
    
    Example:
        from level1.site_config import get_global_attributes
        att = get_global_attributes('site_name','data_type')
       
    """

    specs = get_site_specs(site, data_type)
    add_global_description(specs, GLOBAL_ALL)
    
    return specs
            

def add_global_description(site_specs: dict, 
                           global_description: dict) -> None:
    """Adds global attribute description.
    Args:
        site_specs: Site specific global attributes.
        attributes: Global attribute description.
    """
    
    for key in site_specs:
        if key in global_description:
            site_specs[key] = site_specs[key]+' ['+''.join(global_description[key])+']'
            

def get_site_specs(site: str, 
                   data_type: str) -> dict:
    """Site specific global attributes"""

    if site == 'juelich':

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
        'wigos_station_id' : '',
        'principal_investigator' : 'IGMK',
        'instrument_manufacturer' : 'RPG',
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
        
        specs = global_specs

    else:
        raise RuntimeError(['Data type '+ data_type +' not supported for file writing.'])

    return specs