from level1.global_nc import GLOBAL_1B01, GLOBAL_1B11, GLOBAL_1B21

def get_global_attributes(site: str, 
                          data_type: str) -> dict: 
    """ This function initializes site specific global attributes of RPG MWR Level 1 data for NetCDF file writing.
    Args:
        site: Folder containing one day of a RPG MWR binary file type.
        data_type: Data type of the netCDF file.
        
    Returns:
        Dictionary
      
    Raises:
        RuntimeError: Specified site or data type is not supported.
    
    Example:
        from level1.site_config import get_site_att
        att = get_site_att('site_name','data_type')
       
    """

    specs = get_site_specs(site, data_type)
    add_global_description(specs, eval('GLOBAL_' + data_type))
    
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
        
        if data_type == '1B01':
        
            specs = {
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
            }
            
        elif data_type == '1B11':
            
            specs = {
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
            'ir_instrument_manufacturer' : '',
            'ir_instrument_model' : 'RPG-HATPRO',
            'instrument_generation' : 'G5',
            'instrument_hw_id' : '',
            'network_name' : 'ACTRIS',
            'campaign_name' : '',
            'license' : '',
            'factory_history' : '',
            'ir_instrument_fabrication_year' : '',
            }          
            
        elif data_type == '1B21':
            
            specs = {
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
            'met_instrument_manufacturer' : '',
            'met_instrument_model' : '',
            'network_name' : 'ACTRIS',
            'campaign_name' : '',
            'license' : '',
            'factory_history' : '',
            'air_temperature_accuracy' : '',
            'relative_humidity_accuracy' : '',
            'air_pressure_accuracy' : '',
            'rain_rate_accuracy' : '',
            'wind_direction_accuracy' : '',
            'wind_speed_accuracy' : '',            
            'met_instrument_fabrication_year' : '',
            }             
        
            
        else:
            raise RuntimeError(['Data type '+ data_type +' not supported for file writing.'])

    else:
        raise RuntimeError(['Site '+ site +' not supported'])

    return specs