from level1.global_nc import GLOBAL

def get_site_att(site: str,
                 data_type: str) -> dict: 
    """ This function initializes site specific attributes of RPG MWR Level 1 data for NetCDF file writing.
    Args:
        site: Folder containing one day of a RPG MWR binary file type.
        
    Returns:
        Dictionary
      
    Raises:
        RuntimeError: Specified site is not supported.
    
    Example:
        from level1.site_config import get_site_att
        att = get_site_att('site_name')
       
    """

    spec = get_site_spec(site,data_type)
    add_global_description(spec,GLOBAL)
    
    return spec
            

def add_global_description(site_specs: dict, global_description: dict) -> None:
    """Adds global attribute description.
    Args:
        site_specs: Site specific global attributes.
        attributes: Global attribute description.
    """
    
    for key in site_specs:
        if key in global_description:
            site_specs[key] = site_specs[key]+' ['+''.join(global_description[key])+']'
            

def get_site_spec(site, data_type):
    """Specifiy site specific attributes"""

    if site == 'juelich':
        
        if data_type == '1B01':
        
            spec = {
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
            
        else:
            raise RuntimeError(['Data type '+ data_type +' not supported for file writing.'])

    else:
        raise RuntimeError(['Site '+ site +' not supported'])

    return spec