from global_nc import GLOBAL_ALL
import importlib


def get_site_specs(site: str, 
                   data_type: str) -> dict: 
    """ This function initializes site specific global attributes and parameters of RPG MWR Level 1 + 2 data for NetCDF file writing.
    Args:
        site: Name of site.
        data_type: Data type of the netCDF file.
        
    Returns:
        Dictionary
      
    Raises:
        RuntimeError: Specified site or data type is not supported.
    
    Example:
        from read_specs import get_site_specs
        global_attributes, params = get_site_specs('site', 'data_type')
       
    """

    site_dict = importlib.import_module(f"site_config." + site + ".config")
    global_specs = get_global_specs(site_dict.global_specs, data_type)
    # add_global_description(global_specs, GLOBAL_ALL)

    return global_specs, site_dict.params
            

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
            
            
def get_global_specs(global_specs: dict, 
                     data_type: str) -> dict:  
            
    vv = list(global_specs.values())[:]
    kk = list(global_specs.keys())[:]   
    
    # Level 2 specific attributes will be filled using retrieval file(s):
    global_l2 = {
    'retrieval_type' : '',
    'retrieval_elevation_angles' : '',
    'retrieval_frequencies' : '',
    'retrieval_auxiliary_input' : '',
    'retrieval_description': '',
    }      
        
    if data_type == '1B01':
        
        keys = kk[0:24]
        values = vv[0:24]

    elif data_type == '1B11':
        
        keys = kk[0:11] + kk[24:28] + kk[15:17] + kk[18:19] + kk[23:24]
        values = vv[0:11] + vv[24:28] + vv[15:17] + vv[18:19] + vv[23:24]

    elif data_type == '1B21':
        
        keys = kk[0:11] + kk[28:31] + kk[15:17] + kk[18:19] + kk[23:24] + kk[30:37]
        values = vv[0:11] + vv[28:31] + vv[15:17] + vv[18:19] + vv[23:24] + vv[30:37]

    elif data_type == '1C01':
        
        keys = kk[0:37]
        values = vv[0:37]
        
    elif data_type in ('2P01', '2P02', '2P03', '2P04', '2P07', '2P08', '2I01', '2I02', '2S02'):

        vv2 = list(global_l2.values())[:]
        kk2 = list(global_l2.keys())[:]          
        keys = kk[0:24] + kk2
        values = vv[0:24] + vv2
        
    else:
        raise RuntimeError(['Data type '+ data_type +' not supported for file writing.'])
        
    specs = {k: v for k, v in zip(keys, values)}
    
    return specs