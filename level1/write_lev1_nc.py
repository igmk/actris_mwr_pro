from level1.site_config import get_site_att
from level1.rpg_bin import get_rpg_bin
from level1 import rpg_mwr
from level1.meta_nc import ATTRIBUTES_1B01

def lev1_to_nc(site: str,
               data_type: str,
               path_to_files: str, 
               output_file: str):
    """Converts RPG MWR data into Level 1b netCDF file.
    
    This function reads one day of RPG MWR binary files,
    concatenates the data and writes it into netCDF file.
    
    Args:
        site: Name of site
        path_to_lwp_files: Folder containing one day of RPG MWR binary files.
        output_file: Output file name.
        
    Examples:
        >>> from level1.write_lev1_nc import lev1_to_nc
        >>> lev1_to_nc('site_name', '/path/to/files/', 'rpg_mwr.nc')
    """

    """Load and prepare data for netCDF writing"""
    if data_type == '1B01':
        rpg_bin = get_rpg_bin(path_to_files,'brt')
        rpg_bin.data['frequency'] = rpg_bin.header['_f']
        
    else:
        raise RuntimeError(['Data type '+ data_type +' not supported for file writing.'])
    
    hatpro = rpg_mwr.Rpg(rpg_bin.data)
    rpg_mwr.update_attributes(hatpro.data,eval('ATTRIBUTES_'+ data_type))
    
    global_attributes = get_site_att(site,data_type)    
    rpg_mwr.save_rpg(hatpro, output_file, global_attributes)
