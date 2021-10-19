from level1.rpg_bin import get_rpg_bin
from level1 import rpg_mwr
from level1.meta_nc import get_data_attributes
from level1.site_config import get_global_attributes
import numpy as np

def lev1_to_nc(site: str,
               data_type: str,
               path_to_files: str, 
               output_file: str):
    """This function reads one day of RPG MWR binary files,
    adds attributes and writes it into netCDF file.
    
    Args:
        site: Name of site.
        data_type: Data type of the netCDF file.
        path_to_lwp_files: Folder containing one day of RPG MWR binary files.
        output_file: Output file name.
        
    Examples:
        >>> from level1.write_lev1_nc import lev1_to_nc
        >>> lev1_to_nc('site_name', '1B01', '/path/to/files/', 'rpg_mwr.nc')
    """

    rpg_bin = prepare_data(path_to_files,data_type)
    hatpro = rpg_mwr.Rpg(rpg_bin.data)
    get_data_attributes(hatpro.data,data_type)    
    global_attributes = get_global_attributes(site,data_type)
    rpg_mwr.save_rpg(hatpro, output_file, global_attributes, data_type)
    
    
def prepare_data(path_to_files: str, 
                  data_type: str) -> dict:    
    """Load and prepare data for netCDF writing"""
    
    if data_type == '1B01':
        rpg_bin = get_rpg_bin(path_to_files,'brt')
        rpg_bin.data['frequency'] = rpg_bin.header['_f']
        _append_hkd(path_to_files,rpg_bin,data_type)
        
    elif data_type == '1B11':
        rpg_bin = get_rpg_bin(path_to_files,'irt')
        rpg_bin.data['ir_wavelength'] = rpg_bin.header['_f']
        _append_hkd(path_to_files,rpg_bin,data_type)    

    elif data_type == '1B21':
        rpg_bin = get_rpg_bin(path_to_files,'met')
        if (int(rpg_bin.header['_n_sen'],2) & 1) != 0:
            rpg_bin.data['wind_speed'] = rpg_bin.data['adds'][:,0] / 3.6
        if (int(rpg_bin.header['_n_sen'],2) & 2) != 0:
            rpg_bin.data['wind_direction'] = rpg_bin.data['adds'][:,1]
        if (int(rpg_bin.header['_n_sen'],2) & 4) != 0:
            rpg_bin.data['rain_rate'] = rpg_bin.data['adds'][:,2]            
        _append_hkd(path_to_files,rpg_bin,data_type)
        
    else:
        raise RuntimeError(['Data type '+ data_type +' not supported for file writing.'])
        
    return rpg_bin
    
    
def _append_hkd(path_to_files: str, 
                rpg_bin: dict, 
                data_type: str) -> None:
    """Append hkd data on same time grid"""
    
    hkd = get_rpg_bin(path_to_files,'hkd')    
    rpg_bin.data['station_latitude'] = np.interp(rpg_bin.data['time'],hkd.data['time'],hkd.data['lat'])
    rpg_bin.data['station_longitude'] = np.interp(rpg_bin.data['time'],hkd.data['time'],hkd.data['lon'])    
    
    if data_type == '1B01':
        rpg_bin.data['t_amb'] = np.ones( [len(rpg_bin.data['time']),2], np.float32)*-999.
        rpg_bin.data['t_amb'][:,0] = np.interp(rpg_bin.data['time'],hkd.data['time'],hkd.data['temp'][:,0])
        rpg_bin.data['t_amb'][:,1] = np.interp(rpg_bin.data['time'],hkd.data['time'],hkd.data['temp'][:,1])
        rpg_bin.data['t_rec'] = np.ones( [len(rpg_bin.data['time']),2], np.float32)*-999.        
        rpg_bin.data['t_rec'][:,0] = np.interp(rpg_bin.data['time'],hkd.data['time'],hkd.data['temp'][:,2])
        rpg_bin.data['t_rec'][:,1] = np.interp(rpg_bin.data['time'],hkd.data['time'],hkd.data['temp'][:,3]) 