from level1.rpg_bin import get_rpg_bin
from level1 import rpg_mwr
from level1.meta_nc import get_data_attributes
from level1.site_config import get_global_attributes
import numpy as np
from typing import Optional

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
        
    elif data_type == '1C01':
        rpg_bin = get_rpg_bin(path_to_files,'brt')
        rpg_bin.data['frequency'] = rpg_bin.header['_f']        
        _append_hkd(path_to_files,rpg_bin,data_type)
        
        rpg_irt = get_rpg_bin(path_to_files,'irt')
        rpg_bin.data['ir_wavelength'] = rpg_irt.header['_f']
        _add_interpol1(rpg_bin.data,rpg_irt.data['irt'],rpg_irt.data['time'],'irt')
        
        rpg_met = get_rpg_bin(path_to_files,'met')
        _add_interpol1(rpg_bin.data,rpg_met.data['air_temperature'],rpg_met.data['time'],'air_temperature')
        _add_interpol1(rpg_bin.data,rpg_met.data['relative_humidity'],rpg_met.data['time'],'relative_humidity')
        _add_interpol1(rpg_bin.data,rpg_met.data['air_pressure'],rpg_met.data['time'],'air_pressure')
        if (int(rpg_met.header['_n_sen'],2) & 1) != 0:
            _add_interpol1(rpg_bin.data,rpg_met.data['adds'][:,0],rpg_met.data['time'],'wind_speed')
            rpg_bin.data['wind_speed'] = rpg_bin.data['wind_speed'] / 3.6
        if (int(rpg_met.header['_n_sen'],2) & 2) != 0:
            _add_interpol1(rpg_bin.data,rpg_met.data['adds'][:,1],rpg_met.data['time'],'wind_direction')
        if (int(rpg_met.header['_n_sen'],2) & 4) != 0:
            _add_interpol1(rpg_bin.data,rpg_met.data['adds'][:,2],rpg_met.data['time'],'rain_rate')          
        
    else:
        raise RuntimeError(['Data type '+ data_type +' not supported for file writing.'])
        
    return rpg_bin
    
    
def _append_hkd(path_to_files: str, 
                rpg_bin: dict, 
                data_type: str) -> None:
    """Append hkd data on same time grid"""
    
    hkd = get_rpg_bin(path_to_files,'hkd')    
    _add_interpol1(rpg_bin.data,hkd.data['station_latitude'],hkd.data['time'],'station_latitude')
    _add_interpol1(rpg_bin.data,hkd.data['station_longitude'],hkd.data['time'],'station_longitude')   
    
    if data_type in ('1B01','1C01'):
        _add_interpol1(rpg_bin.data,hkd.data['temp'][:,0:2],hkd.data['time'],'t_amb')
        _add_interpol1(rpg_bin.data,hkd.data['temp'][:,2:4],hkd.data['time'],'t_rec')
        
        
def _add_interpol1(data0: dict, 
                   data1: np.ndarray, 
                   time1: np.ndarray,
                   output_name: str) -> None:
    
    if data1.ndim > 1:
        data0[output_name] = np.ones([len(data0['time']), data1.shape[1]], np.float32)*-999.
        for ndim in range(data1.shape[1]):
            data0[output_name][:,ndim] = np.interp(data0['time'],time1,data1[:,ndim])
    else:
        data0[output_name] = np.interp(data0['time'],time1,data1)