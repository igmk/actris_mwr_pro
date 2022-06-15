from read_specs import get_site_specs
from level1.rpg_bin import get_rpg_bin
from level1.lev1_meta_nc import get_data_attributes
from level1.quality_control import apply_qc
from level1.met_quality_control import apply_met_qc
from utils import isbit
import rpg_mwr
import numpy as np
from typing import Optional
import glob

Fill_Value_Float = -999.
Fill_Value_Int = -99  


def lev1_to_nc(site: str,
               data_type: str,
               path_to_files: str, 
               path_to_prev: str,
               path_to_next: str,
               output_file: str):
    """This function reads one day of RPG MWR binary files,
    adds attributes and writes it into netCDF file.
    
    Args:
        site: Name of site.
        data_type: Data type of the netCDF file.
        path_to_files: Folder containing one day of RPG MWR binary files.
        path_to_prev: Folder containing previous day RPG MWR binary files.
        path_to_next: Folder containing next day RPG MWR binary files.
        output_file: Output file name.
        
    Examples:
        >>> from level1.write_lev1_nc import lev1_to_nc
        >>> lev1_to_nc('site_name', '1B01', '/path/to/files/', '/path/to/previous/day/', '/path/to/next/day/', 'rpg_mwr.nc')
    """


    file_list_hkd = get_file_list(path_to_files, path_to_prev, path_to_next, 'hkd')
    global_attributes, params = get_site_specs(site, data_type)
    rpg_bin = prepare_data(path_to_files, path_to_prev, path_to_next, data_type, params)
    if data_type in ('1B01', '1C01'):
        apply_qc(site, rpg_bin.data, params)    
    if data_type in ('1B21', '1C01'):
        apply_met_qc(rpg_bin.data, params)
    hatpro = rpg_mwr.Rpg(rpg_bin.data)
    hatpro.find_valid_times()
    hatpro.data = get_data_attributes(hatpro.data, data_type)
    rpg_mwr.save_rpg(hatpro, output_file, global_attributes, data_type, params, site)
    
    
def get_file_list(path_to_files: str, 
                  path_to_prev: str,
                  path_to_next: str,
                  extension: str):
    
    f_list = sorted(glob.glob(path_to_files + '*.' + extension))
    f_list_p = sorted(glob.glob(path_to_prev + '*.' + extension))
    f_list_n = sorted(glob.glob(path_to_next + '*.' + extension))
    if len(f_list) == 0:
        raise RuntimeError(['Error: no binary files with extension ' + extension + ' found in directory ' + path_to_files])
    if any(f_list_p) & any(f_list_n):
        f_list = [f_list_p[-1], *f_list, f_list_n[0]]        
    elif any(f_list_p) & len(f_list_n) == 0:    
        f_list = [f_list_p[-1], *f_list]
    elif len(f_list_p) == 0 & any(f_list_n):    
        f_list = [*f_list, f_list_n[0]]  
    return f_list
      
    
def prepare_data(path_to_files: str, 
                 path_to_prev: str,
                 path_to_next: str,
                 data_type: str,
                 params: dict) -> dict:    
    """Load and prepare data for netCDF writing"""
    
    if data_type in ('1B01','1C01'):
        file_list_brt = get_file_list(path_to_files, path_to_prev, path_to_next, 'brt')
        rpg_bin = get_rpg_bin(file_list_brt)
        rpg_bin.data['frequency'] = rpg_bin.header['_f']
        fields = ['bandwidth', 'n_sidebands', 'sideband_IF_separation', 'freq_shift', 'receiver_nb', 'receiver']
        for name in fields:
            rpg_bin.data[name] = params[name]
        rpg_bin.data['time_bnds'] = add_time_bounds(rpg_bin.data['time'], params['int_time'])
        rpg_bin.data['pointing_flag'] = np.zeros(len(rpg_bin.data['time']), np.int32)

        file_list_blb = get_file_list(path_to_files, path_to_prev, path_to_next, 'blb')
        file_list_hkd = get_file_list(path_to_files, path_to_prev, path_to_next, 'hkd')
        rpg_hkd = get_rpg_bin(file_list_hkd)
        rpg_blb = get_rpg_bin(file_list_blb)
        _add_blb(rpg_bin, rpg_blb, rpg_hkd, params)

        if params['azi_cor'] != -999.:
            _azi_correction(rpg_bin.data, params)

        if data_type == '1C01':
            try:
                file_list_irt = get_file_list(path_to_files, path_to_prev, path_to_next, 'irt')               
                rpg_irt = get_rpg_bin(file_list_irt)
                rpg_irt.data['irt'][rpg_irt.data['irt'] < 150.] = Fill_Value_Float
                rpg_bin.data['ir_wavelength'] = rpg_irt.header['_f']
                rpg_bin.data['ir_bandwidth'] = params['ir_bandwidth']
                rpg_bin.data['ir_beamwidth'] = params['ir_beamwidth']                    
                _add_interpol1(rpg_bin.data, rpg_irt.data['irt'], rpg_irt.data['time'], 'irt')
                _add_interpol1(rpg_bin.data, rpg_irt.data['ir_ele'], rpg_irt.data['time'], 'ir_ele')
                _add_interpol1(rpg_bin.data, rpg_irt.data['ir_azi'], rpg_irt.data['time'], 'ir_azi')
            except:
                print(['No binary files with extension irt found in directory ' + path_to_files])

            file_list_met = get_file_list(path_to_files, path_to_prev, path_to_next, 'met')                  
            rpg_met = get_rpg_bin(file_list_met)
            _add_interpol1(rpg_bin.data, rpg_met.data['air_temperature'], rpg_met.data['time'], 'air_temperature')
            _add_interpol1(rpg_bin.data, rpg_met.data['relative_humidity'], rpg_met.data['time'], 'relative_humidity')
            _add_interpol1(rpg_bin.data, rpg_met.data['air_pressure'], rpg_met.data['time'], 'air_pressure')
            if (int(rpg_met.header['_n_sen'],2) & 1) != 0:
                _add_interpol1(rpg_bin.data, rpg_met.data['adds'][:,0], rpg_met.data['time'], 'wind_speed')
                rpg_bin.data['wind_speed'] = rpg_bin.data['wind_speed'] / 3.6
            if (int(rpg_met.header['_n_sen'],2) & 2) != 0:
                _add_interpol1(rpg_bin.data, rpg_met.data['adds'][:,1], rpg_met.data['time'], 'wind_direction')
            if (int(rpg_met.header['_n_sen'],2) & 4) != 0:
                _add_interpol1(rpg_bin.data, rpg_met.data['adds'][:,2], rpg_met.data['time'], 'rain_rate')
        
    elif data_type == '1B11':
        file_list_irt = get_file_list(path_to_files, path_to_prev, path_to_next, 'irt') 
        rpg_bin = get_rpg_bin(file_list_irt)
        rpg_bin.data['ir_wavelength'] = rpg_bin.header['_f']
        rpg_bin.data['ir_bandwidth'] = params['ir_bandwidth']
        rpg_bin.data['ir_beamwidth'] = params['ir_beamwidth']

    elif data_type == '1B21':
        file_list_met = get_file_list(path_to_files, path_to_prev, path_to_next, 'met')
        rpg_bin = get_rpg_bin(file_list_met)
        if (int(rpg_bin.header['_n_sen'],2) & 1) != 0:
            rpg_bin.data['wind_speed'] = rpg_bin.data['adds'][:,0] / 3.6
        if (int(rpg_bin.header['_n_sen'],2) & 2) != 0:
            rpg_bin.data['wind_direction'] = rpg_bin.data['adds'][:,1]
        if (int(rpg_bin.header['_n_sen'],2) & 4) != 0:
            rpg_bin.data['rain_rate'] = rpg_bin.data['adds'][:,2]                     
        
    else:
        raise RuntimeError(['Data type '+ data_type +' not supported for file writing.'])
    
    file_list_hkd = get_file_list(path_to_files, path_to_prev, path_to_next, 'hkd')
    _append_hkd(file_list_hkd, rpg_bin, data_type, params)
    rpg_bin.data['station_altitude'] = np.ones(len(rpg_bin.data['time']), np.float32) * params['station_altitude']
        
    return rpg_bin
    
    
def _append_hkd(file_list_hkd: list,
                rpg_bin: dict, 
                data_type: str,
                params: dict) -> None:
    """Append hkd data on same time grid and perform TB sanity check"""
    
    hkd = get_rpg_bin(file_list_hkd)
    
    if all(hkd.data['station_latitude'] == Fill_Value_Float):
        _add_interpol1(rpg_bin.data, np.ones(len(hkd.data['time'])) * params['station_latitude'], hkd.data['time'], 'station_latitude')                       
    else:
        idx = np.where(hkd.data['station_latitude'] != Fill_Value_Float)[0]
        _add_interpol1(rpg_bin.data, hkd.data['station_latitude'][idx], hkd.data['time'][idx], 'station_latitude')                       
    if all(hkd.data['station_longitude'] == Fill_Value_Float):           
        _add_interpol1(rpg_bin.data, np.ones(len(hkd.data['time'])) * params['station_longitude'], hkd.data['time'], 'station_longitude')        
    else:
        idx = np.where(hkd.data['station_longitude'] != Fill_Value_Float)[0]
        _add_interpol1(rpg_bin.data, hkd.data['station_longitude'][idx], hkd.data['time'][idx], 'station_longitude')   
    
    if data_type in ('1B01','1C01'):
        _add_interpol1(rpg_bin.data, np.mean(hkd.data['temp'][:,0:2], axis=1), hkd.data['time'], 't_amb')
        _add_interpol1(rpg_bin.data, hkd.data['temp'][:,2:4], hkd.data['time'], 't_rec')
    
        """check time intervals of +-15 min of .HKD data for sanity checks"""
        rpg_bin.data['status'] = np.zeros([len(rpg_bin.data['time']), len(rpg_bin.data['tb'][:].T)], np.int32)                
        for i_time, v_time in enumerate(rpg_bin.data['time']):
            ind = np.where((hkd.data['time'] >= v_time - 900) & (hkd.data['time'] <= v_time + 900))
            status = hkd.data['status'][ind]
            for bit in range(7):
                if np.any((status & 2**bit) == 0):
                    rpg_bin.data['status'][i_time, bit] = 1
                if np.any((status & 2**bit + 8) == 0):
                    rpg_bin.data['status'][i_time, bit + 7] = 1  
            if np.any(((status & 2**25) == 1) | ((status & 2**29) == 1)):
                rpg_bin.data['status'][i_time, 0:6] = 1
            if np.any(((status & 2**27) == 1) | ((status & 2**29) == 1)):
                rpg_bin.data['status'][i_time, 7:13] = 1              
                
        
def _add_interpol1(data0: dict, 
                   data1: np.ndarray, 
                   time1: np.ndarray,
                   output_name: str) -> None:
    
    if data1.ndim > 1:
        data0[output_name] = np.ones([len(data0['time']), data1.shape[1]], np.float32) * Fill_Value_Float
        for ndim in range(data1.shape[1]):
            data0[output_name][:,ndim] = np.interp(data0['time'], time1, data1[:,ndim])
    else:
        data0[output_name] = np.interp(data0['time'], time1, data1)
        
        
def add_time_bounds(time: np.ndarray,
                     int_time: int) -> np.ndarray:
    
    time_bounds = np.ones([len(time), 2]) * Fill_Value_Int
    time_bounds[:,0] = time - int_time
    time_bounds[:,1] = time 
    
    return time_bounds        
        
                
def _add_blb(brt: dict,
             blb: dict,
             hkd: dict,
             params: dict) -> None:
    """Add boundary-layer scans using a linear time axis"""
    
    xx = 0
    time_add = np.ones( blb.header['n'] * blb.header['_n_ang'], np.int32) * Fill_Value_Int
    ele_add = np.ones( blb.header['n'] * blb.header['_n_ang'], np.float32) * Fill_Value_Float
    azi_add = np.ones( blb.header['n'] * blb.header['_n_ang'], np.float32) * Fill_Value_Float
    tb_add = np.ones( [blb.header['n'] * blb.header['_n_ang'], blb.header['_n_f']], np.float32) * Fill_Value_Float
    rain_add = np.ones( blb.header['n'] * blb.header['_n_ang'], np.int32) * Fill_Value_Int
    
    bl_mod = np.ones(len(hkd.data['time']))*-99
    mul_ang = np.where(hkd.data['status'][:] & (1<<18))
    bl_mod[mul_ang] = 1      

    for time_ind, time_blb in enumerate(blb.data['time']):
        
        ind_scan = np.where((hkd.data['time'] >= time_blb - 1.5 * params['scan_time']) & (hkd.data['time'] <= time_blb) & (bl_mod == 1))[0]

        if np.any(ind_scan):
            
            if (not isbit(blb.data['rf_mod'][time_ind], 5)) & (not isbit(blb.data['rf_mod'][time_ind], 6)):
                sq = 0.
            elif (not isbit(blb.data['rf_mod'][time_ind], 5)) & (not isbit(blb.data['rf_mod'][time_ind], 6)):
                sq = 180.
            else:
                sq = 0.
            
            time_add[xx:xx + blb.header['_n_ang'] - 1] = np.linspace(hkd.data['time'][ind_scan[0]], hkd.data['time'][ind_scan[-2]], blb.header['_n_ang'] - 1)
            time_add[xx + blb.header['_n_ang'] - 1] = hkd.data['time'][ind_scan[-1]]
            azi_add[xx:xx + blb.header['_n_ang']] = (sq + params['const_azi']) % 360
            rain_add[xx:xx + blb.header['_n_ang']] = blb.data['rf_mod'][time_ind] & 1
            for ang in range(blb.header['_n_ang']):              
                tb_add[xx, :] = np.squeeze(blb.data['tb'][time_ind, :, ang])
                ele_add[xx] = blb.header['_ang'][ang]
                xx += 1

    bnd_add = add_time_bounds(time_add, params['scan_time'] / blb.header['_n_ang'] - 1)
    
    brt.data['time'] = np.concatenate((brt.data['time'], time_add))    
    ind = np.argsort(brt.data['time'])
    brt.data['time'] = brt.data['time'][ind]
    brt.data['time_bnds'] = np.concatenate((brt.data['time_bnds'], bnd_add))
    brt.data['time_bnds'] = brt.data['time_bnds'][ind,:]
    brt.data['ele'] = np.concatenate((brt.data['ele'], ele_add))
    brt.data['ele'] = brt.data['ele'][ind]
    brt.data['azi'] = np.concatenate((brt.data['azi'], azi_add))
    brt.data['azi'] = brt.data['azi'][ind]
    brt.data['tb'] = np.concatenate((brt.data['tb'], tb_add))
    brt.data['tb'] = brt.data['tb'][ind,:]
    brt.data['rain'] = np.concatenate((brt.data['rain'], rain_add))
    brt.data['rain'] = brt.data['rain'][ind]
    brt.data['pointing_flag'] = np.concatenate((brt.data['pointing_flag'], np.ones(len(time_add), np.int32)))
    brt.data['pointing_flag'] = brt.data['pointing_flag'][ind]
    brt.header['n'] = brt.header['n'] + blb.header['n'] * blb.header['_n_ang']
    
    
def _azi_correction(brt: dict, 
                    params: dict) -> None:
    """Azimuth correction (transform to "geographical" coordinates)"""
    
    ind180 = np.where((brt['azi'][:] >= 0) & (brt['azi'][:] <= 180))
    ind360 = np.where((brt['azi'][:] > 180) & (brt['azi'][:] <= 360))
    brt['azi'][ind180] = params['azi_cor'] - brt['azi'][ind180]
    brt['azi'][ind360] = 360. + params['azi_cor'] - brt['azi'][ind360]    
    brt['azi'][brt['azi'][:] < 0] += 360.