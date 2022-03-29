import numpy as np
import pandas as pd
from utils import setbit, get_coeff_list, df_interp
import datetime
import ephem
import netCDF4 as nc
from pandas.tseries.frequencies import to_offset
Fill_Value_Float = -999.
   
def apply_qc(site: str, 
             data: dict, 
             params: dict) -> None: 
    """ This function performs the quality control of level 1 data.
    Args:
        site: Name of site.
        data: Level 1 data.
        params: Site specific parameters.
        
    Returns:
        None
      
    Raises:
        RuntimeError: 
    
    Example:
        from level1.quality_control import apply_qc
        apply_qc('site', 'lev1_data', 'params')
       
    """    

    data['quality_flag'] = np.zeros(data['tb'].shape, dtype = np.int32)
    data['quality_flag_status'] = np.zeros(data['tb'].shape, dtype = np.int32)
    c_list = get_coeff_list(site, 'tbx')
    ind_bit6 = np.where(data['rain'] == 1)
    ind_bit7 = orbpos(data, params)
        
    for freq, _ in enumerate(data['frequency']):

        """ Bit 1: Missing TB-value """
        if params['flag_status'][0] == 1:
            data['quality_flag_status'][:, freq] = setbit(data['quality_flag_status'][:, freq], 0)
        else:
            ind = np.where(data['tb'][:, freq] == Fill_Value_Float)
            data['quality_flag'][ind, freq] = setbit(data['quality_flag'][ind, freq], 0)
        
        """ Bit 2: TB threshold (lower range) """
        if params['flag_status'][1] == 1:
            data['quality_flag_status'][:, freq] = setbit(data['quality_flag_status'][:, freq], 1)
        else:        
            ind = np.where(data['tb'][:, freq] < params['TB_threshold'][0])
            data['quality_flag'][ind, freq] = setbit(data['quality_flag'][ind, freq], 1)  
        
        """ Bit 3: TB threshold (upper range) """
        if params['flag_status'][2] == 1:
            data['quality_flag_status'][:, freq] = setbit(data['quality_flag_status'][:, freq], 2)
        else:        
            ind = np.where(data['tb'][:, freq] > params['TB_threshold'][1])
            data['quality_flag'][ind, freq] = setbit(data['quality_flag'][ind, freq], 2)   
        
        """ Bit 4: Spectral consistency threshold """
        if params['flag_status'][3] == 1:
            data['quality_flag_status'][:, freq] = setbit(data['quality_flag_status'][:, freq], 3)
        else:        
            ind = spectral_consistency(data, c_list[freq], freq, params['th_std'][freq])
            data['quality_flag'][ind, freq] = setbit(data['quality_flag'][ind, freq], 3) 
        
        """ Bit 5: Receiver sanity """        
        if params['flag_status'][4] == 1:
            data['quality_flag_status'][:, freq] = setbit(data['quality_flag_status'][:, freq], 4)
        else:        
            ind = np.where(data['status'][:, freq] == 1)
            data['quality_flag'][ind, freq] = setbit(data['quality_flag'][ind, freq], 4)
        
        """ Bit 6: Rain flag """
        if params['flag_status'][5] == 1:
            data['quality_flag_status'][:, freq] = setbit(data['quality_flag_status'][:, freq], 5)
        else:        
            data['quality_flag'][ind_bit6, freq] = setbit(data['quality_flag'][ind_bit6, freq], 5)
        
        """ Bit 7: Solar/Lunar flag """
        if params['flag_status'][6] == 1:
            data['quality_flag_status'][:, freq] = setbit(data['quality_flag_status'][:, freq], 6)
        else:        
            data['quality_flag'][ind_bit7, freq] = setbit(data['quality_flag'][ind_bit7, freq], 6)
        
        """ Bit 8: TB offset threshold """
        if params['flag_status'][7] == 1:
            data['quality_flag_status'][:, freq] = setbit(data['quality_flag_status'][:, freq], 7)
        # else:        
        
        
        
def orbpos(data: dict,
           params: dict) -> np.ndarray:
    """ Calculates sun & moon elevation/azimuth angles """
    
    sun = dict()
    sun['azi'] = np.zeros(data['time'].shape) * Fill_Value_Float
    sun['ele'] = np.zeros(data['time'].shape) * Fill_Value_Float
    moon = dict()
    moon['azi'] = np.zeros(data['time'].shape) * Fill_Value_Float
    moon['ele'] = np.zeros(data['time'].shape) * Fill_Value_Float

    sol = ephem.Sun()
    lun = ephem.Moon()
    location = ephem.Observer()

    for ind, _ in enumerate(data['time']):
       
        location.lat = str(data['station_latitude'][ind])
        location.lon = str(data['station_longitude'][ind])
        location.date = datetime.datetime.fromtimestamp(data['time'][ind]).strftime('%Y/%m/%d %H:%M:%S')
        sol.compute(location)
        sun['ele'][ind] = np.rad2deg(sol.alt)
        sun['azi'][ind] = np.rad2deg(sol.az)            
        lun.compute(location)
        moon['ele'][ind] = np.rad2deg(lun.alt)
        moon['azi'][ind] = np.rad2deg(lun.az)         
    
    sun['sunrise'] = data['time'][0]
    sun['sunset'] = data['time'][0] + 24. * 3600.
    i_sun = np.where(sun['ele'] > 0.)
    
    if i_sun:
        sun['sunrise'] = data['time'][i_sun[0][0]]
        sun['sunset'] = data['time'][i_sun[-1][-1]]
        
    ind = np.where((data['ele'][:] <= np.max(sun['ele']) + 10.) & (data['time'][:] >= sun['sunrise']) & (data['time'][:] <= sun['sunset']) & (data['ele'][:] >= sun['ele'][:] - params['saf']) & (data['ele'][:] <= sun['ele'][:] + params['saf']) & (data['azi'][:] >= sun['azi'][:] - params['saf']) & (data['azi'][:] <= sun['azi'][:] + params['saf']))
    
    return ind


def spectral_consistency(data: dict, 
                         c_file: str,
                         ind: np.int32,
                         th_std: np.float32) -> np.ndarray:
    """ Applies spectral consistency coefficients for given frequency index and returns indices to be flagged """
    
    coeff = nc.Dataset(c_file)
    _, freq_ind, coeff_ind = np.intersect1d(data['frequency'], coeff['freq'], assume_unique = False, return_indices = True)
    ele_ind = np.where((data['ele'][:] > coeff['elevation_predictand'][:] - .6) & (data['ele'][:] < coeff['elevation_predictand'][:] + .6))[0]
    flag_ind = []
    if (ele_ind.size > 0) & (freq_ind.size > 0):
        
        tb_ret = coeff['offset_mvr'][:] + np.sum(coeff['coefficient_mvr'][coeff_ind].T * data['tb'][:, freq_ind], axis = 1) + np.sum(coeff['coefficient_mvr'][coeff_ind + (len(data['frequency']) - 1)].T * data['tb'][:, freq_ind]**2, axis = 1)
        tb_df = pd.DataFrame({'Tb': np.abs(data['tb'][ele_ind, ind]-tb_ret[ele_ind])}, index = pd.to_datetime(data['time'][ele_ind], unit = 's'))
        
        tb_std = tb_df.resample("5min", origin = 'start', closed = 'left', label = 'left').std()
        tb_std.index = tb_std.index + to_offset('150s')  
        org = pd.DataFrame({'Tb': tb_ret}, index = pd.to_datetime(data['time'][:], unit = 's'))
        tb_std = df_interp(tb_std, org.index)
        
        ind_flag = np.ones(len(data['time'][:])) * np.nan
        ind_flag[((data['ele'][:] > coeff['elevation_predictand'][:] - .6) & (data['ele'][:] < coeff['elevation_predictand'][:] + .6) & ((tb_std['Tb'].values > th_std)))] = 1
        df = pd.DataFrame({'Flag': ind_flag}, index = pd.to_datetime(data['time'][:], unit = 's'))
        df = df.fillna(method = 'bfill', limit = 120)
        df = df.fillna(method = 'ffill', limit = 300)    
        flag_ind = np.where(df['Flag'].values == 1)[0]
        
    coeff.close()
    return flag_ind
    