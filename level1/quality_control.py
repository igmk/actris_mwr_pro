import numpy as np
import pandas as pd
from utils import setbit, get_coeff_list, df_interp
import datetime
import ephem
import netCDF4 as nc
from pandas.tseries.frequencies import to_offset

Fill_Value_Float = -999.
    
def apply_qc(data: dict, 
             params: dict) -> None: 
    """ This function performs the quality control of level 1 data.
    Args:
        data: Level 1 data.
        params: Site specific parameters.
        
    Returns:
        None
      
    Raises:
        RuntimeError: 
    
    Example:
        from level1.quality_control import apply_qc
        apply_qc('lev1_data','params')
       
    """    

    data['quality_flag'] = np.zeros(data['tb'].shape, dtype = np.int32)
    c_list = get_coeff_list(params['path_spec'], params['algo_spec'][:])
    
    for freq, _ in enumerate(data['frequency']):

        """ Bit 1: Missing TB-value """
        ind = np.where(data['tb'][:, freq] == Fill_Value_Float)
        data['quality_flag'][ind, freq] = setbit(data['quality_flag'][ind, freq], 0)
        
        """ Bit 2: TB threshold (lower range) """
        ind = np.where(data['tb'][:, freq] < params['TB_threshold'][0])
        data['quality_flag'][ind, freq] = setbit(data['quality_flag'][ind, freq], 1)  
        
        """ Bit 3: TB threshold (upper range) """
        ind = np.where(data['tb'][:, freq] > params['TB_threshold'][1])
        data['quality_flag'][ind, freq] = setbit(data['quality_flag'][ind, freq], 2)   
        
        """ Bit 4: Spectral consistency threshold """
        if any(c_list):
            ind = spectral_consistency(data, c_list[freq], freq, params['threshold_spec'][freq], params['factor_spec'][freq])
            data['quality_flag'][ind, freq] = setbit(data['quality_flag'][ind, freq], 3) 
        
        """ Bit 5: Receiver sanity """                
        ind = np.where(data['status'][:, freq] == 1)
        data['quality_flag'][ind, freq] = setbit(data['quality_flag'][ind, freq], 4)
        
        """ Bit 6: Rain flag """
        ind = np.where(data['rain'] == 1)
        data['quality_flag'][ind, freq] = setbit(data['quality_flag'][ind, freq], 5)
        
        """ Bit 7: Solar/Lunar flag """
        sun, moon = orbpos(data)
        ind = np.where((data['ele'][:] <= np.max(sun['ele']) + 10.) & (data['time'][:] >= sun['sunrise']) & (data['time'][:] <= sun['sunset']) & (data['ele'][:] >= sun['ele'][:] - params['saf']) & (data['ele'][:] <= sun['ele'][:] + params['saf']) & (data['azi'][:] >= sun['azi'][:] - params['saf']) & (data['azi'][:] <= sun['azi'][:] + params['saf']))
        data['quality_flag'][ind, freq] = setbit(data['quality_flag'][ind, freq], 6)
        
        """ Bit 8: TB offset threshold """
        
        
        
def orbpos(data: dict) -> dict:
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
        sun['ele'][ind] = np.rad2deg(sol.alt + 0.0)
        sun['azi'][ind] = np.rad2deg(sol.az + 0.0)    
        
        lun.compute(location)
        moon['ele'][ind] = np.rad2deg(lun.alt + 0.0)
        moon['azi'][ind] = np.rad2deg(lun.az + 0.0)         
    
    sun['sunrise'] = data['time'][0] + 0.
    sun['sunset'] = data['time'][0] + 24. * 3600.
    i_sun = np.where(sun['ele'] > 0.)
    if i_sun:
        sun['sunrise'] = data['time'][i_sun[0][0]]
        sun['sunset'] = data['time'][i_sun[-1][-1]]
    
    return sun, moon


def spectral_consistency(data: dict, 
                         c_file: str,
                         ind: np.int32,
                         threshold: np.float32,
                         factor: np.float32) -> np.ndarray:
    """ Applies spectral consistency coefficients for given frequency index and returns indices to be flagged """
    
    coeff = nc.Dataset(c_file)
    _, freq_ind, coeff_ind = np.intersect1d(data['frequency'], coeff['freq'], assume_unique = False, return_indices = True)
    ele_ind = np.squeeze(np.where((data['ele'][:] > coeff.variables['elevation_predictand'][:].data - .6) & (data['ele'][:] < coeff.variables['elevation_predictand'][:].data + .6)))
    if (ele_ind.size > 0) & (freq_ind.size > 0):
        
        tb_ret = coeff.variables['offset_mvr'][:].data + np.sum(coeff.variables['coefficient_mvr'][coeff_ind].T * data['tb'][:, freq_ind], axis = 1) + np.sum(coeff.variables['coefficient_mvr'][coeff_ind + (len(data['frequency']) - 1)].T * data['tb'][:, freq_ind]**2, axis = 1)
        tb_df = pd.DataFrame({'Tb': np.abs(data['tb'][ele_ind, ind]-tb_ret[ele_ind])}, index = pd.to_datetime(data['time'][ele_ind], unit = 's'))
        tb_std = tb_df.resample("2min", origin = 'start', closed = 'left', label = 'left').std()
        
        loffset1 = '1min'
        tb_std.index = tb_std.index + to_offset(loffset1)  
        org = pd.DataFrame({'Tb': tb_ret}, index = pd.to_datetime(data['time'][:], unit = 's'))
        tb_std = df_interp(tb_std, org.index)
        abs_diff = df_interp(tb_df, org.index)
        
        ind_flag = np.ones(len(data['time'][:])) * np.nan
        ind_flag[((data['ele'][:] > coeff.variables['elevation_predictand'][:].data - .6) & (data['ele'][:] < coeff.variables['elevation_predictand'][:].data + .6) & ((tb_std['Tb'].values > coeff.variables['predictand_err'][:].data * factor) | (abs_diff['Tb'].values > threshold)))] = 1
        df = pd.DataFrame({'Flag': ind_flag}, index = pd.to_datetime(data['time'][:], unit = 's'))
        df = df.fillna(method = 'bfill', limit = 120)
        df = df.fillna(method = 'ffill', limit = 300)    

        return np.squeeze(np.where(df['Flag'].values == 1))
    else:
        return []
    