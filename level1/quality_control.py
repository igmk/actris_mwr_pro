import numpy as np
from numpy import ma
import pandas as pd
from utils import setbit, get_coeff_list, df_interp
import datetime
import ephem
import netCDF4 as nc
from pandas.tseries.frequencies import to_offset

Fill_Value_Float = -999.
Fill_Value_Int = -99

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
    data['tb'][data['tb'] == 0.] = Fill_Value_Float
    c_list = get_coeff_list(site, 'tbx')
    ind_bit4, _ = spectral_consistency(data, c_list, params['tbx_f'])
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
            ind = np.where(ind_bit4[:, freq] == 1)
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
    # moon = dict()
    # moon['azi'] = np.zeros(data['time'].shape) * Fill_Value_Float
    # moon['ele'] = np.zeros(data['time'].shape) * Fill_Value_Float

    sol = ephem.Sun()
    # lun = ephem.Moon()
    obs_loc = ephem.Observer()

    for ind, tim in enumerate(data['time']):
        
        obs_loc.lat, obs_loc.lon = str(data['station_latitude'][ind]), str(data['station_longitude'][ind])
        obs_loc.elevation = data['station_altitude'][ind]
        obs_loc.date = datetime.datetime.utcfromtimestamp(tim).strftime('%Y/%m/%d %H:%M:%S')
        sol.compute(obs_loc)
        sun['ele'][ind] = np.rad2deg(sol.alt)
        sun['azi'][ind] = np.rad2deg(sol.az)            
        # lun.compute(obs_loc)
        # moon['ele'][ind] = np.rad2deg(lun.alt)
        # moon['azi'][ind] = np.rad2deg(lun.az)         
    
    sun['sunrise'] = data['time'][0]
    sun['sunset'] = data['time'][0] + 24. * 3600.
    i_sun = np.where(sun['ele'] > 0.)

    if i_sun:
        sun['sunrise'] = data['time'][i_sun[0][0]]
        sun['sunset'] = data['time'][i_sun[-1][-1]]
        
    ind = np.where((data['ele'][:] <= np.max(sun['ele']) + 10.) & (data['time'][:] >= sun['sunrise']) & (data['time'][:] <= sun['sunset']) & (data['ele'][:] >= sun['ele'][:] - params['saf']) & (data['ele'][:] <= sun['ele'][:] + params['saf']) & (data['azi'][:] >= sun['azi'][:] - params['saf']) & (data['azi'][:] <= sun['azi'][:] + params['saf']))
    
    return ind


def spectral_consistency(data: dict, 
                         c_file: list,
                         tbx_f: np.ndarray) -> np.ndarray:
    """ Applies spectral consistency coefficients for given frequency index and returns indices to be flagged """
    
    flag_ind = np.zeros(data['tb'].shape)
    flag_tmp = np.ones(data['tb'].shape) * np.nan
    tb_tot = np.ones(data['tb'].shape) * np.nan
    tb_ret = np.ones(data['tb'].shape) * np.nan
    for ifreq, freq in enumerate(data['frequency']):
        with nc.Dataset(c_file[ifreq]) as coeff:
            _, freq_ind, coeff_ind = np.intersect1d(data['frequency'], coeff['freq'], assume_unique = False, return_indices = True)
            ele_ind = np.where((data['ele'][:] > coeff['elevation_predictand'][:] - .6) & (data['ele'][:] < coeff['elevation_predictand'][:] + .6))[0]

            if (ele_ind.size > 0) & (freq_ind.size > 0):
                tb_ret[:, ifreq] = coeff['offset_mvr'][:] + np.sum(coeff['coefficient_mvr'][coeff_ind].T * data['tb'][:, freq_ind], axis = 1) + np.sum(coeff['coefficient_mvr'][coeff_ind + (len(data['frequency']) - 1)].T * data['tb'][:, freq_ind]**2, axis = 1)
                tb_df = pd.DataFrame({'Tb': (data['tb'][ele_ind, ifreq]-tb_ret[ele_ind, ifreq])}, index = pd.to_datetime(data['time'][ele_ind], unit = 's'))
                org = pd.DataFrame({'Tb': tb_ret[:, ifreq]}, index = pd.to_datetime(data['time'][:], unit = 's'))
                tb_med = tb_df.resample("10min", origin = 'start', closed = 'left', label = 'left').mean()
                tb_med.index = tb_med.index + to_offset('5min')  
                tb_df = df_interp(tb_df, org.index)
                tb_med = df_interp(tb_med, org.index)                

                flag_tmp[((data['ele'][:] > coeff['elevation_predictand'][:] - .6) & (data['ele'][:] < coeff['elevation_predictand'][:] + .6) & (np.abs(tb_df['Tb'].values - tb_med['Tb'].values) > coeff['predictand_err'][:]*tbx_f[ifreq])), ifreq] = 1
                tb_tot[ele_ind, ifreq] = np.abs(tb_df['Tb'][ele_ind])
    
    for _, rec in enumerate(data['receiver_nb']):
        flag_tmp[np.ix_(np.sum(tb_tot[:, data['receiver'] == rec], axis=1)/np.nanmean(np.sum(tb_tot[:, data['receiver'] == rec], axis=1)) > 2., data['receiver'] == rec)] = 1   

    for ifreq, _ in enumerate(data['frequency']):
	    df = pd.DataFrame({'Flag': flag_tmp[:, ifreq]}, index = pd.to_datetime(data['time'][:], unit = 's'))
	    df = df.fillna(method = 'bfill', limit = 60)
	    df = df.fillna(method = 'ffill', limit = 60)    
	    flag_ind[((df['Flag'].values == 1) & (data['pointing_flag'] == 0)), ifreq] = 1

    return flag_ind, tb_ret
    
