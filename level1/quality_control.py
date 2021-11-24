import numpy as np
from utils import setbit
import datetime
import ephem

# import pdb
# pdb.set_trace()

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