import numpy as np
from utils import setbit
   
def apply_met_qc(data: dict, 
             params: dict) -> None: 
    """ This function performs the met quality control of level 1 data.
    Args:
        data: Level 1 data.
        params: Site specific parameters.
        
    Returns:
        None
      
    Raises:
        RuntimeError: 
    
    Example:
        from level1.met_quality_control import apply_met_qc
        apply_met_qc('lev1_data', 'params')
       
    """    

    data['met_quality_flag'] = np.zeros(len(data['time']), dtype = np.int32)
    
    """ Bit 1: low_quality_air_temperature """
    ind = np.where((data['air_temperature'][:] < params['met_thresholds'][0][0]) | (data['air_temperature'][:] > params['met_thresholds'][0][1]))
    data['met_quality_flag'][ind] = setbit(data['met_quality_flag'][ind], 0)

    """ Bit 2: low_quality_relative_humidity """
    ind = np.where((data['relative_humidity'][:] < params['met_thresholds'][1][0]) | (data['relative_humidity'][:] > params['met_thresholds'][1][1]))
    data['met_quality_flag'][ind] = setbit(data['met_quality_flag'][ind], 1)        

    """ Bit 3: low_quality_air_pressure """
    ind = np.where((data['air_pressure'][:] < params['met_thresholds'][2][0]) | (data['air_pressure'][:] > params['met_thresholds'][2][1]))
    data['met_quality_flag'][ind] = setbit(data['met_quality_flag'][ind], 2)            

    """ Bit 4: low_quality_rain_rate """
    ind = np.where((data['rain_rate'][:] < params['met_thresholds'][3][0]) | (data['rain_rate'][:] > params['met_thresholds'][3][1]))
    data['met_quality_flag'][ind] = setbit(data['met_quality_flag'][ind], 3)         

    """ Bit 5: low_quality_wind_direction """
    ind = np.where((data['wind_direction'][:] < params['met_thresholds'][4][0]) | (data['wind_direction'][:] > params['met_thresholds'][4][1]))
    data['met_quality_flag'][ind] = setbit(data['met_quality_flag'][ind], 4)          

    """ Bit 6: low_quality_wind_speed """
    ind = np.where((data['wind_speed'][:] < params['met_thresholds'][5][0]) | (data['wind_speed'][:] > params['met_thresholds'][5][1]))
    data['met_quality_flag'][ind] = setbit(data['met_quality_flag'][ind], 5)           