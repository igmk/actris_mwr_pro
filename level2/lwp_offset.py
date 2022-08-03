import os
import numpy as np
import pandas as pd
from time import gmtime
from pandas.tseries.frequencies import to_offset
from utils import df_interp
from level1.write_lev1_nc import find_lwcl_free

Fill_Value_Float = -999.


def correct_lwp_offset(lev1: dict,
                       lwp_org: np.ndarray,
                       index: np.ndarray, 
                       site: str) -> np.ndarray:
    """This function corrects Lwp offset using the
    2min standard deviation of the 31.4 GHz channel and IR temperature
    
    Args:
        lev1: Level 1 data.
        lwp: Lwp array.
        index: Index to use.
        site: site: Name of site.
        
    Examples:
        >>> from level2.lwp_offset import correct_lwp_offset
        >>> correct_lwp_offset(lev1, lwp, index, 'site_name')
    """     
    
     
    lwcl_i, _ = find_lwcl_free(lev1, index)      
    lwp = np.copy(lwp_org)
    lwp[(lwcl_i != 0) | (lwp > .1)] = np.nan
    time = lev1['time'][index] 
    lwp_df = pd.DataFrame({'Lwp': lwp}, index = pd.to_datetime(time, unit = 's'))
    lwp_mn = lwp_df.resample("20min", origin = 'start', closed = 'left', label = 'left').mean()
    lwp_mn.index = lwp_mn.index + to_offset('10min')

    t1=gmtime(time.data[0])
    if not os.path.isfile('site_config/' + site + '/lwp_offset_' + str(t1[0]) + '.csv'):
        df = pd.DataFrame({'date': [], 'offset': []})
        df.to_csv('site_config/' + site + '/lwp_offset_' + str(t1[0]) + '.csv')

    off = pd.read_csv('site_config/' + site + '/lwp_offset_' + str(t1[0]) + '.csv', usecols = ['date', 'offset'])
    ind = np.where(lwp_mn['Lwp'].values > 0)[0]
    if ind.size > 1:
        off = off.append(pd.DataFrame({'date': time[ind[0]], 'offset': lwp_mn['Lwp'][ind[0]]}, index = {0}), ignore_index = True)
        off = off.append(pd.DataFrame({'date': time[ind[-1]], 'offset': lwp_mn['Lwp'][ind[-1]]}, index = {0}), ignore_index = True)
    elif ind.size == 1:
        off = off.append(pd.DataFrame({'date': time[ind], 'offset': lwp_mn['Lwp'][int(ind)]}, index = {0}), ignore_index = True)
    off = off.sort_values(by=['date'])
    off = off.drop_duplicates(subset=['date'])
    off.to_csv('site_config/' + site + '/lwp_offset_' + str(t1[0]) + '.csv', index = False)

    off_ind = np.where((off['date'].values < time[0]) & (time[0] - off['date'].values < 48.*3600.))[0]
    if off_ind.size == 1:
        off_ind = np.array([int(off_ind), int(off_ind)])
    if (off_ind.size > 1) & (np.isnan(lwp_mn['Lwp'][0])):
        lwp_mn['Lwp'][0] = off['offset'][off_ind[-1]]
    off_ind = np.where((off['date'].values > time[-1]) & (off['date'].values - time[-1] < 48.*3600.))[0]
    if off_ind.size == 1:
        off_ind = np.array([int(off_ind), int(off_ind)])
    if (off_ind.size > 1) & (np.isnan(lwp_mn['Lwp'][-1])):
        lwp_mn['Lwp'][-1] = off['offset'][off_ind[0]]        

    lwp_mn = df_interp(lwp_mn, lwp_df.index)
    lwp_mn = lwp_mn.interpolate(method = 'linear')
    lwp_mn = lwp_mn.fillna(method = 'bfill')
    lwp_offset = lwp_mn['Lwp'].values
    lwp_offset[np.isnan(lwp_offset)] = 0
    lwp_org -= lwp_offset
    
    return lwp_org, lwp_offset
