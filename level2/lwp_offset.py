import numpy as np
import pandas as pd
from pandas.tseries.frequencies import to_offset
from utils import df_interp

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
    

    freq_31 = np.where(lev1['frequency'][:] == 31.4)[0]
    if len(freq_31) != 1:
        lwp_offset = np.ones(len(index)) * Fill_Value_Float
    else:
        time = lev1['time'][index]
        tb = np.squeeze(lev1['tb'][index, freq_31])    
        lwp = np.copy(lwp_org)

        tb_df = pd.DataFrame({'Tb': tb}, index = pd.to_datetime(time, unit = 's'))
        tb_cc = tb_df.resample("2min", origin = 'start', closed = 'left', label = 'left').count()
        tb_cc.index = tb_cc.index + to_offset('1min')
        tb_std = tb_df.resample("2min", origin = 'start', closed = 'left', label = 'left').std()
        tb_std.index = tb_std.index + to_offset('1min')
        tb_std[tb_cc['Tb'] < 30] = np.nan    
        tb_mx = tb_std.resample("20min", origin = 'start', closed = 'left', label = 'left').max()
        tb_mx.index = tb_mx.index + to_offset('10min')
        tb_mx = df_interp(tb_mx, tb_df.index)      

        if 'irt' in lev1.variables:
            irt = lev1['irt'][index, 0]
            irt_df = pd.DataFrame({'Irt': irt[:]}, index = pd.to_datetime(time, unit = 's'))
            irt_mx = irt_df.resample("20min", origin = 'start', closed = 'left', label = 'left').max()
            irt_mx.index = irt_mx.index + to_offset('10min')
            irt_mx = df_interp(irt_mx, irt_df.index) 
            lwp[(irt_mx['Irt'] > 233.15) & (tb_mx['Tb'] > .2) | np.isnan(irt_mx['Irt']) | np.isnan(tb_mx['Tb']) | (lwp > .1)] = np.nan
        else:
            lwp[(tb_mx['Tb'] > .2) | np.isnan(tb_mx['Tb']) | (lwp > .1)] = np.nan

        lwp_df = pd.DataFrame({'Lwp': lwp}, index = pd.to_datetime(time, unit = 's'))
        lwp_mn = lwp_df.resample("20min", origin = 'start', closed = 'left', label = 'left').mean()
        lwp_mn.index = lwp_mn.index + to_offset('10min')

        off = pd.read_csv('site_config/' + site + '/lwp_offset.csv', usecols = ['date', 'offset'])
        ind = np.where(lwp_mn['Lwp'].values > 0)[0]
        if ind.size > 1:
            off = off.append(pd.DataFrame({'date': time[ind[0]], 'offset': lwp_mn['Lwp'][ind[0]]}, index = {0}), ignore_index = True)
            off = off.append(pd.DataFrame({'date': time[ind[-1]], 'offset': lwp_mn['Lwp'][ind[-1]]}, index = {0}), ignore_index = True)
        elif ind.size == 1:
            off = off.append(pd.DataFrame({'date': time[ind], 'offset': lwp_mn['Lwp'][int(ind)]}, index = {0}), ignore_index = True)
        off = off.sort_values(by=['date'])
        off = off.drop_duplicates(subset=['date'])
        off.to_csv('site_config/' + site + '/lwp_offset.csv', index = False)

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
