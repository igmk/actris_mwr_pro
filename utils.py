""" This module contains general helper functions. """

from datetime import datetime, timezone
import time
import pytz
import numpy.ma as ma
import numpy as np
import pandas as pd
from typing import Optional, Tuple
import glob


def epoch2unix(epoch_time, time_ref, 
                 epoch: Optional[tuple] = (2001, 1, 1)):    
    """Converts seconds since (2001,1,1,0,0,0) to unix time in UTC.
    
    Args:
        epoch_time (ndarray): 1-D array of seconds since (2001,1,1,0,0,0)
        
    Returns:
        ndarray: Unix time in seconds since (1970,1,1,0,0,0).
        
    """

    delta = (datetime(*epoch) - datetime(1970,1,1,0,0,0)).total_seconds()
    unix_time = epoch_time + int(delta)
    if time_ref == 0:
        for index in range(len(unix_time)):
            unix_time[index] = time.mktime(datetime.fromtimestamp(unix_time[index], tz=timezone.utc).timetuple())
    return unix_time


def isscalar(array: any) -> bool:
    """Tests if input is scalar.
    By "scalar" we mean that array has a single value.
    Examples:
        >>> isscalar(1)
            True
        >>> isscalar([1])
            True
        >>> isscalar(np.array(1))
            True
        >>> isscalar(np.array([1]))
            True
    """
    
    arr = ma.array(array)
    if not hasattr(arr, '__len__') or arr.shape == () or len(arr) == 1:
        return True
    return False


def setbit(array: np.ndarray, 
           nth_bit: int) -> int:
    """Sets nth bit (0, 1, 2..) on number.
    Args:
        array: Integer array.
        nth_bit: Bit to be set.
    Returns:
        Integer where nth bit is set.
    Raises:
        ValueError: negative bit as input.
    Examples:
        >>> setbit(1, 1)
            3
        >>> setbit(0, 2)
            4
    See also:
        utils.isbit()
    """
    
    if nth_bit < 0:
        raise ValueError('Negative bit number')
    mask = 1 << nth_bit
    array |= mask
    return array


def df_interp(df, new_index):
    """Return a new DataFrame with all columns values interpolated
    to the new_index values."""
    
    df_out = pd.DataFrame(index=new_index)
    df_out.index.name = df.index.name

    for colname, col in df.iteritems():
        df_out[colname] = np.interp(new_index, df.index, col)

    return df_out


def seconds2date(time_in_seconds: float, 
                 epoch: Optional[tuple] = (1970, 1, 1)) -> list:
    """Converts seconds since some epoch to datetime (UTC).
    Args:
        time_in_seconds: Seconds since some epoch.
        epoch: Epoch, default is (1970, 1, 1) (UTC).
    Returns:
        [year, month, day, hours, minutes, seconds] formatted as '05' etc (UTC).
    """
    
    epoch_in_seconds = datetime.timestamp(datetime(*epoch, tzinfo=pytz.utc))
    timestamp = time_in_seconds + epoch_in_seconds
    return datetime.utcfromtimestamp(timestamp).strftime('%Y %m %d %H %M %S').split()


def get_coeff_list(site: str,
                   prefix: str):
    "Returns list of .nc coefficient file(s)"
    
    c_list = sorted(glob.glob('site_config/' + site + '/coefficients/' + prefix + '*.nc'))
    if len(c_list) < 1:
        raise RuntimeError(['Error: no coefficient files found in directory ' + path_to_files])
    return c_list