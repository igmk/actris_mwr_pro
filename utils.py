""" This module contains general helper functions. """

from datetime import datetime, timezone
import time
import pytz
from numpy import ma
import numpy as np
import pandas as pd
from typing import Optional, Tuple, Union
import glob
import netCDF4
from scipy import signal, stats

SECONDS_PER_MINUTE = 60
SECONDS_PER_HOUR = 3600
SECONDS_PER_DAY = 86400


def seconds2hours(time_in_seconds: np.ndarray) -> np.ndarray:
    """Converts seconds since some epoch to fraction hour.
    Args:
        time_in_seconds: 1-D array of seconds since some epoch that starts on midnight.
    Returns:
        Time as fraction hour.
    Notes:
        Excludes leap seconds.
    """
    seconds_since_midnight = np.mod(time_in_seconds, SECONDS_PER_DAY)
    fraction_hour = seconds_since_midnight/SECONDS_PER_HOUR
    if fraction_hour[-1] == 0:
        fraction_hour[-1] = 24
    return fraction_hour


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


def isbit(array: np.ndarray, nth_bit: int) -> np.ndarray:
    """Tests if nth bit (0,1,2..) is set.
    Args:
        array: Integer array.
        nth_bit: Investigated bit.
    Returns:
        Boolean array denoting values where nth_bit is set.
    Raises:
        ValueError: negative bit as input.
    Examples:
        >>> isbit(4, 1)
            False
        >>> isbit(4, 2)
            True
    See also:
        utils.setbit()
    """
    if nth_bit < 0:
        raise ValueError('Negative bit number')
    mask = 1 << nth_bit
    return array & mask > 0


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


def binvec(x: Union[np.ndarray, list]) -> np.ndarray:
    """Converts 1-D center points to bins with even spacing.
    Args:
        x: 1-D array of N real values.
    Returns:
        ndarray: N + 1 edge values.
    Examples:
        >>> binvec([1, 2, 3])
            [0.5, 1.5, 2.5, 3.5]
    """
    edge1 = x[0] - (x[1] - x[0]) / 2
    edge2 = x[-1] + (x[-1] - x[-2]) / 2
    return np.linspace(edge1, edge2, len(x) + 1)


def rebin_2d(
    x_in: np.ndarray,
    array: ma.MaskedArray,
    x_new: np.ndarray,
    statistic: str = "mean",
    n_min: int = 1,
) -> Tuple[ma.MaskedArray, list]:
    """Rebins 2-D data in one dimension.
    Args:
        x_in: 1-D array with shape (n,).
        array: 2-D input data with shape (n, m).
        x_new: 1-D target vector (center points) with shape (N,).
        statistic: Statistic to be calculated. Possible statistics are 'mean', 'std'.
            Default is 'mean'.
        n_min: Minimum number of points to have good statistics in a bin. Default is 1.
    Returns:
        tuple: Rebinned data with shape (N, m) and indices of bins without enough data.
    Notes:
        0-values are masked in the returned array.
    """
    edges = binvec(x_new)
    result = np.zeros((len(x_new), array.shape[1]))
    array_screened = ma.masked_invalid(array, copy=True)  # data may contain nan-values
    for ind, values in enumerate(array_screened.T):
        mask = ~values.mask
        if ma.any(values[mask]):
            result[:, ind], _, _ = stats.binned_statistic(
                x_in[mask], values[mask], statistic=statistic, bins=edges
            )
    result[~np.isfinite(result)] = 0
    masked_result = ma.masked_equal(result, 0)

    # Fill bins with not enough profiles
    empty_indices = []
    for ind in range(len(edges) - 1):
        is_data = np.where((x_in > edges[ind]) & (x_in <= edges[ind + 1]))[0]
        if len(is_data) < n_min:
            masked_result[ind, :] = ma.masked
            empty_indices.append(ind)
    if len(empty_indices) > 0:
        logging.info(f"No data in {len(empty_indices)} bins")

    return masked_result, empty_indices


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
        raise RuntimeError(['Error: no coefficient files found in directory ' + 'site_config/' + site + '/coefficients/'])
    return c_list


def read_nc_field_name(nc_file: str, name: str) -> str:
    """Reads selected variable name from a netCDF file.
    Args:
        nc_file: netCDF file name.
        name: Variable to be read, e.g. 'temperature'.
    Returns:
        str
    """
    with netCDF4.Dataset(nc_file) as nc:
        long_name = nc.variables[name].getncattr('long_name')
    return long_name    


def read_nc_fields(nc_file: str, names: Union[str, list]) -> Union[ma.MaskedArray, list]:
    """Reads selected variables from a netCDF file.
    Args:
        nc_file: netCDF file name.
        names: Variables to be read, e.g. 'temperature' or ['ldr', 'lwp'].
    Returns:
        ndarray/list: Array in case of one variable passed as a string.
        List of arrays otherwise.
    """
    names = [names] if isinstance(names, str) else names
    with netCDF4.Dataset(nc_file) as nc:
        data = [nc.variables[name][:] for name in names]
    return data[0] if len(data) == 1 else data


def convolve2DFFT(slab,kernel,max_missing=0.1,verbose=True):
    """2D convolution using fft.
    <slab>: 2d array, with optional mask.
    <kernel>: 2d array, convolution kernel.
    <max_missing>: real, max tolerable percentage of missings within any
                   convolution window.
                   E.g. if <max_missing> is 0.5, when over 50% of values
                   within a given element are missing, the center will be
                   set as missing (<res>=0, <resmask>=1). If only 40% is
                   missing, center value will be computed using the remaining
                   60% data in the element.
                   NOTE that out-of-bound grids are counted as missings, this
                   is different from convolve2D(), where the number of valid
                   values at edges drops as the kernel approaches the edge.
    Return <result>: 2d convolution.
    """
    assert np.ndim(slab)==2, "<slab> needs to be 2D."
    assert np.ndim(kernel)==2, "<kernel> needs to be 2D."
    assert kernel.shape[0]<=slab.shape[0], "<kernel> size needs to <= <slab> size."
    assert kernel.shape[1]<=slab.shape[1], "<kernel> size needs to <= <slab> size."
    #--------------Get mask for missings--------------
    slab[slab == 0.] = np.nan
    slabcount=1-np.isnan(slab)
    # this is to set np.nan to a float, this won't affect the result as
    # masked values are not used in convolution. Otherwise, nans will
    # affect convolution in the same way as scipy.signal.convolve()
    # and the result will contain nans.
    slab=np.where(slabcount==1,slab,0)
    kernelcount=np.where(kernel==0,0,1)
    result=signal.fftconvolve(slab,kernel,mode='same')
    result_mask=signal.fftconvolve(slabcount,kernelcount,mode='same')
    valid_threshold=(1.-max_missing)*np.sum(kernelcount)
    result/=np.sum(kernel)
    result[(result_mask<valid_threshold)] = np.nan

    return result