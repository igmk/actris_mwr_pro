"""Module for general helper functions."""

import glob
import re
import time
from typing import Iterator
import yaml
from yaml.loader import SafeLoader
from datetime import date, datetime, timedelta, timezone

import netCDF4
import numpy as np
from numpy import ma
from scipy import signal

SECONDS_PER_MINUTE = 60
SECONDS_PER_HOUR = 3600
SECONDS_PER_DAY = 86400
Fill_Value_Float = -999.0
Fill_Value_Int = -99
Epoch = tuple[int, int, int]


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
    fraction_hour = seconds_since_midnight / SECONDS_PER_HOUR
    if fraction_hour[-1] == 0:
        fraction_hour[-1] = 24
    return fraction_hour


def epoch2unix(epoch_time, time_ref, epoch: Epoch = (2001, 1, 1)):
    """Converts seconds since (2001,1,1,0,0,0) to unix time in UTC.

    Args:
        epoch_time (ndarray): 1-D array of seconds since (2001,1,1,0,0,0)

    Returns:
        ndarray: Unix time in seconds since (1970,1,1,0,0,0).

    """

    delta = (datetime(*epoch) - datetime(1970, 1, 1, 0, 0, 0)).total_seconds()
    unix_time = epoch_time + int(delta)
    if time_ref == 0:
        for index, _ in enumerate(unix_time):
            unix_time[index] = time.mktime(
                datetime.fromtimestamp(unix_time[index], timezone.utc).timetuple()
            )
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
    if not hasattr(arr, "__len__") or arr.shape == () or len(arr) == 1:
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
        raise ValueError("Negative bit number")
    mask = 1 << nth_bit
    return array & mask > 0


def setbit(array: np.ndarray, nth_bit: int) -> np.ndarray:
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
        raise ValueError("Negative bit number")
    mask = 1 << nth_bit
    array |= mask
    return array


def interpol_2d(
    x_in: np.ndarray,
    array: ma.MaskedArray,
    x_new: np.ndarray,
) -> ma.MaskedArray:
    """Interpolates 2-D data in one dimension.
    Args:
        x_in: 1-D array with shape (n,).
        array: 2-D input data with shape (n, m).
        x_new: 1-D target vector with shape (N,).
    Returns:
        array: Interpolated data with shape (N, m).
    Notes:
        0-values are masked in the returned array.
    """
    result = np.zeros((len(x_new), array.shape[1]))
    array_screened = ma.masked_invalid(array, copy=True)  # data may contain nan-values
    for ind, values in enumerate(array_screened.T):
        if array.mask == True:
            mask = ~values.mask
            if ma.any(values[mask]):
                result[:, ind] = np.interp(x_new, x_in[mask], values[mask])
        else:
            result[:, ind] = np.interp(x_new, x_in, values)
    result[~np.isfinite(result)] = 0
    masked = ma.make_mask(result)

    return ma.array(result, mask=~masked)


def add_interpol1d(data0: dict, data1: np.ndarray, time1: np.ndarray, output_name: str) -> None:
    """Adds interpolated 1d field to dict
    Args:
        data0: Output dict.
        data1: Input field to be added & interpolated.
        time1: Time of input field.
    """
    if data1.ndim > 1:
        data0[output_name] = (
            np.ones([len(data0["time"]), data1.shape[1]], np.float32) * Fill_Value_Float
        )
        for ndim in range(data1.shape[1]):
            data0[output_name][:, ndim] = np.interp(data0["time"], time1, data1[:, ndim])
    else:
        data0[output_name] = np.interp(data0["time"], time1, data1)


def seconds2date(time_in_seconds: float, epoch: Epoch = (1970, 1, 1)) -> list:
    """Converts seconds since some epoch to datetime (UTC).
    Args:
        time_in_seconds: Seconds since some epoch.
        epoch: Epoch, default is (1970, 1, 1) (UTC).
    Returns:
        [year, month, day, hours, minutes, seconds] formatted as '05' etc (UTC).
    """

    epoch_in_seconds = datetime.timestamp(datetime(*epoch, tzinfo=timezone.utc))
    timestamp = time_in_seconds + epoch_in_seconds
    return datetime.utcfromtimestamp(timestamp).strftime("%Y %m %d %H %M %S").split()


def add_time_bounds(time_arr: np.ndarray, int_time: int) -> np.ndarray:
    "Adds time bounds"
    time_bounds = np.ones([len(time_arr), 2], np.int32) * Fill_Value_Int
    time_bounds[:, 0] = time_arr - int_time
    time_bounds[:, 1] = time_arr

    return time_bounds


def get_coeff_list(site: str, prefix: str):
    "Returns list of .nc coefficient file(s)"

    s_list = [
        glob.glob("site_config/" + site + "/coefficients/" + prefix.lower() + "*"),
        glob.glob("site_config/" + site + "/coefficients/" + prefix.upper() + "*"),
    ]
    c_list = [x for x in s_list if x != []]
    if len(c_list[0]) < 1:
        raise RuntimeError(
            [
                "Error: no coefficient files found in directory "
                + "site_config/"
                + site
                + "/coefficients/"
            ]
        )
    return sorted(c_list[0])


def get_file_list(path_to_files: str, path_to_prev: str, path_to_next: str, extension: str):
    """Returns file list for specified path."""
    f_list = sorted(glob.glob(path_to_files + "*." + extension))
    if len(f_list) == 0:
        f_list = sorted(glob.glob(path_to_files + "*." + extension.upper()))
    f_list_p = sorted(glob.glob(path_to_prev + "*." + extension))
    if len(f_list_p) == 0:
        f_list_p = sorted(glob.glob(path_to_prev + "*." + extension.upper()))
    f_list_n = sorted(glob.glob(path_to_next + "*." + extension))
    if len(f_list_n) == 0:
        f_list_n = sorted(glob.glob(path_to_next + "*." + extension.upper()))

    if len(f_list) == 0:
        raise RuntimeError(
            [
                "Error: no binary files with extension "
                + extension
                + " found in directory "
                + path_to_files
            ]
        )
    if (len(f_list_p) > 0) & (len(f_list_n) > 0):
        f_list = [*f_list_p, *f_list, *f_list_n]
    elif (len(f_list_p) > 0) & (len(f_list_n) == 0):
        f_list = [*f_list_p, *f_list]
    elif (len(f_list_p) == 0) & (len(f_list_n) > 0):
        f_list = [*f_list, *f_list_n]
    return f_list


def read_yaml_config(site: str) -> tuple[dict, dict]:
    """Reads config yaml files."""
    site_file = "site_config/" + site + "/config.yaml"
    with open(site_file) as f:
        site_config = yaml.load(f, Loader=SafeLoader)
    inst_file = "site_config/" + site_config["type"] + ".yaml"
    with open(inst_file) as f:
        inst_config = yaml.load(f, Loader=SafeLoader)
    inst_config["global_specs"].update(site_config["global_specs"])
    for name in inst_config["params"].keys():
        site_config["params"][name] = inst_config["params"][name]

    return inst_config["global_specs"], site_config["params"]


def update_lev1_attributes(attributes: dict, data_type: str) -> None:
    """Removes attributes that are not needed for specified Level 1 data type"""
    if data_type == "1B01":
        att_del = ["ir_instrument", "met_instrument", "_accuracy"]
        key = " "
    elif data_type == "1B11":
        att_del = [
            "instrument_manufacturer",
            "instrument_model",
            "instrument_generation",
            "instrument_hw_id",
            "instrument_calibration",
            "receiver",
            "date_of",
            "instrument_history",
            "met",
            "air",
            "relative",
            "wind",
            "rain",
        ]
        key = "ir_"
    elif data_type == "1B21":
        att_del = [
            "instrument_manufacturer",
            "instrument_model",
            "instrument_generation",
            "instrument_hw_id",
            "instrument_calibration",
            "receiver",
            "date_of",
            "instrument_history",
            "ir_instrument",
            "ir_accuracy",
        ]
        key = "met"
        attributes["source"] = "In Situ"

    for name in list(attributes.keys()):
        if any(x in name for x in att_del) & (name[0:3] != key):
            del attributes[name]


def get_ret_ang(nc_file: str) -> list:
    """Returns highest elevation angle used in retrieval."""
    with netCDF4.Dataset(nc_file) as nc:
        ang = []
        ang_str = re.split(r"[^0-9.]", nc.retrieval_elevation_angles)
        for ii in ang_str:
            if ii != "":
                ang.append(float(ii))
    return ang


def get_ret_freq(nc_file: str) -> np.ndarray:
    """Returns frequencies used in retrieval."""
    with netCDF4.Dataset(nc_file) as nc:
        freq = []
        freq_str = re.split(r"[^0-9.]", nc.retrieval_frequencies)
        for ii in freq_str:
            if ii != "":
                freq.append(float(ii))
    return np.array(freq, dtype=np.float32)


def read_nc_field_name(nc_file: str, name: str) -> str:
    """Reads selected variable name from a netCDF file.
    Args:
        nc_file: netCDF file name.
        name: Variable to be read, e.g. 'temperature'.
    Returns:
        str
    """
    with netCDF4.Dataset(nc_file) as nc:
        long_name = nc.variables[name].getncattr("long_name")
    return long_name


def read_nc_fields(nc_file: str, names: str | list) -> ma.MaskedArray | list:
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


def append_data(data_in: dict, key: str, array: ma.MaskedArray) -> dict:
    """Appends data to a dictionary field (creates the field if not yet present).
    Args:
        data_in: Dictionary where data will be appended.
        key: Key of the field.
        array: Numpy array to be appended to data_in[key].
    """
    data = data_in.copy()
    if key not in data:
        if array.ndim == 1:
            data[key] = array[:]
        else:
            data[key] = array[:, :]
    else:
        data[key] = ma.concatenate((data[key], array))
    return data


def convolve2DFFT(slab, kernel, max_missing=0.1):
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
    assert np.ndim(slab) == 2, "<slab> needs to be 2D."
    assert np.ndim(kernel) == 2, "<kernel> needs to be 2D."
    assert kernel.shape[0] <= slab.shape[0], "<kernel> size needs to <= <slab> size."
    assert kernel.shape[1] <= slab.shape[1], "<kernel> size needs to <= <slab> size."
    # --------------Get mask for missings--------------
    slab[slab == 0.0] = np.nan
    slabcount = 1 - np.isnan(slab)
    # this is to set np.nan to a float, this won't affect the result as
    # masked values are not used in convolution. Otherwise, nans will
    # affect convolution in the same way as scipy.signal.convolve()
    # and the result will contain nans.
    slab = np.where(slabcount == 1, slab, 0)
    kernelcount = np.where(kernel == 0, 0, 1)
    result = signal.fftconvolve(slab, kernel, mode="same")
    result_mask = signal.fftconvolve(slabcount, kernelcount, mode="same")
    valid_threshold = (1.0 - max_missing) * np.sum(kernelcount)
    result /= np.sum(kernel)
    result[(result_mask < valid_threshold)] = np.nan

    return result


def date_string_to_date(date_string: str) -> datetime.date:
    """Convert YYYY-MM-DD to Python date."""
    date_arr = [int(x) for x in date_string.split("-")]
    return date(*date_arr)


def get_time() -> str:
    """Returns current UTC-time."""
    return f"{datetime.utcnow().strftime('%Y-%m-%d %H:%M:%S')} +00:00"


def get_date_from_past(n: int, reference_date: str | None = None) -> str:
    """Return date N-days ago.
    Args:
        n: Number of days to skip (can be negative, when it means the future).
        reference_date: Date as "YYYY-MM-DD". Default is the current date.
    Returns:
        str: Date as "YYYY-MM-DD".
    """
    reference = reference_date or get_time().split()[0]
    date = date_string_to_date(reference) - timedelta(n)
    return str(date)


def get_processing_dates(args) -> tuple[str, str]:
    """Returns processing dates."""
    if args.date is not None:
        start_date = args.date
        stop_date = get_date_from_past(-1, start_date)
    else:
        start_date = args.start
        stop_date = args.stop
    start_date = str(date_string_to_date(start_date))
    stop_date = str(date_string_to_date(stop_date))
    return start_date, stop_date


def isodate2date(date_str: str) -> datetime.date:
    return datetime.strptime(date_str, "%Y-%m-%d").date()


def date_range(start_date: datetime.date, end_date: datetime.date) -> Iterator[datetime.date]:
    """Returns range between two dates (datetimes)."""
    for n in range(int((end_date - start_date).days)):
        yield start_date + timedelta(n)
