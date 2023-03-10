"""Module containing additional helper functions for plotting"""
import locale
from datetime import datetime, timezone

import matplotlib.ticker as ticker
import netCDF4
import numpy as np
from numpy import ma, ndarray

from utils import (
    convolve2DFFT,
    get_ret_freq,
    isbit,
    read_nc_fields,
    read_yaml_config,
    seconds2hours,
)


def _get_ret_flag(nc_file: str, time: ndarray) -> ndarray:
    """Returns quality flag for frequencies used in retrieval."""
    flag = np.zeros(len(time), np.int32)
    lev1_file = _get_lev1(nc_file)
    quality_flag = read_nc_fields(lev1_file, "quality_flag")
    t_lev1 = read_nc_fields(lev1_file, "time")
    f_lev1 = read_nc_fields(lev1_file, "frequency")
    freq = get_ret_freq(nc_file)
    _, freq_ind, _ = np.intersect1d(f_lev1.data, freq, assume_unique=False, return_indices=True)
    _, q_ind, t_ind = np.intersect1d(
        seconds2hours(t_lev1.data), time.data, assume_unique=False, return_indices=True
    )
    quality_flag = quality_flag[q_ind, :]
    site = _read_location(nc_file)
    _, params = read_yaml_config(site)

    if params["flag_status"][3] == 0:
        flag[t_ind[np.sum(isbit(quality_flag[:, freq_ind], 3), axis=1) > 0]] = 1
    else:
        flag[t_ind[np.sum(quality_flag[:, freq_ind], axis=1) > 0]] = 1
    return flag


def _get_lev1(nc_file: str) -> str:
    """Returns name of lev1 file."""
    site = _read_location(nc_file)
    global_attributes, params = read_yaml_config(site)
    datef = datetime.strptime(nc_file[-11:-3], "%Y%m%d")
    data_out_l1 = params["data_out"] + "level1/" + datef.strftime("%Y/%m/%d/")
    lev1_file = (
        data_out_l1
        + "MWR_1C01_"
        + global_attributes["wigos_station_id"]
        + "_"
        + datef.strftime("%Y%m%d")
        + ".nc"
    )
    return lev1_file


def _get_freq_flag(data: ndarray, bits: ndarray) -> ndarray:
    """Returns array of flag values for each frequency."""
    flag = np.ones(data.shape) * np.nan
    for i, bit in enumerate(bits):
        flag[isbit(data, bit)] = i + 1
    flag[isbit(data, 0)] = 0
    return flag


def _get_bit_flag(data: ndarray, bits: ndarray) -> ndarray:
    """Returns array of flag values for each bit"""
    flag = np.ones((len(data), len(bits))) * np.nan
    for i, bit in enumerate(bits):
        flag[isbit(data, bit), i] = i
    return flag


def _get_unmasked_values(data: ma.MaskedArray, time: ndarray) -> tuple[ndarray, ndarray]:
    """Returns unmasked time and data"""
    if ma.is_masked(data) is False:
        return data, time
    good_values = ~data.mask
    return data[good_values], time[good_values]


def _nan_time_gaps(time: ndarray, tgap: float = 5.0 / 60.0) -> ndarray:
    """Finds time gaps bigger than 5min (default) and inserts nan."""
    time_diff = np.diff(time)
    gaps = np.where(time_diff > tgap)[0] + 1
    if len(gaps) > 0:
        time[gaps[0 : np.min([len(time), gaps[-1]])]] = np.nan
    return time


def _gap_array(time: ndarray, case_date, tgap: float = 5.0 / 60.0) -> ndarray:
    """Returns edges of time gaps bigger than 5min (default).
    End of gap for current day is current time."""
    locale.setlocale(locale.LC_TIME, "en_US.UTF-8")
    dtnow = datetime.now(tz=timezone.utc)
    day_e = 24.0
    if dtnow.strftime("%d %b %Y") == case_date.strftime("%d %b %Y"):
        day_e = dtnow.hour + dtnow.minute / 60.0 + dtnow.second / 3600.0
        if day_e - time[-1] < 2.0:
            day_e = time[-1]
    time_diff = np.diff(time, prepend=0.0, append=day_e)
    gaps = np.where(time_diff > tgap)[0]
    gtim = np.zeros((len(gaps), 2), np.float32)
    if len(gaps) > 0:
        for i, ind in enumerate(gaps):
            if ind < len(time):
                gtim[i, :] = [time[ind - 1], time[ind]]
    return gtim


def _calculate_rolling_mean(time: ndarray, data: ndarray, win: float = 0.5) -> tuple[ndarray, int]:
    """Returns rolling mean and used window width."""
    width = len(time[time <= time[0] + win])
    if (width % 2) != 0:
        width = width + 1
    if data.ndim == 1:
        rolling_window = np.kaiser(width, 14)
        rolling_mean = np.convolve(data, rolling_window, "valid")
        rolling_mean = rolling_mean / np.sum(rolling_window)
    else:
        rolling_window = np.ones((1, width)) * np.blackman(width)
        rolling_mean = convolve2DFFT(data, rolling_window.T, max_missing=0.1)
    return rolling_mean, width


def _read_location(nc_file: str) -> str:
    """Returns site name."""
    with netCDF4.Dataset(nc_file) as nc:
        site_name = nc.site_location
    return site_name


def heatmap(data, row_labels, col_labels, ax=None, cbar_kw=None, cbarlabel="", **kwargs):
    """
    Create a heatmap from a numpy array and two lists of labels.

    Parameters
    ----------
    data
        A 2D numpy array of shape (M, N).
    row_labels
        A list or array of length M with the labels for the rows.
    col_labels
        A list or array of length N with the labels for the columns.
    ax
        A `matplotlib.axes.Axes` instance to which the heatmap is plotted.  If
        not provided, use current axes or create a new one.  Optional.
    cbar_kw
        A dictionary with arguments to `matplotlib.Figure.colorbar`.  Optional.
    cbarlabel
        The label for the colorbar.  Optional.
    **kwargs
        All other arguments are forwarded to `imshow`.
    """

    if cbar_kw is None:
        cbar_kw = {}

    # Plot the heatmap
    im = ax.imshow(data, **kwargs)

    # Create colorbar
    cbar = ax.figure.colorbar(im, ax=ax, anchor=(1.2, 0.5), **cbar_kw)
    cbar.ax.set_ylabel(cbarlabel, rotation=-90, va="bottom")

    # Show all ticks and label them with the respective list entries.
    ax.set_xticks(np.arange(data.shape[1]), labels=col_labels)
    ax.set_yticks(np.arange(data.shape[0]), labels=row_labels)

    # Turn spines off and create white grid.
    ax.spines[:].set_visible(False)

    ax.set_xticks(np.arange(data.shape[1] + 1) - 0.5, minor=True)
    ax.set_yticks(np.arange(data.shape[0] + 1) - 0.5, minor=True)
    ax.grid(which="minor", color="w", linestyle="-", linewidth=3)
    ax.tick_params(which="minor", bottom=False, left=False)

    return im, cbar


def annotate_heatmap(
    im, data=None, valfmt="{x:.2f}", textcolors=("black", "white"), threshold=None, **textkw
):
    """
    A function to annotate a heatmap.

    Parameters
    ----------
    im
        The AxesImage to be labeled.
    data
        Data used to annotate.  If None, the image's data is used.  Optional.
    valfmt
        The format of the annotations inside the heatmap.  This should either
        use the string format method, e.g. "$ {x:.2f}", or be a
        `matplotlib.ticker.Formatter`.  Optional.
    textcolors
        A pair of colors.  The first is used for values below a threshold,
        the second for those above.  Optional.
    threshold
        Value in data units according to which the colors from textcolors are
        applied.  If None (the default) uses the middle of the colormap as
        separation.  Optional.
    **kwargs
        All other arguments are forwarded to each call to `text` used to create
        the text labels.
    """

    if not isinstance(data, (list, np.ndarray)):
        data = im.get_array()

    im_dat = im.get_array()

    # Normalize the threshold to the images color range.
    if threshold is not None:
        threshold = im.norm(threshold)
    else:
        threshold = im.norm(data.max()) / 2.0

    # Set default alignment to center, but allow it to be
    # overwritten by textkw.
    kw = dict(horizontalalignment="center", verticalalignment="center")
    kw.update(textkw)

    # Get the formatter in case a string is supplied
    if isinstance(valfmt, str):
        valfmt = ticker.StrMethodFormatter(valfmt)

    # Loop over the data and create a `Text` for each "pixel".
    # Change the text's color depending on the data.
    texts = []
    for i in range(data.shape[0]):
        for j in range(data.shape[1]):
            if data[i, j] > 0.0:
                kw.update(color=textcolors[int(im.norm(im_dat[i, j]) < threshold)])
                text = im.axes.text(j, i, valfmt(data[i, j], None), **kw)
                texts.append(text)

    return texts
