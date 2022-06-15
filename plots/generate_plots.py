"""Misc. plotting routines for actris_mwr_pro products."""
from typing import Optional, Tuple
from datetime import datetime, date
import locale
import numpy as np
from numpy import ma
from numpy import ndarray
import netCDF4
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.colors import ListedColormap
from matplotlib.transforms import Affine2D, Bbox, ScaledTranslation
from mpl_toolkits.axes_grid1 import make_axes_locatable
from utils import read_nc_fields, seconds2hours, convolve2DFFT, isbit, read_nc_field_name, get_coeff_list
from plots.plot_meta import ATTRIBUTES, _COLORS
from scipy.signal import filtfilt
from copy import copy
import importlib
from level1.quality_control import spectral_consistency
from read_specs import get_site_specs


class Dimensions:
    """Dimensions of a generated figure in pixels."""

    width: int
    height: int
    margin_top: int
    margin_right: int
    margin_bottom: int
    margin_left: int

    def __init__(self, fig, axes, pad_inches: Optional[float] = None):
        if pad_inches is None:
            pad_inches = rcParams['savefig.pad_inches']

        tightbbox = (
            fig.get_tightbbox(fig.canvas.get_renderer())
            .padded(pad_inches)
            .transformed(Affine2D().scale(fig.dpi))
        )
        self.width = int(tightbbox.width)
        self.height = int(tightbbox.height)

        x0, y0, x1, y1 = (
            Bbox.union([ax.get_window_extent() for ax in axes])
            .translated(-tightbbox.x0, -tightbbox.y0)
            .extents.round()
        )
        self.margin_top = int(self.height - y1)
        self.margin_right = int(self.width - x1 - 1)
        self.margin_bottom = int(y0 - 1)
        self.margin_left = int(x0)


def generate_figure(nc_file: str,
                    field_names: list,
                    show: bool = False,
                    save_path: str = None,
                    max_y: int = 5,
                    ele_range: Tuple = [0., 91.],
                    pointing: int = 0,
                    dpi: int = 120,
                    image_name: Optional[str] = None,
                    sub_title: bool = True,
                    title: bool = True) -> Dimensions:
    """Generates a Cloudnet figure.
    Args:
        nc_file (str): Input file.
        field_names (list): Variable names to be plotted.
        show (bool, optional): If True, shows the figure. Default is True.
        save_path (str, optional): Setting this path will save the figure (in the
            given path). Default is None, when the figure is not saved.
        max_y (int, optional): Upper limit in the plots (km). Default is 12.
        dpi (int, optional): Figure quality (if saved). Higher value means
            more pixels, i.e., better image quality. Default is 120.
        image_name (str, optional): Name (and full path) of the output image.
            Overrides the *save_path* option. Default is None.
        sub_title (bool, optional): Add subtitle to image. Default is True.
        title (bool, optional): Add title to image. Default is True.
    Returns:
        Dimensions of the generated figure in pixels.
    Examples:
        >>> from cloudnetpy.plotting import generate_figure
        >>> generate_figure('categorize_file.nc', ['Z', 'v', 'width', 'ldr', 'beta', 'lwp'])
        >>> generate_figure('iwc_file.nc', ['iwc', 'iwc_error', 'iwc_retrieval_status'])
    """
        
    valid_fields, valid_names = _find_valid_fields(nc_file, field_names)
    is_height = _is_height_dimension(nc_file)
    time = _read_time_vector(nc_file)
    fig, axes = _initialize_figure(len(valid_fields), dpi)
    
    for ax, field, name in zip(axes, valid_fields, valid_names):
        ax.set_facecolor(_COLORS['lightgray'])
        field = _elevation_filter(nc_file, field, ele_range)   
        time = _elevation_filter(nc_file, time, ele_range) 
        if title:
            _set_title(ax, name, nc_file, '')
        if not is_height:
            source = ATTRIBUTES[name].source
            _plot_instrument_data(ax, field, name, source, time, fig, nc_file, ele_range, pointing)
        else:
            ax_value = _read_ax_values(nc_file)            
            ax_value = (time, ax_value[1])
            field, ax_value = _screen_high_altitudes(field, ax_value, max_y)
            _set_ax(ax, max_y)
            
            plot_type = ATTRIBUTES[name].plot_type
            if plot_type == 'segment':
                _plot_segment_data(ax, field, name, ax_value)

            elif plot_type == 'mesh':
                _plot_colormesh_data(ax, field, name, ax_value, max_y, nc_file)
            
    case_date = _set_labels(fig, axes[-1], nc_file, sub_title)
    _handle_saving(image_name, save_path, show, case_date, valid_names)
    return Dimensions(fig, axes)


def _mark_gaps(time: np.ndarray,
               data: ma.MaskedArray,
               max_allowed_gap: float = 1, 
               mask_edge: int = 0) -> tuple:

    assert time[0] >= 0
    assert time[-1] <= 24
    max_gap = max_allowed_gap / 60
    if not ma.is_masked(data):
        mask_new = np.zeros(data.shape)
    elif ma.all(data.mask) is ma.masked:
        mask_new = np.ones(data.shape)
    else:
        mask_new = np.copy(data.mask)
    data_new = ma.copy(data)
    time_new = np.copy(time)
    gap_indices = np.where(np.diff(time) > max_gap)[0]
    for ia in range(mask_edge):
        gap_indices = np.unique(np.sort(np.append(gap_indices, [gap_indices-ia, gap_indices+ia])))
    temp_array = np.zeros((2, data.shape[1]))
    temp_mask = np.ones((2, data.shape[1]))
    time_delta = 0.
    for ind in np.sort(gap_indices)[::-1]:
        ind += 1
        data_new = np.insert(data_new, ind, temp_array, axis=0)
        mask_new = np.insert(mask_new, ind, temp_mask, axis=0)
        time_new = np.insert(time_new, ind, time[ind] - time_delta)
        time_new = np.insert(time_new, ind, time[ind-1] + time_delta)
    if (time[0] - 0) > max_gap:
        data_new = np.insert(data_new, 0, temp_array, axis=0)
        mask_new = np.insert(mask_new, 0, temp_mask, axis=0)
        time_new = np.insert(time_new, 0, time[0] - time_delta)
        time_new = np.insert(time_new, 0, time_delta)
    if (24 - time[-1]) > max_gap:
        ind = mask_new.shape[0]
        data_new = np.insert(data_new, ind, temp_array, axis=0)
        mask_new = np.insert(mask_new, ind, temp_mask, axis=0)
        time_new = np.insert(time_new, ind, 24 - time_delta)
        time_new = np.insert(time_new, ind, time[-1] + time_delta)
    data_new.mask = mask_new
    return time_new, data_new


def _handle_saving(image_name: Optional[str], save_path: Optional[str], show: bool,
                   case_date: date, field_names: list, fix: str = ""):
    if image_name:
        date_string = case_date.strftime("%Y%m%d")
        plt.savefig(f"{save_path}{date_string}_{image_name}.png", bbox_inches='tight')
    elif save_path:
        file_name = _create_save_name(save_path, case_date, field_names, fix)
        plt.savefig(file_name, bbox_inches='tight')
    if show:
        plt.show()
    plt.close()


def _set_labels(fig, ax, nc_file: str, sub_title: bool = True) -> date:
    ax.set_xlabel('Time (UTC)', fontsize=13)
    case_date = _read_date(nc_file)
    site_name = _read_location(nc_file)
    if sub_title:
        _add_subtitle(fig, case_date, site_name)
    return case_date


def _set_title(ax, field_name: str, nc_file, identifier: str = " from actris_mwr_pro"):
    if ATTRIBUTES[field_name].name:
        ax.set_title(f"{ATTRIBUTES[field_name].name}{identifier}", fontsize=14)
    else:
        ax.set_title(f"{read_nc_field_name(nc_file, field_name)}{identifier}", fontsize=14)


def _find_valid_fields(nc_file: str, names: list) -> Tuple[list, list]:
    """Returns valid field names and corresponding data."""
    valid_names, valid_data = names[:], []
    with netCDF4.Dataset(nc_file) as nc:
        for name in names:
            if name in nc.variables:
                valid_data.append(nc.variables[name][:])
            else:
                valid_names.remove(name)
    if not valid_names:
        raise ValueError('No fields to be plotted')
    return valid_data, valid_names


def _is_height_dimension(full_path: str) -> bool:
    with netCDF4.Dataset(full_path) as nc:
        is_height = 'altitude' in nc.variables
    return is_height


def _elevation_filter(full_path: str, data_field: ndarray, ele_range: Tuple) -> ndarray:
    with netCDF4.Dataset(full_path) as nc:
        if 'ele' in nc.variables:
            elevation = read_nc_fields(full_path, 'ele')
            if data_field.ndim > 1:
                data_field = data_field[(elevation >= ele_range[0]) & (elevation <= ele_range[1]), :]
            else:
                data_field = data_field[(elevation >= ele_range[0]) & (elevation <= ele_range[1])]
    return data_field


def _pointing_filter(full_path: str, data_field: ndarray, ele_range: Tuple, status: int) -> ndarray:
    with netCDF4.Dataset(full_path) as nc:
        if 'pointing_flag' in nc.variables:
            pointing = read_nc_fields(full_path, 'pointing_flag')
            pointing = _elevation_filter(full_path, pointing, ele_range)
            if data_field.ndim > 1:
                data_field = data_field[pointing == status, :]
            else:
                data_field = data_field[pointing == status]
    return data_field


def _initialize_figure(n_subplots: int, dpi) -> tuple:
    """Creates an empty figure according to the number of subplots."""
    fig, axes = plt.subplots(n_subplots, figsize=(16, 4 + (n_subplots-1)*4.8), dpi=dpi, facecolor='white')
    fig.subplots_adjust(left=0.06, right=0.73)
    if n_subplots == 1:
        axes = [axes]
    return fig, axes


def _read_ax_values(full_path: str) -> Tuple[ndarray, ndarray]:
    """Returns time and height arrays."""
    fields = ['time', 'altitude']
    time, height = read_nc_fields(full_path, fields)
    height_km = height / 1000
    return time, height_km


def _read_time_vector(nc_file: str) -> ndarray:
    """Converts time vector to fraction hour."""
    with netCDF4.Dataset(nc_file) as nc:
        time = nc.variables['time'][:]
    return seconds2hours(time)


def _screen_high_altitudes(data_field: ndarray, ax_values: tuple, max_y: int) -> tuple:
    """Removes altitudes from 2D data that are not visible in the figure.
    Bug in pcolorfast causing effect to axis not noticing limitation while
    saving fig. This fixes that bug till pcolorfast does fixing themselves.
    Args:
        data_field (ndarray): 2D data array.
        ax_values (tuple): Time and height 1D arrays.
        max_y (int): Upper limit in the plots (km).
    """
    alt = ax_values[-1]
    if data_field.ndim > 1:
        ind = int((np.argmax(alt > max_y) or len(alt)) + 1)
        data_field = data_field[:, :ind]
        alt = alt[:ind]
    return data_field, (ax_values[0], alt)


def _set_ax(ax, max_y: float, ylabel: str = None, min_y: float = 0.0):
    """Sets ticks and tick labels for plt.imshow()."""
    ticks_x_labels = _get_standard_time_ticks()
    ax.set_ylim(min_y, max_y)
    ax.set_xticks(np.arange(0, 25, 4, dtype=int))
    ax.set_xticklabels(ticks_x_labels, fontsize=12)
    ax.set_ylabel('Height (km)', fontsize=13)
    ax.set_xlim(0, 24)
    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=13)


def _get_standard_time_ticks(resolution: int = 4) -> list:
    """Returns typical ticks / labels for a time vector between 0-24h."""
    return [f"{int(i):02d}:00" if 24 > i > 0 else ''
            for i in np.arange(0, 24.01, resolution)]


def _plot_segment_data(ax, data: ma.MaskedArray, name: str, axes: tuple, nc_file: str):
    """Plots categorical 2D variable.
    Args:
        ax (obj): Axes object of subplot (1,2,3,.. [1,1,],[1,2]... etc.)
        data (ndarray): 2D data array.
        name (string): Name of plotted data.
        axes (tuple): Time and height 1D arrays.
    """
    variables = ATTRIBUTES[name]
    clabel = [x[0] for x in variables.clabel]
    cbar = [x[1] for x in variables.clabel]
    cmap = ListedColormap(cbar)
    pl = ax.pcolor(*axes, data.T, cmap=cmap, shading='nearest', vmin=-0.5, vmax=len(cbar) - 0.5)
    ax.grid(axis='y')
    colorbar = _init_colorbar(pl, ax)
    colorbar.set_ticks(np.arange(len(clabel)))  
    if name == 'quality_flag_2':
        site_name = _read_location(nc_file)
        site_dict = importlib.import_module(f"site_config." + site_name + ".config")
        clabel[2] = clabel[2] + ' (' + str(site_dict.params['TB_threshold'][1]) + ' K)'
        clabel[1] = clabel[1] + ' (' + str(site_dict.params['TB_threshold'][0]) + ' K)'
    colorbar.ax.set_yticklabels(clabel, fontsize=13)



def _plot_colormesh_data(ax, data: ma.MaskedArray, name: str, axes: tuple, max_y: int, nc_file: str):
    """Plots continuous 2D variable.
    Creates only one plot, so can be used both one plot and subplot type of figs.
    Args:
        ax (obj): Axes object of subplot (1,2,3,.. [1,1,],[1,2]... etc.)
        data (ndarray): 2D data array.
        name (string): Name of plotted data.
        axes (tuple): Time and height 1D arrays.
    """
    variables = ATTRIBUTES[name]
    if ATTRIBUTES[name].nlev:
        nlev = ATTRIBUTES[name].nlev
    else:
        nlev = 16
        
    assert variables.plot_range is not None
     
    if name == 'relative_humidity':
        data[data > 1.] = 1.
        data*=100.

    if variables.plot_type == 'bit':
        cmap = ListedColormap(variables.cbar)
        pos = ax.get_position()
        ax.set_position([pos.x0, pos.y0, pos.width * 0.965, pos.height])
    else:
        cmap = plt.get_cmap(variables.cbar, nlev)

    vmin, vmax = variables.plot_range
    if np.ma.median(np.diff(axes[0][:])) < 5/60:
        data, width = _calculate_rolling_mean(axes[0], data, win=30/60)    
        time, data = _mark_gaps(axes[0][:], data, 35, 20)
    else:
        time, data = _mark_gaps(axes[0][:], data, 35, 0)
    pl = ax.contourf(time, axes[1], data.T, levels=np.linspace(vmin,vmax,nlev), cmap=cmap, extend=variables.cbar_ext)
    ds = int(np.round(len(time)*.15))
    cp = ax.contour(time[ds:len(time)-ds], axes[1], data[ds:len(time)-ds, :].T, levels=np.linspace(vmin,vmax,nlev), colors='black', linewidths=.0001)
    cbl = plt.clabel(cp, fontsize=8)
    ta = []
    for l in cbl:
        l.set_va('bottom')
        l.set_weight('bold')
        if float(l.get_text()) in ta:
            l.set_visible(False)
        ta = np.append(ta, [float(l.get_text())])
    cp = ax.contour(time, axes[1], data.T, levels=np.linspace(vmin,vmax,nlev), colors='black', linewidths=.8)
   
    if variables.plot_type != 'bit':
        colorbar = _init_colorbar(pl, ax)
        locator = colorbar.ax.yaxis.get_major_locator()
        # locator, formatter = colorbar._get_ticker_locator_formatter()
        locator.set_params(nbins=6)
        colorbar.update_ticks()
        colorbar.set_label(variables.clabel, fontsize=13)


def _plot_instrument_data(ax, data: ma.MaskedArray, name: str, product: Optional[str], time: ndarray, fig, nc_file: str, ele_range: Tuple, pointing: int):
    if product == 'int':
        _plot_int(ax, data, name, time)
    elif product == 'met':
        _plot_met(ax, data, name, time, nc_file)
    elif product == 'tb':
        _plot_tb(ax, data, time, fig, nc_file, ele_range, pointing, name)
    elif product == 'irt':
        _plot_irt(ax, data, name, time, nc_file)        
    elif product == 'qf':
        _plot_qf(ax, data, time, fig, nc_file)
    elif product == 'mqf':
        _plot_mqf(ax, data, time)       
    elif product == 'sen':
        _plot_sen(ax, data, name, time, nc_file)
        
    pos = ax.get_position()
    ax.set_position([pos.x0, pos.y0, pos.width * 0.965, pos.height])
    
    
def _plot_sen(ax, data_in: ndarray, name: str, time: ndarray, nc_file: str):
    variables = ATTRIBUTES[name]
    pointing_flag = read_nc_fields(nc_file, 'pointing_flag')
    quality_flag = read_nc_fields(nc_file, 'quality_flag')
    qf = _get_freq_data(quality_flag, [6])
    vmin, vmax = variables.plot_range    
    
    ax.plot(time[pointing_flag == 0], data_in[pointing_flag == 0], '--.', color=_COLORS['darkgreen'], label='single pointing', linewidth=.8)
    if name == 'ele':
        ax.plot(time[(pointing_flag == 1) & (data_in > 0.)], data_in[(pointing_flag == 1) & (data_in > 0.)], '--.', alpha=.75, color=_COLORS['green'], label='multiple pointing', linewidth=.8)
    else:
        ax.plot(time[pointing_flag == 1], data_in[pointing_flag == 1], '--.', alpha=.75, color=_COLORS['green'], label='multiple pointing', linewidth=.8)
    ax.plot(time[np.any(qf == 1, axis=1)], data_in[np.any(qf == 1, axis=1)], 'r.', linewidth=1, label='sun_in_beam')
    ax.set_ylim((vmin, vmax))
    ax.legend(loc='upper right')
    _set_ax(ax, vmax, variables.ylabel, vmin)

    
def _plot_irt(ax, data_in: ndarray, name: str, time: ndarray, nc_file: str):
    variables = ATTRIBUTES[name]
    vmin, vmax = variables.plot_range
    ir_wavelength = read_nc_fields(nc_file, 'ir_wavelength')
    if not ma.all(ma.is_masked(data_in[:, 0])):
        ax.plot(time, data_in[:,0],'o', markersize=.75, fillstyle='full', color='sienna', label=str(ir_wavelength[0])+' µm')
    if not ma.all(ma.is_masked(data_in[:, 1])):    
        ax.plot(time, data_in[:,1],'o', markersize=.75, fillstyle='full', color=_COLORS['shockred'], label=str(ir_wavelength[1])+' µm')
    ax.set_ylim((vmin, vmax))
    ax.legend(loc='upper right')
    _set_ax(ax, vmax, variables.ylabel, vmin)
    
    
def _plot_mqf(ax, data_in: ndarray, time: ndarray):
    qf = _get_bit_data(data_in, np.arange(6))
    _plot_segment_data(ax, qf, 'met_quality_flag', (time, np.linspace(.5,5.5,6)))
    ax.set_yticks(np.arange(6))
    ax.yaxis.set_ticklabels([]) 
    _set_ax(ax, 6, '')  
    ax.set_title(ATTRIBUTES['met_quality_flag'].name)
    
    
def _plot_qf(ax, data_in: ndarray, time: ndarray, fig, nc_file: str):
    fig.clear()
    fig, axs = plt.subplots(3,1, figsize=(12.52, 16), dpi=120, facecolor='w')        
    frequency = read_nc_fields(nc_file, 'frequency')
 
    qf = _get_bit_data(data_in[:,0], [4,5,6])
    _plot_segment_data(axs[0], qf, 'quality_flag_0', (time, np.linspace(.5,2.5,3)), nc_file)
    axs[0].set_yticks(np.arange(3))
    axs[0].yaxis.set_ticklabels([]) 
    _set_ax(axs[0], 3, '')  
    axs[0].set_title(ATTRIBUTES['quality_flag_0'].name)
    
    qf = _get_freq_data(data_in, [3])
    _plot_segment_data(axs[1], qf, 'quality_flag_1', (time, np.linspace(.5,len(frequency)-.5,len(frequency))), nc_file)    
    axs[1].set_title(ATTRIBUTES['quality_flag_1'].name)
    
    qf = _get_freq_data(data_in, [1,2])
    _plot_segment_data(axs[2], qf, 'quality_flag_2', (time, np.linspace(.5,len(frequency)-.5,len(frequency))), nc_file)
    axs[2].set_title(ATTRIBUTES['quality_flag_2'].name)
    
    offset = ScaledTranslation(0/72, 8/72, fig.dpi_scale_trans)    
    for i in [1,2]:        
        axs[i].set_yticks(range(len(frequency)))
        axs[i].set_yticklabels(frequency)
        for label in axs[i].yaxis.get_majorticklabels():
            label.set_transform(label.get_transform() + offset)        
        _set_ax(axs[i], len(frequency), 'Frequency [GHz]')    
        axs[i].plot(np.linspace(*axs[i].get_xlim(), len(time)), np.ones(len(time))*7, 'k-', linewidth=1)
    _set_labels(fig, axs[-1], nc_file)        
    
    
def _plot_tb(ax, data_in: ndarray, time: ndarray, fig, nc_file: str, ele_range: Tuple, pointing: int, name: str):   
    frequency = read_nc_fields(nc_file, 'frequency')
    if name == 'tb': 
        quality_flag = read_nc_fields(nc_file, 'quality_flag')
        data_in = _pointing_filter(nc_file, data_in, ele_range, pointing)
        time = _pointing_filter(nc_file, time, ele_range, pointing)
        quality_flag = _elevation_filter(nc_file, quality_flag, ele_range)
        quality_flag = _pointing_filter(nc_file, quality_flag, ele_range, pointing)
    else:
        sc = dict()
        sc['tb'] = read_nc_fields(nc_file, 'tb')
        sc['receiver_nb'] = read_nc_fields(nc_file, 'receiver_nb')
        sc['receiver'] = read_nc_fields(nc_file, 'receiver')
        sc['frequency'] = read_nc_fields(nc_file, 'frequency')
        sc['ele'] = read_nc_fields(nc_file, 'ele')
        sc['time'] = read_nc_fields(nc_file, 'time')
        sc['pointing_flag'] = np.zeros(len(sc['time']))
        site = _read_location(nc_file)
        c_list = get_coeff_list(site, 'tbx')
        _, params = get_site_specs(site, '1C01')
        quality_flag, _ = spectral_consistency(sc, c_list, params['tbx_f'])        
        data_in = np.abs(data_in - sc['tb'])

        
    fig.clear()
    fig, axs = plt.subplots(7,2, figsize=(13, 16), facecolor='w', edgecolor='k', sharex=True, dpi=120)
    fig.subplots_adjust(hspace=.035, wspace=.15)
    fig.text(.06, .5, 'Brightness Temperature [K]', va='center', rotation='vertical', fontsize=20)
    fig.text(.445, .09, 'flagged data', va='center', fontsize=20, color='r')
    axs[0,0].set_title('K-Band Channels', fontsize=15, color=_COLORS['darkgray'], loc='right')
    axs[0,1].set_title('V-Band Channels', fontsize=15, color=_COLORS['darkgray'], loc='right')
    trans = ScaledTranslation(10/72, -5/72, fig.dpi_scale_trans)
    tb_m, tb_s = [], []
    
    for i, ax in enumerate(axs.T.flatten()):
        no_flag = np.where(quality_flag[:,i] == 0)[0]
        if len(np.array(no_flag)) == 0:
            no_flag = np.arange(len(time))        
        tb_m = np.append(tb_m, [np.mean(data_in[no_flag, i])]) #TB mean
        tb_s = np.append(tb_s, [np.std(data_in[no_flag, i])]) #TB std        
        ax.plot(time, np.ones(len(time))*tb_m[i], '--', color=_COLORS['darkgray'], linewidth=1)
        ax.plot(time, data_in[:,i],'ko', markersize=.75, fillstyle='full')
        flag = np.where(quality_flag[:,i] > 0)[0]
        ax.plot(time[flag], data_in[flag,i],'ro', markersize=.75, fillstyle='full')
        ax.set_facecolor(_COLORS['lightgray'])
        dif = np.max(data_in[no_flag,i]) - np.min(data_in[no_flag,i])
        _set_ax(ax, np.max(data_in[no_flag,i])+.25*dif, '', np.max([0., np.min(data_in[no_flag,i])-.25*dif]))
        _set_labels(fig, ax, nc_file)
        ax.text(.05,.9,str(frequency[i])+' GHz', transform=ax.transAxes + trans, color=_COLORS['darkgray'], fontweight='bold')
        ax.text(.55,.9,str(round(tb_m[i],2))+' +/- '+str(round(tb_s[i],2))+' K', transform=ax.transAxes+trans, color=_COLORS['darkgray'], fontweight='bold')
        
    if name == 'tb_spectrum':
        tb_m = np.ones((len(sc['time']), len(sc['receiver_nb']))) * np.nan
        for irec, rec in enumerate(sc['receiver_nb']):
            tb_m[:, irec] = np.sum(data_in[:, sc['receiver'] == rec], axis=1)/np.median(np.sum(data_in[:, sc['receiver'] == rec], axis=1))  

    "TB mean"
    axa = fig.add_subplot(121)
    axa.set_position([.125, -.05, .72, .125])
    axa.tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False)
    axa.tick_params(axis='y', which='both', right=False, left=False, labelleft=False)
    for pos in ['right','top','bottom','left']:
        axa.spines[pos].set_visible(False)
    axa.set_facecolor(_COLORS['lightgray'])
    
    if name == 'tb':
        axaK = fig.add_subplot(121)
        axaK.set_position([.125, -.05, .36, .125])
        axaK.plot(frequency, tb_m, 'ko', markerfacecolor='k', markersize=4)
        axaK.errorbar(frequency, tb_m, yerr=tb_s, xerr=None, linestyle='', capsize=8, color='k')
        axaK.set_xticks(frequency)
        axaK.set_xticklabels(axaK.get_xticks(), rotation = 30)
        axaK.set_xlim([22, 32])
        axaK.set_ylim([0, np.nanmax(tb_m+tb_s)+30])
        axaK.tick_params(axis='both', labelsize=12)
        axaK.set_facecolor(_COLORS['lightgray'])

        axaV = fig.add_subplot(122)
        axaV.set_position([.54, -.05, .36, .125])
        axaV.plot(frequency, tb_m, 'ko', markerfacecolor='k', markersize=4)
        axaV.errorbar(frequency, tb_m, yerr=tb_s, xerr=None, linestyle='', capsize=8, color='k')
        axaV.set_xticks(frequency)
        axaV.set_xticklabels(axaV.get_xticks(), rotation = 30)    
        axaV.set_xlim([51, 58.5])    
        axaV.set_ylim([0, np.nanmax(tb_m+tb_s)+30])
        axaV.tick_params(axis='both', labelsize=12)
        axaV.set_facecolor(_COLORS['lightgray'])

        axaK.spines['right'].set_visible(False)
        axaV.spines['left'].set_visible(False)
        axaV.yaxis.tick_right()    
        d = .015
        axaK.plot((1-d,1+d), (-d,+d), transform=axaK.transAxes, color='k', clip_on=False)
        axaK.plot((1-d,1+d),(1-d,1+d), transform=axaK.transAxes, color='k', clip_on=False)
        axaV.plot((-d,+d), (1-d,1+d), transform=axaV.transAxes, color='k', clip_on=False)
        axaV.plot((-d,+d), (-d,+d), transform=axaV.transAxes, color='k', clip_on=False)    
        axaK.set_ylabel('Brightness Temperature [K]', fontsize=12)
        axaV.text(-.08, .9, 'TB daily means +/- standard deviation', fontsize=13, horizontalalignment='center', transform=axaV.transAxes, color=_COLORS['darkgray'], fontweight='bold')
        axaV.text(-.08, -.35, 'Frequency [GHz]', fontsize=13, horizontalalignment='center', transform=axaV.transAxes)
        
    else:
        axa = fig.add_subplot(111)
        axa.set_position([.125, -.05, .775, .125])
        cl = [_COLORS['darksky'], _COLORS['lightpurple']]
        lb = ['K-Band', 'V-Band']
        for irec, _ in enumerate(sc['receiver_nb']):
            axa.plot(time, tb_m[:, irec], color=cl[irec], markersize=.75, fillstyle='full', label=lb[irec])
        axa.legend(loc='upper right')
        axa.set_ylabel('Brightness Temperature [K]', fontsize=12)
        axa.set_xlabel('Time (UTC)', fontsize=12)
        ticks_x_labels = _get_standard_time_ticks()
        axa.set_xticks(np.arange(0, 25, 4, dtype=int))
        axa.set_xticklabels(ticks_x_labels, fontsize=12)
        axa.set_xlim(0, 24)        
        axa.set_ylim([0, np.nanmax(tb_m)+.5])
        axa.text(.5, .9, 'Relative deviations retrieved from observed TB', fontsize=13, horizontalalignment='center', transform=axa.transAxes, color=_COLORS['darkgray'], fontweight='bold')
    
    
def _plot_met(ax, data_in: ndarray, name: str, time: ndarray, nc_file: str):
    ylabel = ATTRIBUTES[name].ylabel
    if name == 'relative_humidity':
        data_in*=100.
        ylabel = '%'
    data, time = _get_unmasked_values(data_in, time)
    if name == 'wind_direction':
        spd = read_nc_fields(nc_file, 'wind_speed')
        rolling_mean, width = _dir_avg(time, spd, data)
    else:        
        rolling_mean, width = _calculate_rolling_mean(time, data)
    time = _nan_time_gaps(time)
    rolling_mean = np.interp(time, time[int(width / 2 - 1):int(-width / 2)], rolling_mean)
    
    ax.plot(time, data, '.', alpha=.8, color=_COLORS['darksky'], markersize=1)
    ax.plot(time, rolling_mean, 'o', fillstyle='full', color='darkblue', markersize=3)
    vmin, vmax = ATTRIBUTES[name].plot_range    
    if (name == 'air_temperature') | (name == 'air_pressure'):
        vmin, vmax = np.nanmin(data)-1., np.nanmax(data)+1.  
    if (name == 'rain_rate') | (name == 'wind_speed'):
        vmin, vmax = 0., np.nanmax([np.nanmax(data)+1., 3.])         
    _set_ax(ax, vmax, ylabel, min_y=vmin)    
    ax.grid(True)
    

def _plot_int(ax, data_in: ma.MaskedArray, name: str, time: ndarray):
    data, time = _get_unmasked_values(data_in, time)
    # data = _filter_noise(data)
    rolling_mean, width = _calculate_rolling_mean(time, data)
    time = _nan_time_gaps(time)
    rolling_mean = np.interp(time, time[int(width / 2 - 1):int(-width / 2)], rolling_mean)
    
    ax.plot(time, data, color='royalblue', lw=.5)
    ax.axhline(linewidth=0.8, color='k')
    ax.plot(time, rolling_mean, color='sienna', linewidth=2.0)
    ax.plot(time, rolling_mean, color='wheat', linewidth=0.6)
    vmin, vmax = ATTRIBUTES[name].plot_range
    if name == 'iwv':
        vmin, vmax = np.nanmin(data)-1., np.nanmax(data)+1.   
    _set_ax(ax, vmax, ATTRIBUTES[name].ylabel, min_y=vmin)
    
    
def _dir_avg(time: ndarray, 
             spd: ndarray, 
             drc: ndarray,
             win: float=.5):
    width = len(time[time <= time[0] + win])
    if (width % 2) != 0:
        width = width + 1
    seq = range(len(time))
    dir_avg = []
    for i in range(len(seq) - width + 1):
        dir_avg.append(_windvec(spd[seq[i: i + width]], drc[seq[i: i + width]]))
    return np.array(dir_avg), width


def _windvec(spd: ndarray,
             drc: ndarray):    
    ve = -np.mean(spd * np.sin(np.deg2rad(drc)))
    vn = -np.mean(spd * np.cos(np.deg2rad(drc)))
    uv = np.sqrt(ve * ve + vn * vn)
    vdir = np.rad2deg(np.arctan2(ve, vn))
    if vdir < 180.0:
        Dv = vdir + 180.0
    elif vdir > 180.0:
        Dv = vdir - 180
    return Dv    

    
def _filter_noise(data: ndarray) -> ndarray:
    """IIR filter"""
    n = int(np.rint(np.nextafter((len(data) / 1000), (len(data) / 1000)+1)))
    if n <= 1:
        n = 2
    b = [1.0 / n] * n
    a = 1
    return filtfilt(b, a, data)


def _get_freq_data(data: ndarray, bits: ndarray) -> ndarray:
    flag = np.ones(data.shape)*np.nan
    for i, bit in enumerate(bits):
        flag[isbit(data, bit)] = i+1
    flag[isbit(data, 0)] = 0
    return flag


def _get_bit_data(data: ndarray, bits: ndarray) -> ndarray:
    flag = np.ones((len(data), len(bits)))*np.nan
    for i, bit in enumerate(bits):
        flag[isbit(data, bit), i] = i
    return flag


def _get_unmasked_values(data: ma.MaskedArray, time: ndarray) -> Tuple[np.ndarray, np.ndarray]:
    if ma.is_masked(data) is False:
        return data, time
    good_values = ~data.mask
    return data[good_values], time[good_values]


def _nan_time_gaps(time: ndarray) -> ndarray:
    """Finds time gaps bigger than 5min and inserts nan."""
    time_diff = np.diff(time)
    dec_hour_5min = 0.085
    gaps = np.where(time_diff > dec_hour_5min)[0]
    time[gaps] = np.nan
    return time


def _calculate_rolling_mean(time: ndarray, data: ndarray, win: float = .5) -> Tuple[ndarray, int]:
    width = len(time[time <= time[0] + win])
    if (width % 2) != 0:
        width = width + 1
    if data.ndim == 1:
        rolling_window = np.blackman(width)
        rolling_mean = np.convolve(data, rolling_window, 'valid')
        rolling_mean = rolling_mean / np.sum(rolling_window)
    else:
        rolling_window = np.ones((1, width))*np.blackman(width)
        rolling_mean = convolve2DFFT(data, rolling_window.T, max_missing=.1)
    return rolling_mean, width


def _init_colorbar(plot, axis):
    divider = make_axes_locatable(axis)
    cax = divider.append_axes("right", size="1%", pad=0.25)
    return plt.colorbar(plot, fraction=1.0, ax=axis, cax=cax)


def _read_location(nc_file: str) -> str:
    """Returns site name."""
    with netCDF4.Dataset(nc_file) as nc:
        site_name = nc.site
    return site_name


def _read_date(nc_file: str) -> date:
    """Returns measurement date."""
    locale.setlocale(locale.LC_TIME, 'en_US.UTF-8')
    with netCDF4.Dataset(nc_file) as nc:
        case_date = datetime.strptime(nc.date, "%Y-%m-%d")
    return case_date


def _add_subtitle(fig, case_date: date, site_name: str):
    """Adds subtitle into figure."""
    text = _get_subtitle_text(case_date, site_name)
    fig.suptitle(text, fontsize=13, y=0.885, x=0.07, horizontalalignment='left',
                 verticalalignment='bottom')


def _get_subtitle_text(case_date: date, site_name: str) -> str:
    site_name = site_name.replace('-', ' ')
    return f"{site_name}, {case_date.strftime('%-d %b %Y')}"


def _create_save_name(save_path: str, case_date: date, field_names: list, fix: str = '') -> str:
    """Creates file name for saved images."""
    date_string = case_date.strftime("%Y%m%d")
    return f"{save_path}{date_string}_{'_'.join(field_names)}{fix}.png"
