"""RpgArray Class"""
import locale
from datetime import datetime, timezone

import netCDF4
import numpy as np
from numpy import ma

import utils
import version
from level1.lev1_meta_nc import MetaData


class RpgArray:
    """Stores netCDF4 variables, numpy arrays and scalars as RpgArrays.
    Args:
        variable: The netCDF4 :class:`Variable` instance, numpy array (masked or regular),
            or scalar (float, int).
        name: Name of the variable.
        units_from_user: Units of the variable.
    Attributes:
        name (str): Name of the variable.
        data (ndarray): The actual data.
        data_type (str): 'i4' for integers, 'f4' for floats.
        units (str): The `units_from_user` argument if it is given. Otherwise
            copied from the original netcdf4 variable. Empty if input is just data.
    """

    def __init__(
        self,
        variable: netCDF4.Variable | np.ndarray | float | int,
        name: str,
        units_from_user: str | None = None,
        dimensions: str | None = None,
    ):
        self.variable = variable
        self.name = name
        self.data = self._init_data()
        self.units = self._init_units(units_from_user)
        self.data_type = self._init_data_type()
        self.dimensions = dimensions

    def fetch_attributes(self) -> list:
        """Returns list of user-defined attributes."""

        attributes = []
        for attr in self.__dict__:
            if attr not in ("name", "data", "data_type", "variable", "dimensions"):
                attributes.append(attr)
        return attributes

    def set_attributes(self, attributes: MetaData) -> None:
        """Overwrites existing instance attributes."""

        for key in attributes._fields:  # To iterate namedtuple fields.
            data = getattr(attributes, key)
            if data:
                setattr(self, key, data)

    def _init_data(self) -> np.ndarray:
        if isinstance(self.variable, netCDF4.Variable):
            return self.variable[:]
        if isinstance(self.variable, np.ndarray):
            return self.variable
        if isinstance(self.variable, (int, float)):
            return np.array(self.variable)
        if isinstance(self.variable, str):
            try:
                numeric_value = utils.str_to_numeric(self.variable)
                return np.array(numeric_value)
            except ValueError:
                pass
        raise ValueError(f"Incorrect RpgArray input: {self.variable}")

    def _init_units(self, units_from_user: str | None) -> str:
        if units_from_user is not None:
            return units_from_user
        return getattr(self.variable, "units", "")

    def _init_data_type(self) -> str:
        if self.data.dtype in (np.float32, np.float64):
            return "f4"
        return "i4"

    def __getitem__(self, ind: tuple) -> np.ndarray:
        return self.data[ind]


class Rpg:
    """Base class for RPG MWR."""

    def __init__(self, raw_data: dict):
        self.raw_data = raw_data
        self.date = self._get_date()
        self.data = {}
        self.data = self._init_data()

    def _init_data(self) -> dict:
        data = {}
        for key in self.raw_data:
            data[key] = RpgArray(self.raw_data[key], key)
        return data

    def _get_date(self):
        epoch = datetime(1970, 1, 1).timestamp()
        time_median = float(ma.median(self.raw_data["time"]))
        time_median += epoch
        return datetime.utcfromtimestamp(time_median).strftime("%Y-%m-%d")

    def find_valid_times(self):
        # sort timestamps
        time = self.data["time"].data[:]
        ind = time.argsort()
        self._screen(ind)

        # remove duplicate timestamps
        time = self.data["time"].data[:]
        _, ind = np.unique(time, return_index=True)
        self._screen(ind)

        # find valid date
        time = self.data["time"].data[:]
        ind = np.zeros(len(time), dtype=int)
        for i, t in enumerate(time):
            if "-".join(utils.seconds2date(t)[:3]) == self.date:
                ind[i] = 1
        self._screen(np.where(ind == 1)[0])

    def _screen(self, ind: np.ndarray):
        if len(ind) < 1:
            raise RuntimeError(["Error: no valid data for date: " + self.date])
        n_time = len(self.data["time"].data)
        for _, array in self.data.items():
            data = array.data
            if data.ndim > 0 and data.shape[0] == n_time:
                if data.ndim == 1:
                    screened_data = data[ind]
                else:
                    screened_data = data[ind, :]
                data = screened_data


def save_rpg(rpg: Rpg, output_file: str, att: dict, data_type: str) -> None:
    """Saves the RPG MWR file."""

    if data_type == "1B01":
        dims = {
            "time": len(rpg.data["time"][:]),
            "frequency": len(rpg.data["tb"][:].T),
            "receiver_nb": len(rpg.data["receiver_nb"][:]),
            "bnds": 2,
        }
    elif data_type == "1B11":
        dims = {
            "time": len(rpg.data["time"][:]),
            "ir_wavelength": len(rpg.data["irt"][:].T),
        }
    elif data_type == "1B21":
        dims = {"time": len(rpg.data["time"][:])}
    elif data_type == "1C01":
        if "irt" in rpg.data:
            dims = {
                "time": len(rpg.data["time"][:]),
                "frequency": len(rpg.data["tb"][:].T),
                "receiver_nb": len(rpg.data["receiver_nb"][:]),
                "ir_wavelength": len(rpg.data["irt"][:].T),
                "bnds": 2,
            }
        else:
            dims = {
                "time": len(rpg.data["time"][:]),
                "frequency": len(rpg.data["tb"][:].T),
                "receiver_nb": len(rpg.data["receiver_nb"][:]),
                "bnds": 2,
            }
    elif data_type in ("2P01", "2P02", "2P03", "2P04", "2P07", "2P08"):
        dims = {
            "time": len(rpg.data["time"][:]),
            "bnds": 2,
            "altitude": len(rpg.data["altitude"][:]),
        }
    elif data_type in ("2I01", "2I02"):
        dims = {"time": len(rpg.data["time"][:]), "bnds": 2}
    elif data_type == "2S02":
        dims = {
            "time": len(rpg.data["time"][:]),
            "bnds": 2,
            "receiver_nb": len(rpg.data["receiver_nb"][:]),
            "frequency": len(rpg.data["tb_spectrum"][:].T),
        }
    else:
        raise RuntimeError(["Data type " + data_type + " not supported for file writing."])

    with init_file(output_file, dims, rpg.data, att) as rootgrp:
        setattr(rootgrp, "date", rpg.date)
        # setattr(rootgrp, 'site', att['site_location'])


def init_file(
    file_name: str, dimensions: dict, rpg_arrays: dict, att_global: dict
) -> netCDF4.Dataset:
    """Initializes a RPG MWR file for writing.
    Args:
        file_name: File name to be generated.
        dimensions: Dictionary containing dimension for this file.
        rpg_arrays: Dictionary containing :class:`RpgArray` instances.
        att_global: Dictionary containing site specific global attributes
        date: Date string of file
    """

    nc = netCDF4.Dataset(file_name, "w", format="NETCDF4_CLASSIC")
    for key, dimension in dimensions.items():
        nc.createDimension(key, dimension)
    _write_vars2nc(nc, rpg_arrays)
    _add_standard_global_attributes(nc, att_global)
    return nc


def _get_dimensions(nc: netCDF4.Dataset, data: np.ndarray) -> tuple:
    """Finds correct dimensions for a variable."""

    if utils.isscalar(data):
        return ()
    variable_size = ()
    file_dims = nc.dimensions
    array_dims = data.shape
    for length in array_dims:
        dim = [key for key in file_dims.keys() if file_dims[key].size == length][0]
        variable_size = variable_size + (dim,)
    return variable_size


def _write_vars2nc(nc: netCDF4.Dataset, mwr_variables: dict) -> None:
    """Iterates over RPG instances and write to netCDF file."""

    for obj in mwr_variables.values():
        if obj.data_type == "f4":
            fill_value = -999.0
        else:
            fill_value = -99

        size = obj.dimensions or _get_dimensions(nc, obj.data)
        if obj.name == "time_bnds":
            size = ("time", "bnds")
        if obj.name == "receiver_nb":
            size = "receiver_nb"
        if obj.name == "irt":
            size = ("time", "ir_wavelength")
        if obj.name == "ir_wavelength":
            size = "ir_wavelength"
        if obj.name == "t_amb":
            size = ("time", "bnds")
        nc_variable = nc.createVariable(
            obj.name, obj.data_type, size, zlib=True, fill_value=fill_value
        )
        nc_variable[:] = obj.data
        for attr in obj.fetch_attributes():
            setattr(nc_variable, attr, getattr(obj, attr))


def _add_standard_global_attributes(nc: netCDF4.Dataset, att_global) -> None:
    nc.mwrpy_version = version.__version__
    locale.setlocale(locale.LC_TIME, "en_US.UTF-8")
    nc.processed = datetime.now(tz=timezone.utc).strftime("%d %b %Y %H:%M:%S") + " UTC"
    for name, value in att_global.items():
        if value is None:
            value = ""
        setattr(nc, name, value)
