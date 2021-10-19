from typing import Optional, Union, Tuple
import numpy as np
import netCDF4
from level1.meta_nc import MetaData
import utils

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

    def __init__(self,
                 variable: Union[netCDF4.Variable, np.ndarray, float, int],
                 name: str,
                 units_from_user: Optional[str] = None,
                 dimensions: Optional[tuple] = None):
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
            if attr not in ('name', 'data', 'data_type', 'variable', 'dimensions'):
                attributes.append(attr)
        return attributes

    def set_attributes(self, 
                       attributes: MetaData) -> None:
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
        raise ValueError(f'Incorrect RpgArray input: {self.variable}')

    def _init_units(self, 
                    units_from_user: Union[str, None]) -> str:
        if units_from_user is not None:
            return units_from_user
        return getattr(self.variable, 'units', '')

    def _init_data_type(self) -> str:
        if self.data.dtype in (np.float32, np.float64):
            return 'f4'
        return 'i4'

    def __getitem__(self, ind: tuple) -> np.ndarray:
        return self.data[ind]


class Rpg:
    """Base class for RPG MWR."""
    
    def __init__(self, raw_data: dict):
        self.raw_data = raw_data
        self.data = {}
        self.data = self._init_data()

    def _init_data(self) -> dict:
        data = {}
        for key in self.raw_data:
            data[key] = RpgArray(self.raw_data[key], key)
        return data
    
    
def save_rpg(rpg: Rpg,
             output_file: str,
             att: dict,
             data_type: str) -> Tuple[str, list]:
    """Saves the RPG MWR file."""
    
    if data_type == '1B01':
        dims = {'time': len(rpg.data['time'][:]),
                'frequency': len(rpg.data['tb'][:].T),
                'n_receivers': len(rpg.data['t_rec'][:].T)}
    elif data_type == '1B11':
        dims = {'time': len(rpg.data['time'][:]),
                'ir_wavelength': len(rpg.data['irt'][:].T)}
    elif data_type == '1B21':
        dims = {'time': len(rpg.data['time'][:])}        
    else:
        raise RuntimeError(['Data type '+ data_type +' not supported for file writing.'])

    rootgrp = init_file(output_file, dims, rpg.data, att)
    rootgrp.close()
    
    
def init_file(file_name: str,
              dimensions: dict,
              rpg_arrays: dict,
              att_global: dict) -> netCDF4.Dataset:
    """Initializes a RPG MWR file for writing.
    Args:
        file_name: File name to be generated.
        dimensions: Dictionary containing dimension for this file.
        rpg_arrays: Dictionary containing :class:`RpgArray` instances.
        att_global: Dictionary containing site specific global attributes
        
    """
    nc = netCDF4.Dataset(file_name, 'w', format='NETCDF4_CLASSIC')
    for key, dimension in dimensions.items():
        nc.createDimension(key, dimension)
    _write_vars2nc(nc, rpg_arrays)
    _add_standard_global_attributes(nc, att_global)
    return nc
    
    
def _get_dimensions(nc: netCDF4.Dataset, 
                    data: np.ndarray) -> tuple:
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


def _write_vars2nc(nc: netCDF4.Dataset, 
                   cloudnet_variables: dict) -> None:
    """Iterates over Cloudnet instances and write to netCDF file."""
    
    for obj in cloudnet_variables.values():
        if obj.data_type == 'f4':
            fill_value = -999.
        else:
            fill_value = -99

        size = obj.dimensions or _get_dimensions(nc, obj.data)
        nc_variable = nc.createVariable(obj.name, obj.data_type, size, zlib=True,
                                        fill_value=fill_value)
        nc_variable[:] = obj.data
        for attr in obj.fetch_attributes():
            setattr(nc_variable, attr, getattr(obj, attr))
            

def _add_standard_global_attributes(nc: netCDF4.Dataset, 
                                    att_global) -> None:
    for name, value in att_global.items():
        setattr(nc, name, value)