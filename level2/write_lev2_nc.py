import rpg_mwr
from level2.lev2_meta_nc import get_data_attributes
from site_config import get_site_specs
import numpy as np
import netCDF4 as nc
from typing import Optional
import glob
import rpg_mwr
from scipy.interpolate import interp1d


Fill_Value_Float = -999.
Fill_Value_Int = -99  

def lev2_to_nc(site: str,
               data_type: str,
               path_to_file: str, 
               output_file: str) -> dict:
    """This function reads Level 1 files,
    applies retrieval coefficients for Level 2 products and writes it into a netCDF file.
    
    Args:
        site: Name of site.
        data_type: Data type of the netCDF file.
        path_to_file: Folder containing Level 1 file.
        output_file: Output file name.
        
    Examples:
        >>> from level2.write_lev2_nc import lev2_to_nc
        >>> lev2_to_nc('site_name', '2P00', '/path/to/file/', 'lev2_data.nc')
    """
    
    lev1 = nc.Dataset(path_to_file)
    global_attributes, params = get_site_specs(site, data_type) 
    rpg_dat = get_products(lev1, data_type, params)
    hatpro = rpg_mwr.Rpg(rpg_dat)
    hatpro.data = get_data_attributes(hatpro.data, data_type)
    rpg_mwr.save_rpg(hatpro, output_file, global_attributes, data_type, params)
    
    
def get_products(lev1: dict, 
                 data_type: str, 
                 params: dict) -> dict:
    "Derive specified Level 2 products"
    
    rpg_dat = dict()
    lev1_vars = ['time', 'time_bounds', 'station_latitude', 'station_longitude', 'azi', 'ele']
    for ivars in lev1_vars:
        rpg_dat[ivars] = lev1.variables[ivars]
    
    if data_type == '2I00':
        
        path_coeff = '/run/user/1000/gvfs/sftp:host=fouis.meteo.uni-koeln.de,user=tmarke/home/hatpro/mwr_pro_jue/mwr_pro/retrievals/'
        offset_lwp, lin_lwp, quad_lwp, c_lwp = get_mvr_coeff(path_coeff, lev1.variables['frequency'][:], params['algo_lwp'][:])
        offset_iwv, lin_iwv, quad_iwv, c_iwv = get_mvr_coeff(path_coeff, lev1.variables['frequency'][:], params['algo_iwv'][:])        

        rpg_dat['Lwp'] = np.ones(len(lev1.variables['time'])) * Fill_Value_Float
        rpg_dat['Iwv'] = np.ones(len(lev1.variables['time'])) * Fill_Value_Float

        for itim, _ in enumerate(lev1.variables['time']):
            
            if (np.min(np.abs((lev1.variables['ele'][itim]) - c_lwp['ele'])) < .6) & (np.sum(lev1.variables['quality_flag'][itim, :]) == 0):

                offset = offset_lwp(lev1.variables['ele'][itim].data)
                coeff_lin = lin_lwp(lev1.variables['ele'][itim].data)
                coeff_quad = quad_lwp(lev1.variables['ele'][itim].data) 
                rpg_dat['Lwp'][itim] = offset + np.sum(coeff_lin * lev1.variables['tb'][itim, :]) + np.sum(coeff_quad * lev1.variables['tb'][itim, :] **2)

                offset = offset_iwv(lev1.variables['ele'][itim].data)
                coeff_lin = lin_iwv(lev1.variables['ele'][itim].data)
                coeff_quad = quad_iwv(lev1.variables['ele'][itim].data) 
                rpg_dat['Iwv'][itim] = offset + np.sum(coeff_lin * lev1.variables['tb'][itim, :]) + np.sum(coeff_quad * lev1.variables['tb'][itim, :] **2)                
                
                
#     elif data_type == '2P00':
        
        
#     elif data_type == '2S00':
        
        
    else:
        raise RuntimeError(['Data type '+ data_type +' not supported for file writing.'])
        
            
    return rpg_dat


def get_coeff_list(path_to_files: str, 
                   params: dict):
    
    c_list = []
    for i, name in enumerate(params):        
        c_list = c_list + sorted(glob.glob((path_to_files + name) + '*.nc'))
    if len(c_list) < 1:
        raise RuntimeError(['Error: no coefficient files found in directory ' + path_to_files])
    return c_list


def get_mvr_coeff(path: str,
                  freq: np.ndarray,
                  params: dict):
    "Extract retrieval coefficients for given file(s) and perform elevation angle interpolation"

    c_list = get_coeff_list(path, params)
    cf = dict()
    cf['ele'] = np.ones(len(c_list)) * Fill_Value_Float    
    cf['freq'] = np.ones([len(freq), len(c_list)]) * Fill_Value_Float
    cf['coeff_lin'] = np.zeros([len(freq), len(c_list)])
    cf['coeff_quad'] = np.zeros([len(freq), len(c_list)])
    cf['offset'] = np.zeros(len(c_list))
    cf['err'] = np.ones(len(c_list)) * Fill_Value_Float
    

    for i, file in enumerate(c_list):
        coeff = nc.Dataset(file)
        cf['ele'][i] = coeff.variables['elevation_predictand'][0]
        _, freq_ind, _ = np.intersect1d(freq[:], coeff.variables['freq'][:], assume_unique = True, return_indices = True)
        cf['freq'][freq_ind, i] = coeff.variables['freq'][freq_ind]
        cf['coeff_lin'][freq_ind, i] = coeff.variables['coefficient_mvr'][freq_ind]
        cf['coeff_quad'][freq_ind, i] = coeff.variables['coefficient_mvr'][freq_ind + len(freq_ind)]
        cf['offset'][i] = coeff.variables['offset_mvr'][0]
        cf['err'][i] = coeff.variables['predictand_err'][0]
            
    if len(c_list) > 1:            
        f_offset = interp1d(cf['ele'], cf['offset'])
        f_lin = interp1d(cf['ele'], cf['coeff_lin'])
        f_quad = interp1d(cf['ele'], cf['coeff_quad'])

    else:
        f_offset = interp1d([cf['ele'][0] - .6, cf['ele'][0] + .6], [cf['offset'][0], cf['offset'][0]])
        x = np.reshape([cf['coeff_lin'][:, 0], cf['coeff_lin'][:, 0]], (2, len(freq)))
        f_lin = interp1d([cf['ele'][0] - .6, cf['ele'][0] + .6], x.T)
        x = np.reshape([cf['coeff_quad'][:, 0], cf['coeff_quad'][:, 0]], (2, len(freq)))
        f_quad = interp1d([cf['ele'][0] - .6, cf['ele'][0] + .6], x.T)

    return f_offset, f_lin, f_quad, cf