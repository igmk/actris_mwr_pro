import rpg_mwr
from level2.lev2_meta_nc import get_data_attributes
import numpy as np
import netCDF4 as nc
import glob
from scipy.interpolate import interp1d
import numpy.ma as ma
import pandas as pd
import importlib
from pandas.tseries.frequencies import to_offset
from utils import df_interp

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
    
    module = importlib.import_module(f"site_config." + site)
    global_attributes, params = getattr(module, f"get_site_specs")(data_type)
    lev1 = nc.Dataset(path_to_file)  
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
            
    if data_type == '2I00':
        
        offset_lwp, lin_lwp, quad_lwp, c_lwp, ran_lwp, sys_lwp = get_mvr_coeff(params['path_coeff'] + 'lwp_', lev1.variables['frequency'][:], params['algo_lwp'][:], data_type)
        offset_iwv, lin_iwv, quad_iwv, c_iwv, ran_iwv, sys_iwv = get_mvr_coeff(params['path_coeff'] + 'iwv_', lev1.variables['frequency'][:], params['algo_iwv'][:], data_type)  
            
        index = np.where((np.abs(lev1.variables['ele'][:] - c_lwp['ele']) < .6) & (lev1.variables['pointing_flag'][:] == 1) & (np.sum(lev1.variables['quality_flag'], axis = 1) == 0) & (np.abs((lev1.variables['ele'][:]) - c_iwv['ele']) < .6))
        
        index = index[0][:]
        
        offset = offset_lwp(lev1.variables['ele'][index].data)
        coeff_lin = lin_lwp(lev1.variables['ele'][index].data)
        coeff_quad = quad_lwp(lev1.variables['ele'][index].data)        
        rpg_dat['Lwp'] = offset + np.sum(coeff_lin.T * lev1.variables['tb'][index, :], axis = 1) + np.sum(coeff_quad.T * lev1.variables['tb'][index, :] **2, axis = 1)
        rpg_dat['Lwp_random_error'] = ran_lwp(lev1.variables['ele'][index].data)
        rpg_dat['Lwp_systematic_error'] = sys_lwp(lev1.variables['ele'][index].data)       
        
        offset = offset_iwv(lev1.variables['ele'][index].data)
        coeff_lin = lin_iwv(lev1.variables['ele'][index].data)
        coeff_quad = quad_iwv(lev1.variables['ele'][index].data)        
        rpg_dat['Iwv'] = offset + np.sum(coeff_lin.T * lev1.variables['tb'][index, :], axis = 1) + np.sum(coeff_quad.T * lev1.variables['tb'][index, :] **2, axis = 1)
        rpg_dat['Iwv_random_error'] = ran_iwv(lev1.variables['ele'][index].data)
        rpg_dat['Iwv_systematic_error'] = sys_iwv(lev1.variables['ele'][index].data)   
        
        for ivars in lev1_vars:
            if lev1.variables[ivars].ndim > 1:
                rpg_dat[ivars] = lev1.variables[ivars][index, :] 
            else:
                rpg_dat[ivars] = lev1.variables[ivars][index]
                
        lwp_offset = get_lwp_offset(rpg_dat['time'], rpg_dat['Lwp'], lev1.variables['irt'][index, 0].data)
        rpg_dat['Lwp'] = rpg_dat['Lwp'] - lwp_offset
        
                
    elif data_type == '2P00':
                
        offset_tze, lin_tze, quad_tze, c_tze, ran_tze, sys_tze = get_mvr_coeff(params['path_coeff'] + 'tze_', lev1.variables['frequency'][:], params['algo_tze'][:], data_type)
        offset_hze, lin_hze, quad_hze, c_hze, ran_hze, sys_hze = get_mvr_coeff(params['path_coeff'] + 'hze_', lev1.variables['frequency'][:], params['algo_hze'][:], data_type)

        index = np.where(np.any(np.abs((np.ones((len(lev1.variables['ele'][:]), len(c_tze['ele']))) * c_tze['ele']) - np.transpose(np.ones((len(c_tze['ele']), len(lev1.variables['ele'][:]))) * lev1.variables['ele'][:])) < .6, axis = 1) & np.any(np.abs((np.ones((len(lev1.variables['ele'][:]), len(c_hze['ele']))) * c_hze['ele']) - np.transpose(np.ones((len(c_hze['ele']), len(lev1.variables['ele'][:]))) * lev1.variables['ele'][:])) < .6, axis = 1) & (np.sum(lev1.variables['quality_flag'], axis = 1) == 0))
        index = index[0][:]
                             
        rpg_dat['altitude'] = c_tze['height_grid']            
        rpg_dat['temperature'] = np.ones((len(index), c_tze['n_height_grid'])) * Fill_Value_Float
        rpg_dat['water_vapor_vmr'] = np.ones((len(index), c_tze['n_height_grid'])) * Fill_Value_Float      
            
        tze_offset = offset_tze(lev1.variables['ele'][index].data)
        tze_coeff_lin = lin_tze(lev1.variables['ele'][index].data)
        tze_coeff_quad = quad_tze(lev1.variables['ele'][index].data) 
        hze_offset = offset_hze(lev1.variables['ele'][index].data)
        hze_coeff_lin = lin_hze(lev1.variables['ele'][index].data)
        hze_coeff_quad = quad_hze(lev1.variables['ele'][index].data)         
        
        for ialt, _ in enumerate(c_tze['height_grid']):

            rpg_dat['temperature'][:, ialt] = tze_offset[ialt, :].T + np.sum(tze_coeff_lin[ialt, :, :].T * lev1.variables['tb'][index, :], axis = 1) + np.sum(tze_coeff_quad[ialt, :, :].T * lev1.variables['tb'][index, :] **2, axis = 1)
            rpg_dat['water_vapor_vmr'][:, ialt] = hze_offset[ialt, :].T + np.sum(hze_coeff_lin[ialt, :, :].T * lev1.variables['tb'][index, :], axis = 1) + np.sum(hze_coeff_quad[ialt, :, :].T * lev1.variables['tb'][index, :] **2, axis = 1)    
        
        rpg_dat['temperature_random_error'] = ran_tze(lev1.variables['ele'][index].data).T
        rpg_dat['temperature_systematic_error'] = sys_tze(lev1.variables['ele'][index].data).T    
        rpg_dat['water_vapor_vmr_random_error'] = ran_hze(lev1.variables['ele'][index].data).T
        rpg_dat['water_vapor_vmr_systematic_error'] = sys_hze(lev1.variables['ele'][index].data).T     
        rpg_dat['relative_humidity'], rpg_dat['relative_humidity_random_error'] = abshum_to_rh(rpg_dat['temperature'], rpg_dat['water_vapor_vmr'], rpg_dat['temperature_random_error'], rpg_dat['water_vapor_vmr_random_error'])
        rpg_dat['relative_humidity'], rpg_dat['relative_humidity_systematic_error'] = abshum_to_rh(rpg_dat['temperature'], rpg_dat['water_vapor_vmr'], rpg_dat['temperature_systematic_error'], rpg_dat['water_vapor_vmr_systematic_error'])          
        
        for ivars in lev1_vars:
            if lev1.variables[ivars].ndim > 1:
                rpg_dat[ivars] = lev1.variables[ivars][index, :] 
            else:
                rpg_dat[ivars] = lev1.variables[ivars][index] 
                

#     elif data_type == '2S00':
        
        
    else:
        raise RuntimeError(['Data type '+ data_type +' not supported for file writing.'])
        
            
    return rpg_dat


def get_coeff_list(path_to_files: str, 
                   params: dict):
    "Returns list of .nc coefficient file(s)"
    
    c_list = []
    for i, name in enumerate(params):        
        c_list = c_list + sorted(glob.glob((path_to_files + name) + '*.nc'))
    if len(c_list) < 1:
        raise RuntimeError(['Error: no coefficient files found in directory ' + path_to_files])
    return c_list


def get_mvr_coeff(path: str,
                  freq: np.ndarray,
                  params: dict, 
                  data_type: str):
    "Extract retrieval coefficients for given file(s) and perform elevation angle interpolation"

    c_list = get_coeff_list(path, params)
    cf = dict()
    
    if data_type == '2I00':    
        
        cf['ele'] = np.ones(len(c_list)) * Fill_Value_Float    
        cf['freq'] = np.ones([len(freq), len(c_list)]) * Fill_Value_Float
        cf['coeff_lin'] = np.zeros([len(freq), len(c_list)])
        cf['coeff_quad'] = np.zeros([len(freq), len(c_list)])
        cf['offset'] = np.zeros(len(c_list))
        cf['err_ran'] = np.ones(len(c_list)) * Fill_Value_Float
        cf['err_sys'] = np.ones(len(c_list)) * Fill_Value_Float

        for i, file in enumerate(c_list):
            coeff = nc.Dataset(file)
            cf['ele'][i] = coeff.variables['elevation_predictand'][0]
            _, freq_ind, freq_cf = np.intersect1d(freq[:], coeff.variables['freq'][:], assume_unique = True, return_indices = True)
            cf['freq'][freq_ind, i] = coeff.variables['freq'][freq_cf]
            cf['coeff_lin'][freq_ind, i] = coeff.variables['coefficient_mvr'][freq_cf]
            cf['coeff_quad'][freq_ind, i] = coeff.variables['coefficient_mvr'][freq_cf + len(freq_cf)]
            cf['offset'][i] = coeff.variables['offset_mvr'][0]
            cf['err_ran'][i] = coeff.variables['predictand_err'][0]
            cf['err_sys'][i] = coeff.variables['predictand_err_sys'][0]

        if len(c_list) > 1:            
            f_offset = interp1d(cf['ele'], cf['offset'])
            f_lin = interp1d(cf['ele'], cf['coeff_lin'])
            f_quad = interp1d(cf['ele'], cf['coeff_quad'])
            e_ran = interp1d(cf['ele'], cf['err_ran'])
            e_sys = interp1d(cf['ele'], cf['err_sys'])

        else:
            f_offset = interp1d([cf['ele'][0] - .6, cf['ele'][0] + .6], [cf['offset'][0], cf['offset'][0]])
            x = np.reshape([cf['coeff_lin'][:, 0], cf['coeff_lin'][:, 0]], (2, len(freq)))
            f_lin = interp1d([cf['ele'][0] - .6, cf['ele'][0] + .6], x.T)
            x = np.reshape([cf['coeff_quad'][:, 0], cf['coeff_quad'][:, 0]], (2, len(freq)))
            f_quad = interp1d([cf['ele'][0] - .6, cf['ele'][0] + .6], x.T)
            e_ran = interp1d([cf['ele'][0] - .6, cf['ele'][0] + .6], [cf['err_ran'][0], cf['err_ran'][0]])
            e_sys = interp1d([cf['ele'][0] - .6, cf['ele'][0] + .6], [cf['err_sys'][0], cf['err_sys'][0]])
            
            
    elif data_type == '2P00':
        
        coeff = nc.Dataset(c_list[0])
        n_height_grid = len(coeff.variables['height_grid'])  

        cf['ele'] = np.ones(len(c_list)) * Fill_Value_Float    
        cf['freq'] = np.ones([len(freq), len(c_list)]) * Fill_Value_Float
        cf['coeff_lin'] = np.zeros([n_height_grid, len(freq), len(c_list)])
        cf['coeff_quad'] = np.zeros([n_height_grid, len(freq), len(c_list)])
        cf['offset'] = np.zeros([n_height_grid, len(c_list)])
        cf['err_ran'] = np.ones([n_height_grid, len(c_list)]) * Fill_Value_Float  
        cf['err_sys'] = np.ones([n_height_grid, len(c_list)]) * Fill_Value_Float  
        cf['n_height_grid'] = n_height_grid
        cf['height_grid'] = coeff.variables['height_grid']
        
        for i, file in enumerate(c_list):
            coeff = nc.Dataset(file)
            cf['ele'][i] = coeff.variables['elevation_predictand'][0]
            _, freq_ind, freq_cf = np.intersect1d(freq[:], coeff.variables['freq'][:], assume_unique = True, return_indices = True)
            cf['freq'][freq_ind, i] = coeff.variables['freq'][freq_cf]
            cf['coeff_lin'][:, freq_ind, i] = np.transpose(coeff.variables['coefficient_mvr'][freq_cf, :])
            cf['coeff_quad'][:, freq_ind, i] = np.transpose(coeff.variables['coefficient_mvr'][freq_cf + len(freq_cf), :])
            cf['offset'][:, i] = coeff.variables['offset_mvr'][:]
            cf['err_ran'][:, i] = coeff.variables['predictand_err'][:]  
            cf['err_sys'][:, i] = coeff.variables['predictand_err_sys'][:]  
            
        if len(c_list) > 1:            
            f_offset = interp1d(cf['ele'], cf['offset'])
            f_lin = interp1d(cf['ele'], cf['coeff_lin'])
            f_quad = interp1d(cf['ele'], cf['coeff_quad'])   
            e_ran = interp1d(cf['ele'], cf['err_ran'])
            e_sys = interp1d(cf['ele'], cf['err_sys'])
            
        else:
            raise RuntimeError(['More than 1 coefficient file needed for data type '+ data_type ])
        
        
#     elif data_type == '2S00':
        
        
    else:
        raise RuntimeError(['Data type '+ data_type +' not supported for file writing.'])        

    return f_offset, f_lin, f_quad, cf, e_ran, e_sys



def abshum_to_rh(T: np.ndarray, 
                 q: np.ndarray,
                 dT: np.ndarray,
                 dq: np.ndarray):
    "Calculate relative humidity and its error from absolute humidity and temperature using analytical approximation of Clausius-Clapeyron equation (+ error propagation)"
    
    "specific gas constant for water vapor (J/kg K)"
    Rw = 462.
    "vapor pressure e0 (Pa) at T0 (K)"
    T0 = 273.15
    e0 = 611.
    "specific heat for evaporation (J/kg)"
    L = (2500.-2.42*(T-T0))*1000.    
    "saturation pressure in Pa "
    es = e0*np.exp((L/(Rw*T0))*((T-T0)/T))
    "water vapor pressure"
    e = q*Rw*T 
    "relative humidity"
    rh = e/es
    
    "error propagation"
    drh_dq = Rw*T/es
    des_dT = es*(T*(T-T0)*(-2420)+T0*L)/(Rw*T0*T**2)   
    drh_dT = q*Rw/es**2*(es-T*des_dT)
    drh = np.sqrt((drh_dq * dq)**2 + (drh_dT * dT)**2)
    
    return rh, drh
    
    
def get_lwp_offset(time: np.ndarray, 
                   lwp: np.ndarray,
                   irt: np.ndarray) -> np.ndarray:
    "Correct Lwp using 2min standard deviation and IRT"

    lwp_df = pd.DataFrame({'Lwp': lwp}, index = pd.to_datetime(time, unit = 's'))
    lwp_std = lwp_df.resample("2min").std()
    lwp_mx = lwp_std.resample("20min").max()
    loffset = '10min'
    lwp_mx.index = lwp_mx.index + to_offset(loffset)
    lwp_mx = df_interp(lwp_mx, lwp_df.index)

    irt_df = pd.DataFrame({'Irt': irt[:]}, index = pd.to_datetime(time, unit = 's'))
    irt_mn = irt_df.resample("2min").mean()
    irt_mx = irt_mn.resample("20min").max()
    loffset = '10min'
    irt_mx.index = irt_mx.index + to_offset(loffset)
    irt_mx = df_interp(irt_mx, irt_df.index)    
    
    cld = np.ones(len(lwp))
    cld[(irt_mx['Irt'].values > 233.15) & (lwp_mx['Lwp'].values > .002)] = 0
    lwp_n = np.copy(lwp)
    lwp_n[cld == 0] = np.nan
    lwp_df = pd.DataFrame({'Lwp': lwp_n}, index = pd.to_datetime(time, unit = 's'))
    lwp_mn = lwp_df.resample("20min").mean()
    loffset = '10min'
    lwp_mn.index = lwp_mn.index + to_offset(loffset)
    lwp_mn = df_interp(lwp_mn, lwp_df.index)
    lwp_mn = lwp_mn.interpolate(method = 'linear')
    lwp_mn = lwp_mn.fillna(method = 'bfill')
    lwp_offset = lwp_mn['Lwp'].values
    lwp_offset[np.isnan(lwp_offset)] = 0
    
    return lwp_offset
