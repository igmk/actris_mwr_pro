from site_config.read_specs import get_site_specs
from level2.lev2_meta_nc import get_data_attributes
import rpg_mwr
import numpy as np
import netCDF4 as nc
import glob
from scipy.interpolate import interp1d
import numpy.ma as ma
import pandas as pd
from utils import df_interp
from typing import Optional

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
    
    global_attributes, params = get_site_specs(site, data_type)
    lev1 = nc.Dataset(path_to_file)  
    rpg_dat, att = get_products(lev1, data_type, params, global_attributes)
    hatpro = rpg_mwr.Rpg(rpg_dat)
    hatpro.data = get_data_attributes(hatpro.data, data_type)
    rpg_mwr.save_rpg(hatpro, output_file, att, data_type, params)
    
    
def get_products(lev1: dict, 
                 data_type: str, 
                 params: dict, 
                 glob_att: dict) -> dict:
    "Derive specified Level 2 products"
    
    rpg_dat = dict()
    lev1_vars = ['time', 'time_bounds', 'station_latitude', 'station_longitude', 'azi', 'ele']    
    
            
    if data_type == '2I06':
        
        offset_lwp, lin_lwp, quad_lwp, c_lwp, ran_lwp, sys_lwp, att = get_mvr_coeff(params['path_coeff'] + 'lwp_', lev1.variables['frequency'][:], params['algo_lwp'][:], data_type, glob_att)

        index = np.where((np.abs(lev1.variables['ele'][:] - c_lwp['ele']) < .6) & (lev1.variables['pointing_flag'][:] == 1) & np.logical_or(np.equal(np.sum(lev1.variables['quality_flag'], axis = 1), 0), np.equal(np.sum(lev1.variables['quality_flag'], axis = 1), 32)))
        index = index[0][:]
        if any(index):
            offset = offset_lwp(lev1.variables['ele'][index].data)
            coeff_lin = lin_lwp(lev1.variables['ele'][index].data)
            coeff_quad = quad_lwp(lev1.variables['ele'][index].data)        
            rpg_dat['Lwp'] = offset + np.sum(coeff_lin.T * lev1.variables['tb'][index, :], axis = 1) + np.sum(coeff_quad.T * lev1.variables['tb'][index, :] **2, axis = 1)
            rpg_dat['Lwp_random_error'] = ran_lwp(lev1.variables['ele'][index].data)
            rpg_dat['Lwp_systematic_error'] = sys_lwp(lev1.variables['ele'][index].data)         

            lwp_offset = get_lwp_offset(lev1.variables['time'][index], rpg_dat['Lwp'], lev1, index)
            rpg_dat['Lwp'] = rpg_dat['Lwp'] - lwp_offset
            rpg_dat['Lwp_off'] = lwp_offset
        
        
    elif data_type == '2I07':
        
        offset_iwv, lin_iwv, quad_iwv, c_iwv, ran_iwv, sys_iwv, att = get_mvr_coeff(params['path_coeff'] + 'iwv_', lev1.variables['frequency'][:], params['algo_iwv'][:], data_type, glob_att)
        
        index = np.where((lev1.variables['pointing_flag'][:] == 1) & np.logical_or(np.equal(np.sum(lev1.variables['quality_flag'], axis = 1), 0), np.equal(np.sum(lev1.variables['quality_flag'], axis = 1), 32)) & (np.abs((lev1.variables['ele'][:]) - c_iwv['ele']) < .6))
        index = index[0][:]
        if any(index):
            offset = offset_iwv(lev1.variables['ele'][index].data)
            coeff_lin = lin_iwv(lev1.variables['ele'][index].data)
            coeff_quad = quad_iwv(lev1.variables['ele'][index].data)        
            rpg_dat['Iwv'] = offset + np.sum(coeff_lin.T * lev1.variables['tb'][index, :], axis = 1) + np.sum(coeff_quad.T * lev1.variables['tb'][index, :] **2, axis = 1)
            rpg_dat['Iwv_random_error'] = ran_iwv(lev1.variables['ele'][index].data)
            rpg_dat['Iwv_systematic_error'] = sys_iwv(lev1.variables['ele'][index].data)          
                
                
    elif data_type == '2P02':
                
        offset_tze, lin_tze, quad_tze, c_tze, ran_tze, sys_tze, att = get_mvr_coeff(params['path_coeff'] + 'tze_', lev1.variables['frequency'][:], params['algo_tze'][:], data_type, glob_att)

        index = np.where(np.any(np.abs((np.ones((len(lev1.variables['ele'][:]), len(c_tze['ele']))) * c_tze['ele']) - np.transpose(np.ones((len(c_tze['ele']), len(lev1.variables['ele'][:]))) * lev1.variables['ele'][:])) < .6, axis = 1) & (lev1.variables['pointing_flag'][:] == 0) & np.logical_or(np.equal(np.sum(lev1.variables['quality_flag'], axis = 1), 0), np.equal(np.sum(lev1.variables['quality_flag'], axis = 1), 32)))
        index = index[0][:]
        if any(index):                     
            rpg_dat['altitude'] = c_tze['height_grid']            
            rpg_dat['temperature'] = np.ones((len(index), c_tze['n_height_grid'])) * Fill_Value_Float

            tze_offset = offset_tze(lev1.variables['ele'][index].data)
            tze_coeff_lin = lin_tze(lev1.variables['ele'][index].data)
            tze_coeff_quad = quad_tze(lev1.variables['ele'][index].data)        

            for ialt, _ in enumerate(c_tze['height_grid']):

                rpg_dat['temperature'][:, ialt] = tze_offset[ialt, :].T + np.sum(tze_coeff_lin[ialt, :, :].T * lev1.variables['tb'][index, :], axis = 1) + np.sum(tze_coeff_quad[ialt, :, :].T * lev1.variables['tb'][index, :] **2, axis = 1) 

            rpg_dat['temperature_random_error'] = ran_tze(lev1.variables['ele'][index].data).T
            rpg_dat['temperature_systematic_error'] = sys_tze(lev1.variables['ele'][index].data).T            
                

    elif data_type == '2P03':
                
        offset_hze, lin_hze, quad_hze, c_hze, ran_hze, sys_hze, att = get_mvr_coeff(params['path_coeff'] + 'hze_', lev1.variables['frequency'][:], params['algo_hze'][:], data_type, glob_att)

        index = np.where(np.any(np.abs((np.ones((len(lev1.variables['ele'][:]), len(c_hze['ele']))) * c_hze['ele']) - np.transpose(np.ones((len(c_hze['ele']), len(lev1.variables['ele'][:]))) * lev1.variables['ele'][:])) < .6, axis = 1) & np.logical_or(np.equal(np.sum(lev1.variables['quality_flag'], axis = 1), 0), np.equal(np.sum(lev1.variables['quality_flag'], axis = 1), 32)))
        index = index[0][:]
        if any(index):                     
            rpg_dat['altitude'] = c_hze['height_grid']            
            rpg_dat['water_vapor_vmr'] = np.ones((len(index), c_hze['n_height_grid'])) * Fill_Value_Float      

            hze_offset = offset_hze(lev1.variables['ele'][index].data)
            hze_coeff_lin = lin_hze(lev1.variables['ele'][index].data)
            hze_coeff_quad = quad_hze(lev1.variables['ele'][index].data)         

            for ialt, _ in enumerate(c_hze['height_grid']):

                rpg_dat['water_vapor_vmr'][:, ialt] = hze_offset[ialt, :].T + np.sum(hze_coeff_lin[ialt, :, :].T * lev1.variables['tb'][index, :], axis = 1) + np.sum(hze_coeff_quad[ialt, :, :].T * lev1.variables['tb'][index, :] **2, axis = 1)    

            rpg_dat['water_vapor_vmr_random_error'] = ran_hze(lev1.variables['ele'][index].data).T
            rpg_dat['water_vapor_vmr_systematic_error'] = sys_hze(lev1.variables['ele'][index].data).T       
             
                
    elif data_type == '2P04':
        
        offset_tze, lin_tze, quad_tze, c_tze, ran_tze, sys_tze, att = get_mvr_coeff(params['path_coeff'] + 'tze_', lev1.variables['frequency'][:], params['algo_tze'][:], data_type, glob_att)
        offset_hze, lin_hze, quad_hze, c_hze, ran_hze, sys_hze, att = get_mvr_coeff(params['path_coeff'] + 'hze_', lev1.variables['frequency'][:], params['algo_hze'][:], data_type, glob_att)

        index = np.where(np.any(np.abs((np.ones((len(lev1.variables['ele'][:]), len(c_tze['ele']))) * c_tze['ele']) - np.transpose(np.ones((len(c_tze['ele']), len(lev1.variables['ele'][:]))) * lev1.variables['ele'][:])) < .6, axis = 1) & (lev1.variables['pointing_flag'][:] == 0) & np.logical_or(np.equal(np.sum(lev1.variables['quality_flag'], axis = 1), 0), np.equal(np.sum(lev1.variables['quality_flag'], axis = 1), 32)) & np.any(np.abs((np.ones((len(lev1.variables['ele'][:]), len(c_hze['ele']))) * c_hze['ele']) - np.transpose(np.ones((len(c_hze['ele']), len(lev1.variables['ele'][:]))) * lev1.variables['ele'][:])) < .6, axis = 1))
        index = index[0][:]
        if any(index):                     
            rpg_dat['altitude'] = c_tze['height_grid']            
            rpg_dat['temperature'] = np.ones((len(index), c_tze['n_height_grid'])) * Fill_Value_Float
            rpg_dat['water_vapor_vmr'] = np.ones((len(index), c_hze['n_height_grid'])) * Fill_Value_Float      

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
                
                
    else:
        raise RuntimeError(['Data type '+ data_type +' not supported for file writing.']) 
        
    for ivars in lev1_vars:
        if lev1.variables[ivars].ndim > 1:
            rpg_dat[ivars] = lev1.variables[ivars][index, :] 
        else:
            rpg_dat[ivars] = lev1.variables[ivars][index]           
                    
    return rpg_dat, att


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
                  data_type: str, 
                  glob_att: dict):
    "Extract retrieval coefficients for given file(s) and perform elevation angle interpolation"

    c_list = get_coeff_list(path, params)
    cf = dict()
    
    if data_type in ('2I06', '2I07'):    
        
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
            
            
    elif data_type in ('2P02', '2P03', '2P04'):
        
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

    glob_att['retrieval_type'] = coeff.regression_type
    glob_att['retrieval_elevation_angles'] = str(cf['ele'] )
    glob_att['retrieval_frequencies'] = str(coeff['freq'][:].data)
    glob_att['retrieval_auxiliary_input'] = coeff.surface_mode
    glob_att['retrieval_description'] = coeff.retrieval_version

    return f_offset, f_lin, f_quad, cf, e_ran, e_sys, glob_att


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
                   lev1: dict,
                   index: np.ndarray) -> np.ndarray:
    "Correct Lwp offset using 2min standard deviation and IRT"
    
    lwp_df = pd.DataFrame({'Lwp': lwp}, index = pd.to_datetime(time, unit = 's'))
    lwp_std = lwp_df.resample("2min", origin = 'start', closed = 'left', offset = '1min', label = 'left').std()
    lwp_mx = lwp_std.resample("20min", origin = 'start', closed = 'left', offset = '10min', label = 'left').max()
    lwp_mx = df_interp(lwp_mx, lwp_df.index)   
    
    lwp_n = np.copy(lwp)
    if 'Irt' in lev1.variables:
        irt = lev1.variables['irt'][index, 0].data
        irt_df = pd.DataFrame({'Irt': irt[:]}, index = pd.to_datetime(time, unit = 's'))
        irt_mx = irt_df.resample("20min", origin = 'start', closed = 'left', offset = '10min', label = 'left').max()
        irt_mx = df_interp(irt_mx, irt_df.index) 
        lwp_n[(irt_mx['Irt'] > 233.15) & (lwp_mx['Lwp'] > .002) | np.isnan(irt_mx['Irt']) | np.isnan(lwp_mx['Lwp'])] = np.nan
    else:
        lwp_n[(lwp_mx['Lwp'] > .002) | np.isnan(lwp_mx['Lwp'])] = np.nan
    lwp_df = pd.DataFrame({'Lwp': lwp_n}, index = pd.to_datetime(time, unit = 's'))
    lwp_mn = lwp_df.resample("20min", origin = 'start', closed = 'left', offset = '10min', label = 'left').mean()
    lwp_mn = df_interp(lwp_mn, lwp_df.index)
    lwp_mn = lwp_mn.interpolate(method = 'linear')
    lwp_mn = lwp_mn.fillna(method = 'bfill')
    lwp_offset = lwp_mn['Lwp'].values
    lwp_offset[np.isnan(lwp_offset)] = 0
    
    return lwp_offset
