from read_specs import get_site_specs
from level2.lev2_meta_nc import get_data_attributes
from utils import df_interp, get_coeff_list
import rpg_mwr
import numpy as np
from numpy import ma
import netCDF4 as nc
from scipy.interpolate import interp1d
import pandas as pd
from pandas.tseries.frequencies import to_offset

Fill_Value_Float = -999.
Fill_Value_Int = -99  


def lev2_to_nc(site: str,
               data_type: str,
               lev1_path: str, 
               output_file: str) -> dict:
    """This function reads Level 1 files,
    applies retrieval coefficients for Level 2 products and writes it into a netCDF file.
    
    Args:
        site: Name of site.
        data_type: Data type of the netCDF file.
        lev1_path: Path of Level 1 file.
        output_file: Output file name.
        
    Examples:
        >>> from level2.write_lev2_nc import lev2_to_nc
        >>> lev2_to_nc('site_name', '2P00', '/path/to/lev1_file/lev1_data.nc', 'lev2_data.nc')
    """
    
    global_attributes, params = get_site_specs(site, data_type)
    lev1 = nc.Dataset(lev1_path)  
    rpg_dat, att = get_products(site, lev1, data_type, params, global_attributes)
    lev1.close()
    hatpro = rpg_mwr.Rpg(rpg_dat)
    hatpro.data = get_data_attributes(hatpro.data, data_type)
    rpg_mwr.save_rpg(hatpro, output_file, att, data_type, params, site)
    
    
def get_products(site: str, 
                 lev1: dict, 
                 data_type: str, 
                 params: dict, 
                 glob_att: dict) -> dict:
    "Derive specified Level 2 products"
    
    rpg_dat = dict()
    lev1_vars = ['time', 'time_bounds', 'station_latitude', 'station_longitude', 'azi', 'ele']    
    
            
    if data_type == '2I06':
        
        offset_lwp, lin_lwp, quad_lwp, c_lwp, ran_lwp, sys_lwp, att = get_mvr_coeff(site, 'lwp', lev1['frequency'][:], glob_att)
        
        _, freq_ind, _ = np.intersect1d(lev1['frequency'][:], c_lwp['freq'][:, 0], assume_unique = False, return_indices = True)
        index = np.where((lev1['pointing_flag'][:] == 0) & (np.abs((lev1['ele'][:]) - c_lwp['ele']) < .6))[0]
        flag = np.where(np.sum(lev1['quality_flag'][:, freq_ind], axis = 1) > 0)[0]

        if any(index):
            offset = offset_lwp(lev1['ele'][index])
            coeff_lin = lin_lwp(lev1['ele'][index])
            coeff_quad = quad_lwp(lev1['ele'][index])        
            rpg_dat['Lwp'] = offset + np.sum(coeff_lin.T * lev1['tb'][index, :], axis = 1) + np.sum(coeff_quad.T * lev1['tb'][index, :] **2, axis = 1)
            rpg_dat['Lwp_random_error'] = ran_lwp(lev1['ele'][index])
            rpg_dat['Lwp_systematic_error'] = sys_lwp(lev1['ele'][index])         

            lwp_offset = get_lwp_offset(lev1['time'][index], rpg_dat['Lwp'], lev1, index, site)
            rpg_dat['Lwp'] = rpg_dat['Lwp'] - lwp_offset
            rpg_dat['Lwp_offset'] = lwp_offset
        else:
            raise RuntimeError(['No suitable data found for processing.'])             
        
        
    elif data_type == '2I07':
        
        offset_iwv, lin_iwv, quad_iwv, c_iwv, ran_iwv, sys_iwv, att = get_mvr_coeff(site, 'iwv', lev1['frequency'][:], glob_att)
        
        _, freq_ind, _ = np.intersect1d(lev1['frequency'][:], c_iwv['freq'][:, 0], assume_unique = False, return_indices = True)
        index = np.where((lev1['pointing_flag'][:] == 0) & (np.abs((lev1['ele'][:]) - c_iwv['ele']) < .6))[0]
        flag = np.where(np.sum(lev1['quality_flag'][:, freq_ind], axis = 1) > 0)[0]
        
        if any(index):
            offset = offset_iwv(lev1['ele'][index])
            coeff_lin = lin_iwv(lev1['ele'][index])
            coeff_quad = quad_iwv(lev1['ele'][index])        
            rpg_dat['Iwv'] = offset + np.sum(coeff_lin.T * lev1['tb'][index, :], axis = 1) + np.sum(coeff_quad.T * lev1['tb'][index, :] **2, axis = 1)
            rpg_dat['Iwv_random_error'] = ran_iwv(lev1['ele'][index])
            rpg_dat['Iwv_systematic_error'] = sys_iwv(lev1['ele'][index])   
        else:
            raise RuntimeError(['No suitable data found for processing.'])     
            
            
    elif data_type == '2P01':
                
        offset_tze, lin_tze, quad_tze, c_tze, ran_tze, sys_tze, att = get_mvr_coeff(site, 'tze', lev1['frequency'][:], glob_att)

        _, freq_ind, _ = np.intersect1d(lev1['frequency'][:], c_tze['freq'][:, 0], assume_unique = False, return_indices = True)
        index = np.where((lev1['pointing_flag'][:] == 0) & np.any(np.abs((np.ones((len(lev1['ele'][:]), len(c_tze['ele']))) * c_tze['ele']) - np.transpose(np.ones((len(c_tze['ele']), len(lev1['ele'][:]))) * lev1['ele'][:])) < .6, axis = 1))[0]
        flag = np.where(np.sum(lev1['quality_flag'][:, freq_ind], axis = 1) > 0)[0]
        
        if any(index):
            rpg_dat['altitude'] = c_tze['height_grid']            
            rpg_dat['temperature'] = np.ones((len(index), c_tze['n_height_grid'])) * Fill_Value_Float      

            tze_offset = offset_tze(lev1['ele'][index])
            tze_coeff_lin = lin_tze(lev1['ele'][index])
            tze_coeff_quad = quad_tze(lev1['ele'][index])               

            for ialt, _ in enumerate(c_tze['height_grid']):

                rpg_dat['temperature'][:, ialt] = tze_offset[ialt, :].T + np.sum(tze_coeff_lin[ialt, :, :].T * lev1['tb'][index, :], axis = 1) + np.sum(tze_coeff_quad[ialt, :, :].T * lev1['tb'][index, :] **2, axis = 1)    

            rpg_dat['temperature_random_error'] = ran_tze(lev1['ele'][index]).T
            rpg_dat['temperature_systematic_error'] = sys_tze(lev1['ele'][index]).T       
            
        else:
            raise RuntimeError(['No suitable data found for processing.'])                
                
                
    elif data_type == '2P02':
                
        offset_tel, lin_tel, quad_tel, c_tel, ran_tel, sys_tel, att = get_mvr_coeff(site, 'tel', lev1['frequency'][:], glob_att)

        _, freq_ind, _ = np.intersect1d(lev1['frequency'][:], c_tel['freq'], assume_unique = False, return_indices = True)
        _, freq_bl, _ = np.intersect1d(c_tel['freq'], c_tel['freq_bl'], assume_unique = False, return_indices = True)
                    
        ix0 = np.where((lev1['ele'][:] > c_tel['ele'][0] - .6) & (lev1['ele'][:] < c_tel['ele'][0] + .6) & (lev1['pointing_flag'][:] == 1))[0]

        if ix0.size > 1:
            ele = []
            flag = []
            for i in range(len(c_tel['ele'])):
                ele.append(lev1['ele'][ix0 + i])
                flag.append(np.sum(lev1['quality_flag'][ix0 + i, freq_ind], axis = 1))
            ele = np.stack(ele) * np.ones((len(c_tel['ele']), len(lev1['ele'][ix0]))) 
            flag = np.stack(flag) * np.ones((len(c_tel['ele']), len(lev1['ele'][ix0]))) 

            index = np.where(np.any(np.abs((np.ones((len(lev1['ele'][ix0]), len(c_tel['ele']))) * c_tel['ele']) - np.transpose(ele)) < .6, axis = 1))[0]
            index = ix0[index]
            flag = np.where(pd.DataFrame((flag > 0)).all())[0]

            if any(index):
                tb = np.ones((len(c_tel['freq']), len(c_tel['ele']), len(index))) * Fill_Value_Float            
                for ii, iv in enumerate(index):
                    tb[:, :, ii] = lev1['tb'][iv:iv + len(c_tel['ele']), freq_ind].T

                tb_alg = []
                if len(freq_ind) - len(freq_bl) > 0:
                    tb_alg = np.squeeze(tb[0:len(freq_ind)-len(freq_bl), 0, :])
                for ifq, _ in enumerate(c_tel['freq_bl']):
                    tb_alg = np.append(tb_alg, np.squeeze(tb[freq_bl[ifq], :, :]), axis = 0)

                rpg_dat['altitude'] = c_tel['height_grid'] 
                rpg_dat['temperature'] = np.ones((len(index), c_tel['n_height_grid'])) * Fill_Value_Float   
                for ialt, _ in enumerate(c_tel['height_grid']):
                    rpg_dat['temperature'][:, ialt] = offset_tel[ialt] + np.sum(lin_tel[:, ialt] * tb_alg[:, :].T, axis = 1)

                rpg_dat['temperature_random_error'] = np.repeat(ran_tel, len(index)).reshape((len(ran_tel), len(index))).T
                rpg_dat['temperature_systematic_error'] = np.repeat(ran_tel, len(index)).reshape((len(ran_tel), len(index))).T 

            else:
                raise RuntimeError(['No suitable data found for processing.'])   
        else:
            raise RuntimeError(['No suitable data found for processing.'])                 
                

    elif data_type == '2P03':
                
        offset_hze, lin_hze, quad_hze, c_hze, ran_hze, sys_hze, att = get_mvr_coeff(site, 'hze', lev1['frequency'][:], glob_att)

        _, freq_ind, _ = np.intersect1d(lev1['frequency'][:], c_hze['freq'][:, 0], assume_unique = False, return_indices = True)
        index = np.where((lev1['pointing_flag'][:] == 0) & np.any(np.abs((np.ones((len(lev1['ele'][:]), len(c_hze['ele']))) * c_hze['ele']) - np.transpose(np.ones((len(c_hze['ele']), len(lev1['ele'][:]))) * lev1['ele'][:])) < .6, axis = 1))[0]
        flag = np.where(np.sum(lev1['quality_flag'][:, freq_ind], axis = 1) > 0)[0]
        
        if any(index):
            rpg_dat['altitude'] = c_hze['height_grid']            
            rpg_dat['water_vapor_vmr'] = np.ones((len(index), c_hze['n_height_grid'])) * Fill_Value_Float      

            hze_offset = offset_hze(lev1['ele'][index])
            hze_coeff_lin = lin_hze(lev1['ele'][index])
            hze_coeff_quad = quad_hze(lev1['ele'][index])               

            for ialt, _ in enumerate(c_hze['height_grid']):

                rpg_dat['water_vapor_vmr'][:, ialt] = hze_offset[ialt, :].T + np.sum(hze_coeff_lin[ialt, :, :].T * lev1['tb'][index, :], axis = 1) + np.sum(hze_coeff_quad[ialt, :, :].T * lev1['tb'][index, :] **2, axis = 1)    

            rpg_dat['water_vapor_vmr_random_error'] = ran_hze(lev1['ele'][index]).T
            rpg_dat['water_vapor_vmr_systematic_error'] = sys_hze(lev1['ele'][index]).T       
            
        else:
            raise RuntimeError(['No suitable data found for processing.'])               
             
                
    elif data_type == '2P04':
        
        offset_tze, lin_tze, quad_tze, c_tze, ran_tze, sys_tze, att = get_mvr_coeff(site, 'tze', lev1['frequency'][:], glob_att)
        offset_hze, lin_hze, quad_hze, c_hze, ran_hze, sys_hze, att = get_mvr_coeff(site, 'hze', lev1['frequency'][:], glob_att)

        _, freq_ind_tze, _ = np.intersect1d(lev1['frequency'][:], c_tze['freq'][:,0], assume_unique = False, return_indices = True)
        _, freq_ind_hze, _ = np.intersect1d(lev1['frequency'][:], c_hze['freq'][:,0], assume_unique = False, return_indices = True)
        index = np.where((lev1['pointing_flag'][:] == 0) & np.any(np.abs((np.ones((len(lev1['ele'][:]), len(c_tze['ele']))) * c_tze['ele']) - np.transpose(np.ones((len(c_tze['ele']), len(lev1['ele'][:]))) * lev1['ele'][:])) < .6, axis = 1) & np.any(np.abs((np.ones((len(lev1['ele'][:]), len(c_hze['ele']))) * c_hze['ele']) - np.transpose(np.ones((len(c_hze['ele']), len(lev1['ele'][:]))) * lev1['ele'][:])) < .6, axis = 1))[0]
        flag = np.where((np.sum(lev1['quality_flag'][:, freq_ind_tze], axis = 1) > 0) | (np.sum(lev1['quality_flag'][:, freq_ind_hze], axis = 1) > 0))[0]
        
        if any(index):                     
            rpg_dat['altitude'] = c_tze['height_grid']            
            rpg_dat['temperature'] = np.ones((len(index), c_tze['n_height_grid'])) * Fill_Value_Float
            rpg_dat['water_vapor_vmr'] = np.ones((len(index), c_hze['n_height_grid'])) * Fill_Value_Float      

            tze_offset = offset_tze(lev1['ele'][index])
            tze_coeff_lin = lin_tze(lev1['ele'][index])
            tze_coeff_quad = quad_tze(lev1['ele'][index])   
            hze_offset = offset_hze(lev1['ele'][index])
            hze_coeff_lin = lin_hze(lev1['ele'][index])
            hze_coeff_quad = quad_hze(lev1['ele'][index])   

            for ialt, _ in enumerate(c_tze['height_grid']):

                rpg_dat['temperature'][:, ialt] = tze_offset[ialt, :].T + np.sum(tze_coeff_lin[ialt, :, :].T * lev1['tb'][index, :], axis = 1) + np.sum(tze_coeff_quad[ialt, :, :].T * lev1['tb'][index, :] **2, axis = 1) 
                rpg_dat['water_vapor_vmr'][:, ialt] = hze_offset[ialt, :].T + np.sum(hze_coeff_lin[ialt, :, :].T * lev1['tb'][index, :], axis = 1) + np.sum(hze_coeff_quad[ialt, :, :].T * lev1['tb'][index, :] **2, axis = 1)   

            rpg_dat['temperature_random_error'] = ran_tze(lev1['ele'][index]).T
            rpg_dat['temperature_systematic_error'] = sys_tze(lev1['ele'][index]).T               
            rpg_dat['water_vapor_vmr_random_error'] = ran_hze(lev1['ele'][index]).T
            rpg_dat['water_vapor_vmr_systematic_error'] = sys_hze(lev1['ele'][index]).T           

            rpg_dat['relative_humidity'], rpg_dat['relative_humidity_random_error'] = abshum_to_rh(rpg_dat['temperature'], rpg_dat['water_vapor_vmr'], rpg_dat['temperature_random_error'], rpg_dat['water_vapor_vmr_random_error'])
            _, rpg_dat['relative_humidity_systematic_error'] = abshum_to_rh(rpg_dat['temperature'], rpg_dat['water_vapor_vmr'], rpg_dat['temperature_systematic_error'], rpg_dat['water_vapor_vmr_systematic_error'])             
            
        else:
            raise RuntimeError(['No suitable data found for processing.'])                 
                
    else:
        raise RuntimeError(['Data type '+ data_type +' not supported for file writing.']) 
    
    "add level1 data"    
    for ivars in lev1_vars:
        if lev1[ivars].ndim > 1:
            rpg_dat[ivars] = lev1[ivars][index, :] 
        else:
            rpg_dat[ivars] = lev1[ivars][index]  
            
    "mask flagged data"        
    for ivars in rpg_dat:
        if rpg_dat[ivars].ndim > 1:
            rpg_dat[ivars][flag, :] = ma.masked
        elif (rpg_dat[ivars].ndim == 1) & (len(rpg_dat[ivars]) == len(rpg_dat['time'])):
            rpg_dat[ivars][flag] = ma.masked 
          
                    
    return rpg_dat, att


def get_mvr_coeff(site: str,
                  prefix: str,
                  freq: np.ndarray,
                  glob_att: dict):
    "Extract retrieval coefficients for given file(s) and perform elevation angle interpolation for integrated quantities"

    c_list = get_coeff_list(site, prefix)
    cf = dict()    
    
    if prefix in ('lwp', 'iwv'):    
        
        cf['ele'] = np.ones(len(c_list)) * Fill_Value_Float    
        cf['freq'] = np.ones([len(freq), len(c_list)]) * Fill_Value_Float
        cf['coeff_lin'] = np.zeros([len(freq), len(c_list)])
        cf['coeff_quad'] = np.zeros([len(freq), len(c_list)])
        cf['offset'] = np.zeros(len(c_list))
        cf['err_ran'] = np.ones(len(c_list)) * Fill_Value_Float
        cf['err_sys'] = np.ones(len(c_list)) * Fill_Value_Float

        for i, file in enumerate(c_list):
            coeff = nc.Dataset(file)
            cf['ele'][i] = coeff['elevation_predictor'][0]
            _, freq_ind, freq_cf = np.intersect1d(freq[:], coeff['freq'][:], assume_unique = False, return_indices = True)
            if any(freq_ind):
                cf['freq'][freq_ind, i] = coeff['freq'][freq_cf]
                cf['coeff_lin'][freq_ind, i] = coeff['coefficient_mvr'][freq_cf]
                cf['coeff_quad'][freq_ind, i] = coeff['coefficient_mvr'][freq_cf + len(freq_cf)]
                cf['offset'][i] = coeff['offset_mvr'][0]
                cf['err_ran'][i] = coeff['predictand_err'][0]
                cf['err_sys'][i] = coeff['predictand_err_sys'][0]
            else:
                raise RuntimeError(['Instrument and retrieval frequencies do not match.']) 

        if len(c_list) > 1:            
            f_offset = interp1d(cf['ele'], cf['offset'])
            f_lin = interp1d(cf['ele'], cf['coeff_lin'])
            f_quad = interp1d(cf['ele'], cf['coeff_quad'])
            e_ran = interp1d(cf['ele'], cf['err_ran'])
            e_sys = interp1d(cf['ele'], cf['err_sys'])

        else:
            def f_offset(x): return cf['offset']
            def f_lin(x): return cf['coeff_lin']  
            def f_quad(x): return cf['coeff_quad']     
            def e_ran(x): return np.ones(len(x)) * cf['err_ran']
            def e_sys(x): return np.ones(len(x)) * cf['err_sys']    
        
            
    elif prefix in ('tze', 'hze'):
        
        coeff = nc.Dataset(c_list[0])
        n_height_grid = coeff.dimensions['n_height_grid'].size
        n_angles = len(c_list)

        cf['ele'] = np.ones(n_angles) * Fill_Value_Float    
        cf['freq'] = np.ones([len(freq), n_angles]) * Fill_Value_Float
        cf['coeff_lin'] = np.zeros([n_height_grid, len(freq), n_angles])
        cf['coeff_quad'] = np.zeros([n_height_grid, len(freq), n_angles])
        cf['offset'] = np.zeros([n_height_grid, n_angles])
        cf['err_ran'] = np.ones([n_height_grid, n_angles]) * Fill_Value_Float  
        cf['err_sys'] = np.ones([n_height_grid, n_angles]) * Fill_Value_Float  
        cf['n_height_grid'] = n_height_grid
        cf['height_grid'] = coeff['height_grid']
        
        
        for i, file in enumerate(c_list):
            coeff = nc.Dataset(file)
            cf['ele'][i] = coeff['elevation_predictor'][i]
            _, freq_ind, freq_cf = np.intersect1d(freq[:], coeff['freq'][:], assume_unique = False, return_indices = True)
            if any(freq_ind):
                cf['freq'][freq_ind, i] = coeff['freq'][freq_cf]
                cf['coeff_lin'][:, freq_ind, i] = np.transpose(coeff['coefficient_mvr'][freq_cf, :])
                cf['coeff_quad'][:, freq_ind, i] = np.transpose(coeff['coefficient_mvr'][freq_cf + len(freq_cf), :])
                cf['offset'][:, i] = coeff['offset_mvr'][:]
                cf['err_ran'][:, i] = coeff['predictand_err'][:]  
                cf['err_sys'][:, i] = coeff['predictand_err_sys'][:]  
            else:
                raise RuntimeError(['Instrument and retrieval frequencies do not match.'])  

        f_offset = lambda x: np.transpose(np.array([cf['offset'][:, (np.abs(cf['ele']-v)).argmin()] for v in x]))
        f_lin = lambda x: np.transpose(np.array([cf['coeff_lin'][:, :, (np.abs(cf['ele']-v)).argmin()] for v in x]), (1, 2, 0))
        f_quad = lambda x: np.transpose(np.array([cf['coeff_quad'][:, :, (np.abs(cf['ele']-v)).argmin()] for v in x]), (1, 2, 0))
        e_ran = lambda x: np.transpose(np.array([cf['err_ran'][:, (np.abs(cf['ele']-v)).argmin()] for v in x]))
        e_sys = lambda x: np.transpose(np.array([cf['err_sys'][:, (np.abs(cf['ele']-v)).argmin()] for v in x]))  
        
        
    elif prefix == 'tel':
        
        coeff = nc.Dataset(c_list[0])
        _, freq_ind, freq_cf = np.intersect1d(freq[:], coeff['freq'][:], assume_unique = False, return_indices = True)
        if any(freq_ind):
            cf['ele'] = coeff['elevation_predictor'][:]
            cf['height_grid'] = coeff['height_grid']
            cf['freq'] = coeff['freq']
            cf['freq_bl'] = coeff['freq_bl']
            cf['n_height_grid'] = coeff.dimensions['n_height_grid'].size
            f_offset = coeff['offset_mvr'][:]
            f_lin = coeff['coefficient_mvr'][:, :]
            f_quad = np.nan
            e_ran =  coeff['predictand_err'][:] 
            e_sys = coeff['predictand_err_sys'][:]                  
        else:
            raise RuntimeError(['Instrument and retrieval frequencies do not match.'])  
            
            
#     elif data_type == '2S00':
        
        
    else:
        raise RuntimeError(['Data type '+ data_type +' not supported for file writing.'])

    glob_att['retrieval_type'] = coeff.regression_type
    glob_att['retrieval_elevation_angles'] = str(cf['ele'])
    glob_att['retrieval_frequencies'] = str(coeff['freq'][:])
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
                   index: np.ndarray, 
                   site: str) -> np.ndarray:
    "Correct Lwp offset using 2min standard deviation and IRT"
    
    lwp_df = pd.DataFrame({'Lwp': lwp}, index = pd.to_datetime(time, unit = 's'))
    lwp_cc = lwp_df.resample("2min", origin = 'start', closed = 'left', label = 'left').count()
    lwp_cc.index = lwp_cc.index + to_offset('1min')
    lwp_std = lwp_df.resample("2min", origin = 'start', closed = 'left', label = 'left').std()
    lwp_std.index = lwp_std.index + to_offset('1min')
    lwp_std[lwp_cc['Lwp'] < 30] = np.nan    
    lwp_mx = lwp_std.resample("20min", origin = 'start', closed = 'left', label = 'left').max()
    lwp_mx.index = lwp_mx.index + to_offset('10min')
    lwp_mx = df_interp(lwp_mx, lwp_df.index)   
    
    lwp_n = np.copy(lwp)
    if 'irt' in lev1.variables:
        irt = lev1['irt'][index, 0]
        irt_df = pd.DataFrame({'Irt': irt[:]}, index = pd.to_datetime(time, unit = 's'))
        irt_mx = irt_df.resample("20min", origin = 'start', closed = 'left', label = 'left').max()
        irt_mx.index = irt_mx.index + to_offset('10min')
        irt_mx = df_interp(irt_mx, irt_df.index) 
        lwp_n[(irt_mx['Irt'] > 233.15) & (lwp_mx['Lwp'] > .002) | np.isnan(irt_mx['Irt']) | np.isnan(lwp_mx['Lwp'])] = np.nan
    else:
        lwp_n[(lwp_mx['Lwp'] > .002) | np.isnan(lwp_mx['Lwp'])] = np.nan
        
    lwp_df = pd.DataFrame({'Lwp': lwp_n}, index = pd.to_datetime(time, unit = 's'))
    lwp_mn = lwp_df.resample("20min", origin = 'start', closed = 'left', label = 'left').mean()
    lwp_mn.index = lwp_mn.index + to_offset('10min')
    
    off = pd.read_csv('site_config/' + site + '/lwp_offset.csv', usecols = ['date', 'offset'])
    ind = np.where(lwp_mn['Lwp'].values > 0)[0]
    if ind.size > 1:
        off = off.append(pd.DataFrame({'date': time[ind[0]], 'offset': lwp_mn['Lwp'][ind[0]]}, index = {0}), ignore_index = True)
        off = off.append(pd.DataFrame({'date': time[ind[-1]], 'offset': lwp_mn['Lwp'][ind[-1]]}, index = {0}), ignore_index = True)
    elif ind.size == 1:
        off = off.append(pd.DataFrame({'date': time[ind], 'offset': lwp_mn['Lwp'][int(ind)]}, index = {0}), ignore_index = True)
    off = off.sort_values(by=['date'])
    off = off.drop_duplicates(subset=['date'])
    off.to_csv('site_config/' + site + '/lwp_offset.csv', index = False)
                
    off_ind = np.where((off['date'].values < time[0]) & (time[0] - off['date'].values < 48.*3600.))[0]
    if off_ind.size == 1:
        off_ind = np.array([int(off_ind), int(off_ind)])
    if (off_ind.size > 1) & (np.isnan(lwp_mn['Lwp'][0])):
        lwp_mn['Lwp'][0] = off['offset'][off_ind[-1]]
    off_ind = np.where((off['date'].values > time[-1]) & (off['date'].values - time[-1] < 48.*3600.))[0]
    if off_ind.size == 1:
        off_ind = np.array([int(off_ind), int(off_ind)])
    if (off_ind.size > 1) & (np.isnan(lwp_mn['Lwp'][-1])):
        lwp_mn['Lwp'][-1] = off['offset'][off_ind[0]]        

    lwp_mn = df_interp(lwp_mn, lwp_df.index)
    lwp_mn = lwp_mn.interpolate(method = 'linear')
    lwp_mn = lwp_mn.fillna(method = 'bfill')
    lwp_offset = lwp_mn['Lwp'].values
    lwp_offset[np.isnan(lwp_offset)] = 0   
    
    return lwp_offset
