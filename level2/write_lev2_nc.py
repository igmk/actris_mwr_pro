from read_specs import get_site_specs
from level2.lev2_meta_nc import get_data_attributes
from level1.quality_control import spectral_consistency
from level2.get_ret_coeff import get_mvr_coeff
from level2.lwp_offset import correct_lwp_offset
from utils import rebin_2d, get_coeff_list
import rpg_mwr
import numpy as np
from numpy import ma
import netCDF4 as nc
import pandas as pd
import metpy.calc as mpcalc
from metpy.units import units, masked_array

Fill_Value_Float = -999.
Fill_Value_Int = -99  
"""specific gas constant for water vapor (J/kg K)
and vapor pressure e0 (Pa) at T0 (K)"""
Rw, e0, T0 = 462., 611., 273.15


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
    
    if data_type not in ('2P01', '2P02', '2P03', '2P04', '2P07', '2P08', '2I01', '2I02', '2S02'):
        raise RuntimeError(['Data type '+ data_type +' not supported for file writing.'])
    global_attributes, params = get_site_specs(site, data_type)
    with nc.Dataset(lev1_path) as lev1:
        rpg_dat, coeff, index, flag = get_products(site, lev1, data_type, params)
        _combine_lev1(lev1, rpg_dat, index)
        _mask_flag(rpg_dat, flag)
    if data_type in ('2P01', '2P02', '2P03', '2I01', '2I02', '2S02'):
        _add_att(global_attributes, coeff)
    hatpro = rpg_mwr.Rpg(rpg_dat)
    hatpro.data = get_data_attributes(hatpro.data, data_type)
    rpg_mwr.save_rpg(hatpro, output_file, global_attributes, data_type, params)
    
    
def get_products(site: str, 
                 lev1: dict, 
                 data_type: str, 
                 params: dict) -> dict:
    "Derive specified Level 2 products"
    
    rpg_dat = dict()  
    
            
    if data_type == '2I01':
        
        offset, lin, quad, coeff, ran_err, sys_err = get_mvr_coeff(site, 'lwp', lev1['frequency'][:])
        
        _, freq_ind, _ = np.intersect1d(lev1['frequency'][:], coeff['freq'][:, 0], assume_unique = False, return_indices = True)
        index = np.where((lev1['pointing_flag'][:] == 0) & (np.abs((lev1['ele'][:]) - coeff['ele']) < .6))[0]
        flag = np.where(np.sum(lev1['quality_flag'][index, freq_ind], axis = 1) > 0)[0]

        if len(index) == 0:
            raise RuntimeError(['No suitable data found for processing for data type: ' + data_type])
            
        coeff_offset = offset(lev1['ele'][index])
        coeff_lin = lin(lev1['ele'][index])
        coeff_quad = quad(lev1['ele'][index])        
        lwp = coeff_offset + np.sum(coeff_lin.T * lev1['tb'][index, :], axis = 1) + np.sum(coeff_quad.T * lev1['tb'][index, :] **2, axis = 1)
        rpg_dat['lwp_random_error'] = ran_err(lev1['ele'][index])
        rpg_dat['lwp_systematic_error'] = sys_err(lev1['ele'][index])
        freq_31 = np.where(lev1['frequency'][:] == 31.4)[0]
        if len(freq_31) != 1:
            rpg_dat['lwp'], rpg_dat['lwp_offset'] = lwp, np.ones(len(index)) * Fill_Value_Float
        else:
            rpg_dat['lwp'], rpg_dat['lwp_offset'] = correct_lwp_offset(lev1.variables, lwp, index, site)        
        
        
    elif data_type == '2I02':
        
        offset, lin, quad, coeff, ran_err, sys_err = get_mvr_coeff(site, 'iwv', lev1['frequency'][:])
        
        _, freq_ind, _ = np.intersect1d(lev1['frequency'][:], coeff['freq'][:, 0], assume_unique = False, return_indices = True)
        index = np.where((lev1['pointing_flag'][:] == 0) & (np.abs((lev1['ele'][:]) - coeff['ele']) < .6))[0]
        flag = np.where(np.sum(lev1['quality_flag'][index, freq_ind], axis = 1) > 0)[0]
        
        if len(index) == 0:
            raise RuntimeError(['No suitable data found for processing for data type: ' + data_type])
            
        coeff_offset = offset(lev1['ele'][index])
        coeff_lin = lin(lev1['ele'][index])
        coeff_quad = quad(lev1['ele'][index])        
        rpg_dat['iwv'] = coeff_offset + np.sum(coeff_lin.T * lev1['tb'][index, :], axis = 1) + np.sum(coeff_quad.T * lev1['tb'][index, :] **2, axis = 1)
        rpg_dat['iwv_random_error'] = ran_err(lev1['ele'][index])
        rpg_dat['iwv_systematic_error'] = sys_err(lev1['ele'][index])        
            
            
    elif data_type == '2P01':
                
        offset, lin, quad, coeff, ran_err, sys_err = get_mvr_coeff(site, 'tze', lev1['frequency'][:])

        _, freq_ind, _ = np.intersect1d(lev1['frequency'][:], coeff['freq'][:, 0], assume_unique = False, return_indices = True)
        index = np.where((lev1['pointing_flag'][:] == 0) & np.any(np.abs((np.ones((len(lev1['ele'][:]), len(coeff['ele']))) * coeff['ele']) - np.transpose(np.ones((len(coeff['ele']), len(lev1['ele'][:]))) * lev1['ele'][:])) < .6, axis = 1))[0]
        flag = np.where(np.sum(lev1['quality_flag'][index, freq_ind], axis = 1) > 0)[0]
        
        if len(index) == 0:
            raise RuntimeError(['No suitable data found for processing for data type: ' + data_type])
            
        rpg_dat['altitude'] = coeff['height_grid']            
        rpg_dat['temperature'] = ma.masked_all((len(index), coeff['n_height_grid']))     

        coeff_offset = offset(lev1['ele'][index])
        coeff_lin = lin(lev1['ele'][index])
        coeff_quad = quad(lev1['ele'][index])               

        for ialt, _ in enumerate(rpg_dat['altitude']):
            rpg_dat['temperature'][:, ialt] = coeff_offset[:, ialt] + np.sum(coeff_lin[:, ialt, :] * lev1['tb'][index, :], axis=1) + np.sum(coeff_quad[:, ialt, :] * lev1['tb'][index, :]**2, axis=1)    
        rpg_dat['temperature_random_error'] = ran_err(lev1['ele'][index])
        rpg_dat['temperature_systematic_error'] = sys_err(lev1['ele'][index])
                
                
    elif data_type == '2P02':
                
        offset, lin, quad, coeff, ran_err, sys_err = get_mvr_coeff(site, 'tel', lev1['frequency'][:])

        _, freq_ind, _ = np.intersect1d(lev1['frequency'][:], coeff['freq'], assume_unique = False, return_indices = True)
        _, freq_bl, _ = np.intersect1d(coeff['freq'], coeff['freq_bl'], assume_unique = False, return_indices = True)
                    
        ix0 = np.where((lev1['ele'][:] > coeff['ele'][0] - .6) & (lev1['ele'][:] < coeff['ele'][0] + .6) & (lev1['pointing_flag'][:] == 1))[0]
        if len(ix0) == 0:
            raise RuntimeError(['No suitable data found for processing for data type: ' + data_type])
            
        ele, flag = [], []
        for i in range(len(coeff['ele'])):
            ele.append(lev1['ele'][ix0 + i])
            flag.append(np.sum(lev1['quality_flag'][ix0 + i, freq_ind], axis = 1))
        ele = np.stack(ele) * np.ones((len(coeff['ele']), len(lev1['ele'][ix0]))) 
        flag = np.stack(flag) * np.ones((len(coeff['ele']), len(lev1['ele'][ix0]))) 

        index = np.where(np.any(np.abs((np.ones((len(lev1['ele'][ix0]), len(coeff['ele']))) * coeff['ele']) - np.transpose(ele)) < .6, axis = 1))[0]
        index = ix0[index]
        flag = np.where(pd.DataFrame((flag > 0)).all())[0]

        if len(index) == 0:
            raise RuntimeError(['No suitable data found for processing for data type: ' + data_type])

        tb = np.ones((len(coeff['freq']), len(coeff['ele']), len(index))) * Fill_Value_Float            
        for ii, iv in enumerate(index):
            tb[:, :, ii] = lev1['tb'][iv:iv + len(coeff['ele']), freq_ind].T

        tb_alg = []
        if len(freq_ind) - len(freq_bl) > 0:
            tb_alg = np.squeeze(tb[0:len(freq_ind)-len(freq_bl), 0, :])
        for ifq, _ in enumerate(coeff['freq_bl']):
            tb_alg = np.append(tb_alg, np.squeeze(tb[freq_bl[ifq], :, :]), axis = 0)

        rpg_dat['altitude'] = coeff['height_grid'] 
        rpg_dat['temperature'] = ma.masked_all((len(index), coeff['n_height_grid']))  
        for ialt, _ in enumerate(rpg_dat['altitude']):
            rpg_dat['temperature'][:, ialt] = offset[ialt] + np.sum(lin[:, ialt] * tb_alg[:, :].T, axis = 1)

        rpg_dat['temperature_random_error'] = np.repeat(ran_err, len(index)).reshape((len(ran_err), len(index))).T
        rpg_dat['temperature_systematic_error'] = np.repeat(sys_err, len(index)).reshape((len(sys_err), len(index))).T             
                

    elif data_type == '2P03':
                
        offset, lin, quad, coeff, ran_err, sys_err = get_mvr_coeff(site, 'hze', lev1['frequency'][:])

        _, freq_ind, _ = np.intersect1d(lev1['frequency'][:], coeff['freq'][:, 0], assume_unique = False, return_indices = True)
        index = np.where((lev1['pointing_flag'][:] == 0) & np.any(np.abs((np.ones((len(lev1['ele'][:]), len(coeff['ele']))) * coeff['ele']) - np.transpose(np.ones((len(coeff['ele']), len(lev1['ele'][:]))) * lev1['ele'][:])) < .6, axis = 1))[0]
        flag = np.where(np.sum(lev1['quality_flag'][index, freq_ind], axis = 1) > 0)[0]
        
        if (len(index) == 0 | len(freq_ind) == 0):
            raise RuntimeError(['No suitable data found for processing for data type: ' + data_type])
            
        rpg_dat['altitude'] = coeff['height_grid']            
        rpg_dat['water_vapor_vmr'] = ma.masked_all((len(index), coeff['n_height_grid']))  
        
        coeff_offset = offset(lev1['ele'][index])
        coeff_lin = lin(lev1['ele'][index])
        coeff_quad = quad(lev1['ele'][index])  

        for ialt, _ in enumerate(rpg_dat['altitude']):
            rpg_dat['water_vapor_vmr'][:, ialt] = coeff_offset[:, ialt] + np.sum(coeff_lin[:, ialt, :] * lev1['tb'][index, :], axis=1) + np.sum(coeff_quad[:, ialt, :] * lev1['tb'][index, :]**2, axis=1)    
        rpg_dat['water_vapor_vmr_random_error'] = ran_err(lev1['ele'][index])
        rpg_dat['water_vapor_vmr_systematic_error'] = sys_err(lev1['ele'][index])
             
                
    elif data_type == '2P04':
                
        tem_dat, coeff, index, flag = get_products(site, lev1, '2P02', params)
        hum_dat, _, ix, flg = get_products(site, lev1, '2P03', params) 
        for ivars in hum_dat:
            if hum_dat[ivars].ndim > 1:
                hum_dat[ivars][flg, :] = ma.masked
        hum = rebin_2d(lev1['time'][ix], hum_dat['water_vapor_vmr'][:,:], lev1['time'][index])
        hum_re = rebin_2d(lev1['time'][ix], hum_dat['water_vapor_vmr_random_error'][:,:], lev1['time'][index])
        hum_se = rebin_2d(lev1['time'][ix], hum_dat['water_vapor_vmr_systematic_error'][:,:], lev1['time'][index])        

        rpg_dat['altitude'] = tem_dat['altitude']
        rpg_dat['relative_humidity'] = vap_pres(hum[0], tem_dat['temperature'].data) / mpcalc.saturation_vapor_pressure(masked_array(tem_dat['temperature'], data_units='K')).magnitude
        rpg_dat['relative_humidity_random_error'] = rh_err(tem_dat['temperature'], hum[0], tem_dat['temperature_random_error'], hum_re[0])
        rpg_dat['relative_humidity_systematic_error'] = rh_err(tem_dat['temperature'], hum[0], tem_dat['temperature_systematic_error'], hum_se[0])     
        
        
    elif data_type == '2P07':
        
        tem_dat, coeff, index, flag = get_products(site, lev1, '2P02', params)      
        hum_dat, _, ix, flg = get_products(site, lev1, '2P03', params)   
        for ivars in hum_dat:
            if hum_dat[ivars].ndim > 1:
                hum_dat[ivars][flg, :] = ma.masked        
        hum = rebin_2d(lev1['time'][ix], hum_dat['water_vapor_vmr'][:,:], lev1['time'][index])
        hum_re = rebin_2d(lev1['time'][ix], hum_dat['water_vapor_vmr_random_error'][:,:], lev1['time'][index])
        hum_se = rebin_2d(lev1['time'][ix], hum_dat['water_vapor_vmr_systematic_error'][:,:], lev1['time'][index]) 

        rpg_dat['altitude'] = tem_dat['altitude']
        p_baro = calc_p_baro(tem_dat['temperature'], hum[0], lev1['air_pressure'][index], rpg_dat['altitude'])
        rpg_dat['potential_temperature'] = mpcalc.potential_temperature(masked_array(p_baro, data_units='Pa'), masked_array(tem_dat['temperature'], data_units='K')).magnitude
            
    elif data_type == '2P08':
        
        tem_dat, coeff, index, flag = get_products(site, lev1, '2P02', params)
        hum_dat, _, ix, flg = get_products(site, lev1, '2P03', params)   
        for ivars in hum_dat:
            if hum_dat[ivars].ndim > 1:
                hum_dat[ivars][flg, :] = ma.masked   
        hum = rebin_2d(lev1['time'][ix], hum_dat['water_vapor_vmr'][:,:], lev1['time'][index])
        hum_re = rebin_2d(lev1['time'][ix], hum_dat['water_vapor_vmr_random_error'][:,:], lev1['time'][index])
        hum_se = rebin_2d(lev1['time'][ix], hum_dat['water_vapor_vmr_systematic_error'][:,:], lev1['time'][index])        
        
        rpg_dat['altitude'] = tem_dat['altitude']
        p_baro = calc_p_baro(tem_dat['temperature'], hum[0], lev1['air_pressure'][index], rpg_dat['altitude'])
        Theta = mpcalc.potential_temperature(masked_array(p_baro, data_units='Pa'), masked_array(tem_dat['temperature'], data_units='K')).magnitude
        e = vap_pres(hum[0], tem_dat['temperature'].data)
        rpg_dat['equivalent_potential_temperature'] = Theta+(spec_heat(tem_dat['temperature'].data)*.622*e/(p_baro-e)/1004.)*Theta/tem_dat['temperature'].data
        
        
    elif data_type == '2S02':
        
        c_list = get_coeff_list(site, 'tbx')
        _, _, _, coeff, _, _ = get_mvr_coeff(site, 'tbx', lev1['frequency'][:])
        index = np.where((lev1['ele'][:] > coeff['ele'] - .6) & (lev1['ele'][:] < coeff['ele'] + .6) & (lev1['pointing_flag'][:] == 0))[0]
        _, tb_ret = spectral_consistency(lev1, c_list)
        rpg_dat['tb_spectrum'] = tb_ret[index, :]
        rpg_dat['tb'] = lev1['tb'][index, :]
        fields = ['frequency', 'receiver_nb', 'receiver']
        for name in fields:
            rpg_dat[name] = lev1[name][:]
        flag = np.array(np.zeros(len(index)), dtype = bool)
        
    
    return rpg_dat, coeff, index, flag


def _combine_lev1(lev1: dict,
                  rpg_dat: dict,
                  index: np.ndarray) -> None:
    "add level1 data"    
    lev1_vars = ['time', 'time_bnds', 'station_latitude', 'station_longitude', 'azi', 'ele'] 
    for ivars in lev1_vars:
        if lev1[ivars].ndim > 1:
            rpg_dat[ivars] = lev1[ivars][index, :] 
        else:
            rpg_dat[ivars] = lev1[ivars][index]


def _mask_flag(rpg_dat: dict,
               flag: np.ndarray) -> None:
    "mask flagged data"        
    for ivars in rpg_dat:
        if (rpg_dat[ivars].ndim > 1) & (len(rpg_dat[ivars]) == len(rpg_dat['time'])):
            rpg_dat[ivars][flag, :] = ma.masked
        elif (rpg_dat[ivars].ndim == 1) & (len(rpg_dat[ivars]) == len(rpg_dat['time'])) & (ivars != 'time'):
            rpg_dat[ivars][flag] = ma.masked 
            
            
def _add_att(global_attributes: dict, 
             coeff: dict) -> None:
    "add retrieval attributes"
    fields = ['retrieval_type', 'retrieval_elevation_angles', 'retrieval_frequencies', 'retrieval_auxiliary_input', 'retrieval_description']
    for name in fields:
        global_attributes[name] = coeff[name]            


"specific heat for evaporation (J/kg)"
def spec_heat(T): return (2500.-2.42*(T-T0))*1000. 
"water vapor pressure"
def vap_pres(qs, T): return qs*Rw*T 


def rh_err(T: np.ndarray, 
           qs: np.ndarray,
           dT: np.ndarray,
           dq: np.ndarray) -> np.ndarray:
    "Calculates relative humidity error from absolute humidity and temperature"
    
    L = spec_heat(T)
    e = vap_pres(qs, T)
    es = mpcalc.saturation_vapor_pressure(masked_array(T, data_units='K')).magnitude
    # es = e0*np.exp((L/(Rw*T0))*((T-T0)/T))    
    
    "error propagation"
    drh_dq = Rw*T/es
    # des_dT = es*(T*(T-T0)*(-2420)+T0*L)/(Rw*T0*T**2)   
    des_dT = es * 17.67*243.5 / ((T-T0) + 243.5)**2
    drh_dT = qs*Rw/es**2*(es-T*des_dT)
    drh = np.sqrt((drh_dq * dq)**2 + (drh_dT * dT)**2)
    
    return drh
    

def calc_p_baro(T: np.ndarray,
           qs: np.ndarray,
           p: np.ndarray,
           z: np.ndarray) -> np.ndarray:
    "Calculate pressure in each level using barometric height formula"
    
    Tv = T*(1.+0.608*qs)
    p_baro = ma.masked_all(T.shape) 
    p_baro[(~qs.mask.any(axis=1)) & (~T.mask.any(axis=1)), 0] = p[(~qs.mask.any(axis=1)) & (~T.mask.any(axis=1))] * 100.
    for ialt in (np.arange(len(z) - 1) + 1): 
        p_baro[:, ialt] = p_baro[:, ialt-1] * np.exp(-9.81*(z[ialt] - z[ialt-1]) / (287. * (Tv[:, ialt] + Tv[:, ialt-1]) / 2.))
        
    return p_baro
