from read_specs import get_site_specs
from level2.lev2_meta_nc import get_data_attributes
from level1.quality_control import spectral_consistency
from utils import df_interp, get_coeff_list, rebin_2d
import rpg_mwr
import numpy as np
from numpy import ma
import netCDF4 as nc
from scipy.interpolate import interp1d
import pandas as pd
from pandas.tseries.frequencies import to_offset
import metpy.calc as mpcalc
from metpy.units import units

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
    rpg_mwr.save_rpg(hatpro, output_file, global_attributes, data_type, params, site)
    
    
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
        rpg_dat['lwp'] = coeff_offset + np.sum(coeff_lin.T * lev1['tb'][index, :], axis = 1) + np.sum(coeff_quad.T * lev1['tb'][index, :] **2, axis = 1)
        rpg_dat['lwp_random_error'] = ran_err(lev1['ele'][index])
        rpg_dat['lwp_systematic_error'] = sys_err(lev1['ele'][index])         

        freq_31 = np.where(lev1['frequency'][:] == 31.4)[0]
        if len(freq_31) == 1:
            lwp_offset = get_lwp_offset(lev1['time'][index], np.copy(rpg_dat['lwp']), np.squeeze(lev1['tb'][index, freq_31]), lev1, index, site)
            rpg_dat['lwp'] = rpg_dat['lwp'] - lwp_offset
            rpg_dat['lwp_offset'] = lwp_offset
        else:
            rpg_dat['lwp_offset'] = np.ones(len(index)) * Fill_Value_Float          
        
        
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
            rpg_dat['temperature'][:, ialt] = coeff_offset[ialt, :].T + np.sum(coeff_lin[ialt, :, :].T * lev1['tb'][index, :], axis = 1) + np.sum(coeff_quad[ialt, :, :].T * lev1['tb'][index, :] **2, axis = 1)    

        rpg_dat['temperature_random_error'] = ran_err(lev1['ele'][index]).T
        rpg_dat['temperature_systematic_error'] = sys_err(lev1['ele'][index]).T                     
                
                
    elif data_type == '2P02':
                
        offset, lin, quad, coeff, ran_err, sys_err = get_mvr_coeff(site, 'tel', lev1['frequency'][:])

        _, freq_ind, _ = np.intersect1d(lev1['frequency'][:], coeff['freq'], assume_unique = False, return_indices = True)
        _, freq_bl, _ = np.intersect1d(coeff['freq'], coeff['freq_bl'], assume_unique = False, return_indices = True)
                    
        ix0 = np.where((lev1['ele'][:] > coeff['ele'][0] - .6) & (lev1['ele'][:] < coeff['ele'][0] + .6) & (lev1['pointing_flag'][:] == 1))[0]

        if len(ix0) == 0:
            raise RuntimeError(['No suitable data found for processing for data type: ' + data_type])
            
        ele = []
        flag = []
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
        
        if len(index) == 0:
            raise RuntimeError(['No suitable data found for processing for data type: ' + data_type])
            
        rpg_dat['altitude'] = coeff['height_grid']            
        rpg_dat['water_vapor_vmr'] = ma.masked_all((len(index), coeff['n_height_grid']))  

        coeff_offset = offset(lev1['ele'][index])
        coeff_lin = lin(lev1['ele'][index])
        coeff_quad = quad(lev1['ele'][index])               

        for ialt, _ in enumerate(rpg_dat['altitude']):
            rpg_dat['water_vapor_vmr'][:, ialt] = coeff_offset[ialt, :].T + np.sum(coeff_lin[ialt, :, :].T * lev1['tb'][index, :], axis = 1) + np.sum(coeff_quad[ialt, :, :].T * lev1['tb'][index, :] **2, axis = 1)    

        rpg_dat['water_vapor_vmr_random_error'] = ran_err(lev1['ele'][index]).T
        rpg_dat['water_vapor_vmr_systematic_error'] = sys_err(lev1['ele'][index]).T              
             
                
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
        rpg_dat['relative_humidity'] = (hum[0]*tem_dat['temperature'][:,:]*462.) / mpcalc.saturation_vapor_pressure(tem_dat['temperature'].data*units.K)
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
        Tv = tem_dat['temperature'][:,:]*(1.+0.608*hum[0])
        p_baro = ma.masked_all((len(index), len(rpg_dat['altitude']))) 
        p_baro[:, 0] = lev1['air_pressure'][index] * 100.
        for ialt in (np.arange(len(rpg_dat['altitude']) - 1) + 1): 
            p_baro[:, ialt] = p_baro[:, ialt-1] * np.exp(-9.81*(tem_dat['altitude'][ialt] - tem_dat['altitude'][ialt-1]) / (287. * (Tv[:, ialt] + Tv[:, ialt-1]) / 2.))
        pdb.set_trace()
        rpg_dat['potential_temperature'] = np.array(mpcalc.potential_temperature(p_baro.data*units.Pa, tem_dat['temperature'].data*units.K))        
  
            
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
        Tv = tem_dat['temperature'][:,:]*(1.+0.608*hum[0])
        p_baro = ma.masked_all((len(index), len(rpg_dat['altitude']))) 
        p_baro[:, 0] = lev1['air_pressure'][index] * 100.
        for ialt in (np.arange(len(rpg_dat['altitude']) - 1) + 1): 
            p_baro[:, ialt] = p_baro[:, ialt-1] * np.exp(-9.81*(tem_dat['altitude'][ialt] - tem_dat['altitude'][ialt-1]) / (287. * (Tv[:, ialt] + Tv[:, ialt-1]) / 2.))
        dp = mpcalc.dewpoint(hum[0].data*tem_dat['temperature'].data*462.*units.Pa)
        rpg_dat['equivalent_potential_temperature'] = np.array(mpcalc.equivalent_potential_temperature(p_baro.data*units.Pa, tem_dat['temperature'].data*units.K, dp))
        
        
    elif data_type == '2S02':
        
        c_list = get_coeff_list(site, 'tbx')
        _, _, _, coeff, _, _ = get_mvr_coeff(site, 'tbx', lev1['frequency'][:])
        index = np.where((lev1['ele'][:] > coeff['ele'] - .6) & (lev1['ele'][:] < coeff['ele'] + .6) & (lev1['pointing_flag'][:] == 0))[0]
        _, tb_ret = spectral_consistency(lev1, c_list, params['tbx_f'])
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
        elif (rpg_dat[ivars].ndim == 1) & (len(rpg_dat[ivars]) == len(rpg_dat['time'])):
            rpg_dat[ivars][flag] = ma.masked 
            
            
def _add_att(global_attributes: dict, 
             coeff: dict) -> None:
    "add retrieval attributes"
    fields = ['retrieval_type', 'retrieval_elevation_angles', 'retrieval_frequencies', 'retrieval_auxiliary_input', 'retrieval_description']
    for name in fields:
        global_attributes[name] = coeff[name]            


def get_mvr_coeff(site: str,
                  prefix: str,
                  freq: np.ndarray):
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
            if len(freq_ind) == 0:
                raise RuntimeError(['Instrument and retrieval frequencies do not match.']) 
                
            cf['freq'][freq_ind, i] = coeff['freq'][freq_cf]
            cf['coeff_lin'][freq_ind, i] = coeff['coefficient_mvr'][freq_cf]
            cf['coeff_quad'][freq_ind, i] = coeff['coefficient_mvr'][freq_cf + len(freq_cf)]
            cf['offset'][i] = coeff['offset_mvr'][0]
            cf['err_ran'][i] = coeff['predictand_err'][0]
            cf['err_sys'][i] = coeff['predictand_err_sys'][0]

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
        cf['err_ran'] = ma.masked_all((n_height_grid, n_angles)) 
        cf['err_sys'] = ma.masked_all((n_height_grid, n_angles))  
        cf['n_height_grid'] = n_height_grid
        cf['height_grid'] = coeff['height_grid']
        
        
        for i, file in enumerate(c_list):
            coeff = nc.Dataset(file)
            cf['ele'][i] = coeff['elevation_predictor'][i]
            _, freq_ind, freq_cf = np.intersect1d(freq[:], coeff['freq'][:], assume_unique = False, return_indices = True)
            if len(freq_ind) == 0:
                raise RuntimeError(['Instrument and retrieval frequencies do not match.']) 
                
            cf['freq'][freq_ind, i] = coeff['freq'][freq_cf]
            cf['coeff_lin'][:, freq_ind, i] = np.transpose(coeff['coefficient_mvr'][freq_cf, :])
            cf['coeff_quad'][:, freq_ind, i] = np.transpose(coeff['coefficient_mvr'][freq_cf + len(freq_cf), :])
            cf['offset'][:, i] = coeff['offset_mvr'][:]
            cf['err_ran'][:, i] = coeff['predictand_err'][:]  
            cf['err_sys'][:, i] = coeff['predictand_err_sys'][:]   

        f_offset = lambda x: np.transpose(np.array([cf['offset'][:, (np.abs(cf['ele']-v)).argmin()] for v in x]))
        f_lin = lambda x: np.transpose(np.array([cf['coeff_lin'][:, :, (np.abs(cf['ele']-v)).argmin()] for v in x]), (1, 2, 0))
        f_quad = lambda x: np.transpose(np.array([cf['coeff_quad'][:, :, (np.abs(cf['ele']-v)).argmin()] for v in x]), (1, 2, 0))
        e_ran = lambda x: np.transpose(np.array([cf['err_ran'][:, (np.abs(cf['ele']-v)).argmin()] for v in x]))
        e_sys = lambda x: np.transpose(np.array([cf['err_sys'][:, (np.abs(cf['ele']-v)).argmin()] for v in x]))  
        
        
    elif prefix == 'tel':
        
        coeff = nc.Dataset(c_list[0])
        _, freq_ind, freq_cf = np.intersect1d(freq[:], coeff['freq'][:], assume_unique = False, return_indices = True)
        if len(freq_ind) == 0:
            raise RuntimeError(['Instrument and retrieval frequencies do not match.']) 
            
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
            
            
    elif prefix == 'tbx':        
        
        coeff = nc.Dataset(c_list[0])
        cf['ele'] = coeff['elevation_predictor'][:]
        cf['freq'] = freq[:]
        e_ran =  coeff['predictand_err'][:] 
        e_sys = coeff['predictand_err_sys'][:]         
        f_offset, f_lin, f_quad, = [], [], []
        
        
    else:
        raise RuntimeError(['Prefix '+ prefix +' not recognized for retrieval coefficient file(s).'])
        
    cf['retrieval_type'] = coeff.regression_type
    cf['retrieval_elevation_angles'] = str(cf['ele'])
    cf['retrieval_frequencies'] = str(coeff['freq'][:])
    cf['retrieval_auxiliary_input'] = coeff.surface_mode
    cf['retrieval_description'] = coeff.retrieval_version    
    
    return f_offset, f_lin, f_quad, cf, e_ran, e_sys


def rh_err(T: np.ndarray, 
           q: np.ndarray,
           dT: np.ndarray,
           dq: np.ndarray):
    "Calculate relative humidity error from absolute humidity and temperature using analytical approximation of Clausius-Clapeyron equation (+ error propagation)"
    
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
    
    "error propagation"
    drh_dq = Rw*T/es
    des_dT = es*(T*(T-T0)*(-2420)+T0*L)/(Rw*T0*T**2)   
    drh_dT = q*Rw/es**2*(es-T*des_dT)
    drh = np.sqrt((drh_dq * dq)**2 + (drh_dT * dT)**2)
    
    return drh
    
    
def get_lwp_offset(time: np.ndarray,
                   lwp: np.ndarray,
                   tb: np.ndarray,
                   lev1: dict,
                   index: np.ndarray, 
                   site: str) -> np.ndarray:
    "Correct Lwp offset using 2min standard deviation and IRT"
    
    tb_df = pd.DataFrame({'Tb': tb}, index = pd.to_datetime(time, unit = 's'))
    tb_cc = tb_df.resample("2min", origin = 'start', closed = 'left', label = 'left').count()
    tb_cc.index = tb_cc.index + to_offset('1min')
    tb_std = tb_df.resample("2min", origin = 'start', closed = 'left', label = 'left').std()
    tb_std.index = tb_std.index + to_offset('1min')
    tb_std[tb_cc['Tb'] < 30] = np.nan    
    tb_mx = tb_std.resample("20min", origin = 'start', closed = 'left', label = 'left').max()
    tb_mx.index = tb_mx.index + to_offset('10min')
    tb_mx = df_interp(tb_mx, tb_df.index)      
    
    if 'irt' in lev1.variables:
        irt = lev1['irt'][index, 0]
        irt_df = pd.DataFrame({'Irt': irt[:]}, index = pd.to_datetime(time, unit = 's'))
        irt_mx = irt_df.resample("20min", origin = 'start', closed = 'left', label = 'left').max()
        irt_mx.index = irt_mx.index + to_offset('10min')
        irt_mx = df_interp(irt_mx, irt_df.index) 
        lwp[(irt_mx['Irt'] > 233.15) & (tb_mx['Tb'] > .2) | np.isnan(irt_mx['Irt']) | np.isnan(tb_mx['Tb']) | (lwp > .1)] = np.nan
    else:
        lwp[(tb_mx['Tb'] > .2) | np.isnan(tb_mx['Tb']) | (lwp > .1)] = np.nan
        
    lwp_df = pd.DataFrame({'Lwp': lwp}, index = pd.to_datetime(time, unit = 's'))
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
