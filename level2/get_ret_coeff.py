import numpy as np
from numpy import ma
from utils import get_coeff_list
from scipy.interpolate import interp1d
import netCDF4 as nc

Fill_Value_Float = -999.
Fill_Value_Int = -99  


def get_mvr_coeff(site: str,
                  prefix: str,
                  freq: np.ndarray):
    """This function extracts retrieval coefficients for given files 
    and performs elevation angle interpolation for integrated quantities
    
    Args:
        site: Name of site.
        prefix: Identifier for type of product.
        freq: Frequencies of observations.
        
    Examples:
        >>> from level2.get_ret_coeff import get_mvr_coeff
        >>> get_mvr_coeff('site_name', 'lwp', [22, 31.4])
    """    

    c_list = get_coeff_list(site, prefix)
    cf = dict()    

    if prefix in ('lwp', 'iwv'):    
        
        cf['ele'] = np.ones(len(c_list)) * Fill_Value_Float    
        cf['freq'] = np.ones([len(freq), len(c_list)]) * Fill_Value_Float
        cf['coeff_lin'] = np.zeros([len(c_list), len(freq)])
        cf['coeff_quad'] = np.zeros([len(c_list), len(freq)])
        cf['offset'] = np.zeros(len(c_list))
        cf['err_ran'] = np.ones(len(c_list)) * Fill_Value_Float
        cf['err_sys'] = np.ones(len(c_list)) * Fill_Value_Float
        
        for i, file in enumerate(c_list):
            coeff = nc.Dataset(file)
            cf['ele'][i] = coeff['elevation_predictor'][i]
            _, freq_ind, freq_cf = np.intersect1d(freq[:], coeff['freq'][:], assume_unique = False, return_indices = True)
            if len(freq_ind) == 0:
                raise RuntimeError(['Instrument and retrieval frequencies do not match.']) 
                
            cf['freq'][freq_ind, i] = coeff['freq'][freq_cf]
            cf['coeff_lin'][i, freq_ind] = coeff['coefficient_mvr'][freq_cf]
            if coeff.regression_type == 'quadratic':
                cf['coeff_quad'][i, freq_ind] = coeff['coefficient_mvr'][freq_cf + len(freq_cf)]
            cf['offset'][i] = coeff['offset_mvr'][0]
            cf['err_ran'][i] = coeff['predictand_err'][0]  
            cf['err_sys'][i] = coeff['predictand_err_sys'][0]         
            
        def f_offset(x): return np.array([cf['offset'][(np.abs(cf['ele']-v)).argmin()] for v in x])
        def f_lin(x): return np.array([cf['coeff_lin'][(np.abs(cf['ele']-v)).argmin(), :] for v in x])
        def f_quad(x): return np.array([cf['coeff_quad'][(np.abs(cf['ele']-v)).argmin(), :] for v in x])
        def e_ran(x): return np.array([cf['err_ran'][(np.abs(cf['ele']-v)).argmin()] for v in x])
        def e_sys(x): return np.array([cf['err_sys'][(np.abs(cf['ele']-v)).argmin()] for v in x])         
            
    elif prefix in ('tze', 'hze'):
        
        coeff = nc.Dataset(c_list[0])
        n_height_grid = coeff.dimensions['n_height_grid'].size
        n_angles = len(c_list)

        cf['ele'] = np.ones(n_angles) * Fill_Value_Float    
        cf['freq'] = np.ones([len(freq), n_angles]) * Fill_Value_Float
        cf['coeff_lin'] = np.zeros([n_angles, n_height_grid, len(freq)])
        cf['coeff_quad'] = np.zeros([n_angles, n_height_grid, len(freq)])
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
            cf['coeff_lin'][i, :, freq_ind] = coeff['coefficient_mvr'][freq_cf, :]
            if coeff.regression_type == 'quadratic':
                cf['coeff_quad'][i, :, freq_ind] = coeff['coefficient_mvr'][freq_cf + len(freq_cf), :]
            cf['offset'][:, i] = coeff['offset_mvr'][:]
            cf['err_ran'][:, i] = coeff['predictand_err'][:]  
            cf['err_sys'][:, i] = coeff['predictand_err_sys'][:]   
            
        def f_offset(x): return np.array([cf['offset'][:, (np.abs(cf['ele']-v)).argmin()] for v in x])
        def f_lin(x): return np.array([cf['coeff_lin'][(np.abs(cf['ele']-v)).argmin(), :, :] for v in x])
        def f_quad(x): return np.array([cf['coeff_quad'][(np.abs(cf['ele']-v)).argmin(), :, :] for v in x])
        def e_ran(x): return np.array([cf['err_ran'][:, (np.abs(cf['ele']-v)).argmin()] for v in x])
        def e_sys(x): return np.array([cf['err_sys'][:, (np.abs(cf['ele']-v)).argmin()] for v in x])
        
        
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
        f_lin, f_quad = coeff['coefficient_mvr'][:, :], []
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
    if prefix == 'tel':
        cf['retrieval_frequencies'] = str(coeff['freq_bl'][:])
    else:
        cf['retrieval_frequencies'] = str(coeff['freq'][:])
    cf['retrieval_auxiliary_input'] = coeff.surface_mode
    cf['retrieval_description'] = coeff.retrieval_version    
    
    return f_offset, f_lin, f_quad, cf, e_ran, e_sys