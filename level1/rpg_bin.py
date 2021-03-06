""" This module contains all functions to read in RPG MWR binary files """
import numpy as np
import datetime
import utils
import logging

def get_rpg_bin(file_list: list) -> np.ndarray:    
    """ This function reads one day of a RPG MWR binary file type and concatenates the data. 
    Args:
        file_list: List of files for one day of a RPG MWR binary file type.
        
    Returns:
        Data array
    
    Example:
        >>> from level1.rpg_mwr import get_rpg_bin
        >>> get_rpg_bin('file_list') 
       
    """
    
    rpg_bin = RpgBin(file_list)
    return rpg_bin
    

def stack_files(file_list):    
    """ This function calls extension specific reader and stacks data and header. """
    
    def _stack_data(source, target, fun):
        for name, value in source.items():
            target[name] = (fun((target[name], value))
                            if name in target else value)  
                
    def _stack_header(source, target, fun):
        for name, value in source.items():
            if not name.startswith('_'):
                target[name] = (fun(target[name], value)
                                if name in target else value)  
            else:
                target[name] = value                                
                
    
    reader_name = str('read_' + file_list[0][-3:])
    data = {}
    header = {}
    
    for file in file_list:
        try:
            header_tmp, data_tmp = eval(reader_name + '(file)')
        except (TypeError, ValueError) as err:
            logging.warning(err)
            continue            
        _stack_header(header_tmp, header, np.add)        
        _stack_data(data_tmp,data, np.concatenate)

    return header, data


class RpgBin:
    def __init__(self, file_list):
        self.header, self.raw_data = stack_files(file_list)
        self.raw_data['time'] = utils.epoch2unix(self.raw_data['time'], self.header['_time_ref'])
        self.date = self._get_date()
        self.data = {}
        self._init_data()
        
    def _init_data(self):
        for key in self.raw_data:
            self.data[key] = self.raw_data[key]
        
    def _get_date(self):
        epoch = datetime.datetime(1970, 1, 1).timestamp()
        time_median = float(np.ma.median(self.raw_data['time']))
        time_median += epoch
        return datetime.datetime.utcfromtimestamp(time_median).strftime('%Y %m %d').split()    


def read_brt(file_name: str) -> dict:    
    """ This function reads RPG MWR .BRT binary files. """

    angle_calc = 0 
    """ Default is 0, 1 means calculate angles in old manner (for BRT files with file codes 666666) """

    with open(file_name, "rb") as file :
        
        code = np.fromfile( file, np.int32, 1)
        if code in (666000 , 666666):
            
            def _get_header():                
                """ Read header info """
                           
                n = int(np.fromfile( file, np.uint32, 1))
                time_ref = np.fromfile( file, np.uint32, 1)
                n_f = int(np.fromfile( file, np.int32, 1))
                f = np.fromfile( file, np.float32, n_f)
                xmin = np.fromfile( file, np.float32, n_f)
                xmax = np.fromfile( file, np.float32, n_f)
                
                header_names = ['_code', 'n', '_time_ref', '_n_f', '_f', '_xmin', '_xmax']                 
                header_values = [code, n, time_ref, n_f, f, xmin, xmax]
                header = dict(zip(header_names, header_values))
                return header
            
            def _create_variables():                
                """ Initialize data arrays """
                
                Fill_Value_Float = -999.
                Fill_Value_Int = -99
                vrs = {'time' : np.ones( header['n'], np.int32)*Fill_Value_Int,
                       'rain' : np.ones( header['n'], np.byte)*Fill_Value_Int,
                       'tb' : np.ones( [header['n'], header['_n_f']], np.float32)*Fill_Value_Float,
                       'ele' : np.ones( header['n'], np.float32)*Fill_Value_Float,
                       'azi' : np.ones( header['n'], np.float32)*Fill_Value_Float}
                return vrs
            
            def _angle_calc(ang, code):                
                """ Convert angle """
                
                if code == 666666:
                    if angle_calc == 1:
                        """ Angle calculation until 24.10.2014 """
                        sign = 1
                        if ang < 0:
                            sign = -1
                        az = sign*((ang/100.).astype(np.int32))/10.
                        el = ang - (sign*az*1000.)

                    elif angle_calc == 0:
                        els = ang-100.*((ang/100.).astype(np.int32))
                        azs = (ang-els)/1000.
                        if azs <= 360.:
                            el = els
                            az = azs
                        elif azs > 1000.:
                            az = azs - 1000.
                            el = 100. + els
                        else:
                            raise RuntimeError(['Error: Inconsistency in angle calculation!'])

                elif code == 666000:
                    a_str = str(ang[0])
                    el = float(a_str[0:-5])/100.
                    az = float(a_str[-5:])/100.
                return el, az
            
            def _get_data():
                """ Loop over file to read data """
                
                data = _create_variables()
                for sample in range(header['n']):  

                    data['time'][sample] = np.fromfile( file, np.int32, 1)
                    data['rain'][sample] = np.fromfile( file, np.byte, 1)
                    data['tb'][sample,] = np.fromfile( file, np.float32, header['_n_f'])
                    if code == 666666:
                        ang = np.fromfile( file, np.float32, 1)
                    elif code == 666000:
                        ang = np.fromfile( file, np.int32, 1)
                    data['ele'][sample], data['azi'][sample] = _angle_calc(ang,code)
                file.close()
                return data
            
            header = _get_header()                    
            data = _get_data()   
            return header, data
                        
        else:
            raise RuntimeError(['Error: BRT file code ' + str(code) + ' not suported'])


def read_met(file_name: str) -> dict:    
    """ This function reads RPG MWR .MET binary files. """

    with open(file_name, "rb") as file:
        
        code = np.fromfile( file, np.int32, 1)
        if code in (599658943 , 599658944):
            
            def _get_header():
                """ Read header info """
                
                Fill_Value_Float = -999.
                Fill_Value_Int = -99
                n = int(np.fromfile( file, np.uint32, 1))
                n_add = np.ones( 1, np.byte)
                if code == 599658944:
                    n_add = np.fromfile( file, np.byte, 1)
                n_sen = bin(int(n_add))
                xmin = np.ones( 3+n_sen.count('1'), np.float32)*Fill_Value_Float
                xmax = np.ones( 3+n_sen.count('1'), np.float32)*Fill_Value_Float
                for index in range(3+n_sen.count('1')):
                    xmin[index] = np.fromfile( file, np.float32, 1)
                    xmax[index] = np.fromfile( file, np.float32, 1)
                time_ref = np.fromfile( file, np.uint32, 1)
                
                header_names = ['_code','n', '_n_add', '_n_sen', '_xmin', '_xmax', '_time_ref']                
                header_values = [code, n, n_add, n_sen, xmin, xmax, time_ref]
                header = dict(zip(header_names, header_values))
                return header                
            
            def _create_variables():                
                """ Initialize data arrays """
                
                Fill_Value_Float = -999.
                Fill_Value_Int = -99
                vrs = {'time' : np.ones( header['n'], np.int32)*Fill_Value_Int,
                       'rain' : np.ones( header['n'], np.byte)*Fill_Value_Int,
                       'air_pressure' : np.ones( header['n'], np.float32)*Fill_Value_Float,
                       'air_temperature' : np.ones( header['n'], np.float32)*Fill_Value_Float,            
                       'relative_humidity' : np.ones( header['n'], np.float32)*Fill_Value_Float,
                       'adds' : np.ones( [header['n'], 3], np.float32)*Fill_Value_Float}
                       #'adds' : np.ones( [header['n'], header['_n_sen'].count('1')], np.float32)*Fill_Value_Float}
                return vrs            

            def _get_data():
                """ Loop over file to read data """
                
                data = _create_variables()
                for sample in range(header['n']): 
                    data['time'][sample] = np.fromfile( file, np.int32, 1)
                    data['rain'][sample] = np.fromfile( file, np.byte, 1)
                    data['air_pressure'][sample] = np.fromfile( file, np.float32,1)
                    data['air_temperature'][sample] = np.fromfile( file, np.float32, 1)
                    data['relative_humidity'][sample] = np.fromfile( file, np.float32, 1) / 100.
                    for add in range(header['_n_sen'].count('1')):
                        data['adds'][sample,add] = np.fromfile( file, np.float32, 1)
                file.close()
                return data
            
            header = _get_header()                    
            data = _get_data()   
            return header, data            
            
        else:
            raise RuntimeError(['Error: MET file code ' + str(code) + ' not suported'])
        
            
def read_irt(file_name: str) -> dict:    
    """ This function reads RPG MWR .IRT binary files. """

    with open(file_name, "rb") as file:
        
        code = np.fromfile( file, np.int32, 1)
        if code in (671112495,671112496,671112000):
            
            def _get_header():
                """ Read header info """
                
                n = int(np.fromfile( file, np.uint32, 1))
                xmin = np.fromfile( file, np.float32, 1)
                xmax = np.fromfile( file, np.float32, 1)
                time_ref = np.fromfile( file, np.uint32, 1)
                if code == 671112495:
                    n_f = 1
                    f = 11.1
                else:
                    n_f = int(np.fromfile( file, np.uint32, 1))
                    f = np.fromfile( file, np.float32, n_f )
                    
                header_names = ['_code','n','_xmin','_xmax','_time_ref','_n_f','_f']                    
                header_values = [code, n, xmin, xmax, time_ref, n_f, f]
                header = dict(zip(header_names, header_values))
                return header 
            
            def _create_variables():                
                """ Initialize data arrays """
                
                Fill_Value_Float = -999.
                Fill_Value_Int = -99
                vrs = {'time' : np.ones( header['n'], np.int32)*Fill_Value_Int,
                       'rain' : np.ones( header['n'], np.byte)*Fill_Value_Int,
                       'irt' : np.ones( [header['n'], header['_n_f']], np.float32)*Fill_Value_Float,
                       'ir_ele' : np.ones( header['n'], np.float32)*Fill_Value_Float,
                       'ir_azi' : np.ones( header['n'], np.float32)*Fill_Value_Float}
                return vrs
            
            def _angle_calc(ang, code):                
                """ Convert angle """
                
                if code == 671112496:
                    els = ang-100.*((ang/100.).astype(np.int32))
                    azs = (ang-els)/1000.
                    if azs <= 360.:
                        el = els
                        az = azs
                    elif azs > 1000.:
                        az = azs-1000.
                        el = 100.+els                    
                elif code == 671112000:
                    a_str = str(ang[0])
                    el = float(a_str[0:-5])/100.
                    az = float(a_str[-5:])/100.
                return el, az
                    

            def _get_data():
                """ Loop over file to read data """
                
                data = _create_variables()
                for sample in range(header['n']): 
                    data['time'][sample] = np.fromfile( file, np.int32, 1)
                    data['rain'][sample] = np.fromfile( file, np.byte, 1)
                    data['irt'][sample,] = np.fromfile( file, np.float32, header['_n_f']) + 273.15
                    if code == 671112496:
                        ang = np.fromfile( file, np.float32, 1)
                    elif code == 671112000:
                        ang = np.fromfile( file, np.int32, 1)
                    data['ir_ele'][sample], data['ir_azi'][sample] = _angle_calc(ang,code)                    
                file.close()
                return data
            
            header = _get_header()                    
            data = _get_data()   
            return header, data 
        
        else:
            raise RuntimeError(['Error: IRT file code ' + str(code) + ' not suported'])
     
        
def read_blb(file_name: str) -> dict:    
    """ This function reads RPG MWR .BLB binary files. """

    with open(file_name, "rb") as file:
        
        code = np.fromfile( file, np.int32, 1)
        if code in (567845847 , 567845848):
            
            def _get_header():
                """ Read header info """
                
                n = int(np.fromfile( file, np.uint32, 1))
                if code == 567845848:
                    n_f = int(np.fromfile( file, np.int32, 1))
                else:
                    n_f = 14
                xmin = np.fromfile( file, np.float32, n_f)
                xmax = np.fromfile( file, np.float32, n_f)
                time_ref = np.fromfile( file, np.uint32, 1)
                if code == 567845847:
                    n_f = int(np.fromfile( file, np.int32, 1))
                f = np.fromfile( file, np.float32, n_f)
                n_ang = int(np.fromfile( file, np.int32, 1)) + 1
                ang = np.append(np.fromfile( file, np.float32, n_ang-1), 0)
                
                header_names = ['_code','n','_xmin','_xmax','_time_ref','_n_f','_f','_n_ang','_ang']
                header_values = [code, n, xmin, xmax, time_ref, n_f, f, n_ang, ang]
                header = dict(zip(header_names, header_values))
                return header                 

            def _create_variables():                
                """ Initialize data arrays """
                
                Fill_Value_Float = -999.
                Fill_Value_Int = -99
                vrs = {'time' : np.ones( header['n'], np.int32)*Fill_Value_Int,
                       'rf_mod' : np.ones( header['n'], np.byte)*Fill_Value_Int,
                       'tb' : np.ones( [header['n'], header['_n_f'], header['_n_ang']], np.float32)*Fill_Value_Float}
                return vrs

            def _get_data():
                """ Loop over file to read data """
                
                data = _create_variables()
                for sample in range(header['n']): 
                    data['time'][sample] = np.fromfile( file, np.int32, 1)
                    data['rf_mod'][sample] = np.fromfile( file, np.byte, 1)
                    for freq in range(header['_n_f']):
                        data['tb'][sample, freq, ] = np.fromfile( file, np.float32, header['_n_ang'])   
                file.close()
                return data
            
            header = _get_header()                    
            data = _get_data()   
            return header, data 
        
        else:
            raise RuntimeError(['Error: BLB file code ' + str(code) + ' not suported'])
        
        
def read_hkd(file_name: str) -> dict:    
    """ This function reads RPG MWR .HKD binary files. """

    with open(file_name, "rb") as file:
        
        code = np.fromfile( file, np.int32, 1)
        if code == 837854832:
            
            def _get_header():
                """ Read header info """
                
                n = int(np.fromfile( file, np.uint32, 1))
                time_ref = np.fromfile( file, np.uint32, 1)
                sel = np.fromfile( file, np.uint32, 1)
                
                header_names = ['_code','n','_time_ref','_sel']
                header_values = [code, n, time_ref, sel]
                header = dict(zip(header_names,header_values))
                return header                   

            def _create_variables():                
                """ Initialize data arrays """
                
                Fill_Value_Float = -999.
                Fill_Value_Int = -99
                vrs = {'time' : np.ones( header['n'], np.int32)*Fill_Value_Int,
                       'alarm' : np.ones( header['n'], np.byte)*Fill_Value_Int, 
                       'station_longitude' : np.ones( header['n'], np.float32)*Fill_Value_Float,
                       'station_latitude' : np.ones( header['n'], np.float32)*Fill_Value_Float,          
                       'temp' : np.ones( [header['n'],4], np.float32)*Fill_Value_Float,
                       'stab' : np.ones( [header['n'],2], np.float32)*Fill_Value_Float,           
                       'flash' : np.ones( header['n'], np.int32)*Fill_Value_Int,      
                       'qual' : np.ones( header['n'], np.int32)*Fill_Value_Int,          
                       'status' : np.ones( header['n'], np.int32)*Fill_Value_Int}
                return vrs

            def _get_data():
                """ Loop over file to read data """
                
                data = _create_variables()
                for sample in range(header['n']): 
                    data['time'][sample] = np.fromfile( file, np.int32, 1)
                    data['alarm'][sample] = np.fromfile( file, np.byte, count=1)
                    if (header['_sel'] & 1):
                        data['station_longitude'][sample] = np.fromfile( file, np.float32, count=1)
                        data['station_latitude'][sample] = np.fromfile( file, np.float32, count=1)
                    if (header['_sel'] & 2):
                        data['temp'][sample,] = np.fromfile( file, np.float32, count=4)
                    if (header['_sel'] & 4):
                        data['stab'][sample,] = np.fromfile( file, np.float32, count=2)
                    if (header['_sel'] & 8):
                        data['flash'][sample] = np.fromfile( file, np.int32, count=1)
                    if (header['_sel'] & 16):
                        data['qual'][sample] = np.fromfile( file, np.int32, count=1)
                    if (header['_sel'] & 32):
                        data['status'][sample] = np.fromfile( file, np.int32, count=1)
                file.close()
                return data
            
            header = _get_header()                    
            data = _get_data()   
            return header, data 
        
        else:
            raise RuntimeError(['Error: HKD file code ' + str(code) + ' not suported'])

            
def read_spc(file_name: str) -> dict:    
    """ This function reads RPG MWR .SPC binary files. """

    angle_calc = 0 
    """ Default is 0, 1 means calculate angles in old manner (for SPC files with file codes 666667) """

    with open(file_name, "rb") as file :
        
        code = np.fromfile( file, np.int32, 1)
        if code in (666667 , 667000):
            
            def _get_header():                
                """ Read header info """
                           
                n = int(np.fromfile( file, np.uint32, 1))
                time_ref = np.fromfile( file, np.uint32, 1)
                n_f = int(np.fromfile( file, np.int32, 1))
                f = np.fromfile( file, np.float32, n_f)
                xmin = np.fromfile( file, np.float32, n_f)
                xmax = np.fromfile( file, np.float32, n_f)
                
                header_names = ['_code', 'n', '_time_ref', '_n_f', '_f', '_xmin', '_xmax']                 
                header_values = [code, n, time_ref, n_f, f, xmin, xmax]
                header = dict(zip(header_names, header_values))
                return header
            
            def _create_variables():                
                """ Initialize data arrays """
                
                Fill_Value_Float = -999.
                Fill_Value_Int = -99
                vrs = {'time' : np.ones( header['n'], np.int32)*Fill_Value_Int,
                       'rain' : np.ones( header['n'], np.byte)*Fill_Value_Int,
                       'tb' : np.ones( [header['n'], header['_n_f']], np.float32)*Fill_Value_Float,
                       'ele' : np.ones( header['n'], np.float32)*Fill_Value_Float,
                       'azi' : np.ones( header['n'], np.float32)*Fill_Value_Float}
                return vrs
            
            def _angle_calc(ang, code):                
                """ Convert angle """
                
                if code == 666667:
                    if angle_calc == 1:
                        """ Angle calculation until 24.10.2014 """
                        sign = 1
                        if ang < 0:
                            sign = -1
                        az = sign*((ang/100.).astype(np.int32))/10.
                        el = ang - (sign*az*1000.)

                    elif angle_calc == 0:
                        els = ang-100.*((ang/100.).astype(np.int32))
                        azs = (ang-els)/1000.
                        if azs <= 360.:
                            el = els
                            az = azs
                        elif azs > 1000.:
                            az = azs - 1000.
                            el = 100. + els
                        else:
                            raise RuntimeError(['Error: Inconsistency in angle calculation!'])

                elif code == 667000:
                    a_str = str(ang[0])
                    el = float(a_str[0:-5])/100.
                    az = float(a_str[-5:])/100.
                return el, az
            
            def _get_data():
                """ Loop over file to read data """
                
                data = _create_variables()
                for sample in range(header['n']):  

                    data['time'][sample] = np.fromfile( file, np.int32, 1)
                    data['rain'][sample] = np.fromfile( file, np.byte, 1)
                    data['tb'][sample,] = np.fromfile( file, np.float32, header['_n_f'])
                    if code == 666667:
                        ang = np.fromfile( file, np.float32, 1)
                    elif code == 667000:
                        ang = np.fromfile( file, np.int32, 1)
                    data['ele'][sample], data['azi'][sample] = _angle_calc(ang,code)
                file.close()
                return data
            
            header = _get_header()                    
            data = _get_data()   
            return header, data
                        
        else:
            raise RuntimeError(['Error: SPC file code ' + str(code) + ' not suported'])