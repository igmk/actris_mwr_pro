from collections import namedtuple

def get_data_attributes(rpg_variables: dict, 
                        data_type: str) -> dict:
    """Adds Metadata for RPG MWR Level 2 variables for NetCDF file writing.
    Args:
        rpg_variables: RpgArray instances.
        data_type: Data type of the netCDF file.
        
    Returns:
        Dictionary
        
    Raises:
        RuntimeError: Specified data type is not supported.
    
    Example:
        from level2.lev2_meta_nc import get_data_attributes
        att = get_data_attributes('data','data_type')    
    """
    
    if data_type in ('2P00', '2I00', '2S00'):
        
        attributes = eval('ATTRIBUTES_'+ data_type)
        for key in list(rpg_variables):
            if key in attributes:
                rpg_variables[key].set_attributes(attributes[key])
            else:
                del rpg_variables[key]
                
        index_map = {v: i for i, v in enumerate(attributes)}
        rpg_variables = dict(sorted(rpg_variables.items(), key=lambda pair: index_map[pair[0]]))
        
        return rpg_variables
                
    else:
        raise RuntimeError(['Data type '+ data_type +' not supported for file writing.'])


FIELDS = (
    'long_name',
    'standard_name',
    'units',
    'definition',
    'comment')

MetaData = namedtuple('MetaData', FIELDS)
MetaData.__new__.__defaults__ = (None,) * len(MetaData._fields)


ATTRIBUTES_2P00 = {
    'time': MetaData(
        long_name='Time (UTC) of the measurement',
        standard_name='time',
        units='seconds since 1970-01-01 00:00:00.000',
        comment='Time indication of samples is at end of integration-time',
    ),
    'time_bounds': MetaData(
        long_name='Start and end time (UTC) of the measurements',
        units='seconds since 1970-01-01 00:00:00.000',
    ),
    'station_latitude': MetaData(
        long_name='Latitude of measurement station',
        standard_name='latitude',
        units='degrees_north',
    ),
    'station_longitude': MetaData(
        long_name='Longitude of measurement station',
        standard_name='longitude',
        units='degrees_east',
    ),
    'station_altitude': MetaData(
        long_name='Altitude above mean sea level of measurement station',
        standard_name='altitude',
        units='m',
    ),
    'azi': MetaData(
        long_name='Sensor azimuth angle',
        standard_name='sensor_azimuth_angle',
        units='degrees',
        comment='0=North, 90=East, 180=South, 270=West',
    ),
    'ele': MetaData(
        long_name='Sensor elevation angle',
        standard_name='sensor_elevation_angle',
        units='degrees',
        comment='0=horizon, 90=zenith',
    ),
    'altitude': MetaData(
        long_name='Altitude above sea level',
        units='m',
    ),
    'temperature': MetaData(
        long_name='Retrieved temperature profile',
        units='K',
    ),
    'temperature_random_error': MetaData(
        long_name='Random uncertainty of retrieved temperature profile',
        units='K',
        comment='specify here source of this variable',
    ),
    'temperature_systematic_error': MetaData(
        long_name='Systematic uncertainty of retrieved temperature profile',
        units='K',
        comment='specify here source of this variable',
    ),    
    'water_vapor_vmr': MetaData(
        long_name='Retrieved water vapour profile',
        units='ppm',
    ),  
    'water_vapor_vmr_random_error': MetaData(
        long_name='Random uncertainty of retrieved water',
        units='ppm',
        comment='specify here source of this variable',
    ),     
    'water_vapor_vmr_systematic_error': MetaData(
        long_name='Systematic uncertainty of retrieved water vapour profile',
        units='ppm',
        comment='specify here source of this variable',
    ),     
    'relative_humidity': MetaData(
        long_name='Relative Humidity',
        units='1',
    ),      
    'relative_humidity_random_error': MetaData(
        long_name='Random uncertainty of relative humidity',
        units='1',
        comment='specify here source of this variable',
    ),  
    'relative_humidity_systematic_error': MetaData(
        long_name='Systematic uncertainty of relative humidity',
        units='1',
        comment='specify here source of this variable',
    ),    
    'liquid_content': MetaData(
        long_name='Liquid water content',
        units='g/m^3'
    ),      
    'liquid_content_random_error': MetaData(
        long_name='Random uncertainty of liquid water content',
        units='g/m^3',
        comment='specify here source of this variable',
    ),  
    'liquid_content_systematic_error': MetaData(
        long_name='Systematic uncertainty of liquid water content',
        units='g/m^3',
        comment='specify here source of this variable',
    ),        
}


ATTRIBUTES_2I00 = {
    'time': MetaData(
        long_name='Time (UTC) of the measurement',
        standard_name='time',
        units='seconds since 1970-01-01 00:00:00.000',
        comment='Time indication of samples is at end of integration-time',
    ),
    'time_bounds': MetaData(
        long_name='Start and end time (UTC) of the measurements',
        units='seconds since 1970-01-01 00:00:00.000',
    ),
    'station_latitude': MetaData(
        long_name='Latitude of measurement station',
        standard_name='latitude',
        units='degrees_north',
    ),
    'station_longitude': MetaData(
        long_name='Longitude of measurement station',
        standard_name='longitude',
        units='degrees_east',
    ),
    'station_altitude': MetaData(
        long_name='Altitude above mean sea level of measurement station',
        standard_name='altitude',
        units='m',
    ),
    'azi': MetaData(
        long_name='Sensor azimuth angle',
        standard_name='sensor_azimuth_angle',
        units='degrees',
        comment='0=North, 90=East, 180=South, 270=West',
    ),
    'ele': MetaData(
        long_name='Sensor elevation angle',
        standard_name='sensor_elevation_angle',
        units='degrees',
        comment='0=horizon, 90=zenith',
    ),    
    'Iwv': MetaData(
        long_name='Retrieved column-integrated water vapour',
        units='mm',
    ),     
    'Iwv_random_error': MetaData(
        long_name='Random uncertainty of retrieved column-integrated water vapour',
        units='mm',
        comment='specify here source of this variable',
    ),  
    'Iwv_systematic_error': MetaData(
        long_name='Systematic uncertainty of retrieved column-integrated water vapour',
        units='mm',
        comment='specify here source of this variable',
    ),  
    'Lwp': MetaData(
        long_name='Retrieved column-integrated liquid water path',
        units='mm',
    ),     
    'Lwp_random_error': MetaData(
        long_name='Random uncertainty of retrieved column-integrated liquid water path',
        units='mm',
        comment='specify here source of this variable',
    ),  
    'Lwp_systematic_error': MetaData(
        long_name='Systematic uncertainty of retrieved column-integrated liquid water path',
        units='mm',
        comment='specify here source of this variable',
    ),     
}


ATTRIBUTES_2S00 = {
    'time': MetaData(
        long_name='Time (UTC) of the measurement',
        standard_name='time',
        units='seconds since 1970-01-01 00:00:00.000',
        comment='Time indication of samples is at end of integration-time',
    ),
    'time_bounds': MetaData(
        long_name='Start and end time (UTC) of the measurements',
        units='seconds since 1970-01-01 00:00:00.000',
    ),
    'station_latitude': MetaData(
        long_name='Latitude of measurement station',
        standard_name='latitude',
        units='degrees_north',
    ),
    'station_longitude': MetaData(
        long_name='Longitude of measurement station',
        standard_name='longitude',
        units='degrees_east',
    ),
    'station_altitude': MetaData(
        long_name='Altitude above mean sea level of measurement station',
        standard_name='altitude',
        units='m',
    ),
    'azi': MetaData(
        long_name='Sensor azimuth angle',
        standard_name='sensor_azimuth_angle',
        units='degrees',
        comment='0=North, 90=East, 180=South, 270=West',
    ),
    'ele': MetaData(
        long_name='Sensor elevation angle',
        standard_name='sensor_elevation_angle',
        units='degrees',
        comment='0=horizon, 90=zenith',
    ),    
}