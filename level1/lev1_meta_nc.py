from collections import namedtuple


def get_data_attributes(rpg_variables: dict, 
                        data_type: str) -> dict:
    """Adds Metadata for RPG MWR Level 1 variables for NetCDF file writing.
    Args:
        rpg_variables: RpgArray instances.
        data_type: Data type of the netCDF file.
        
    Returns:
        Dictionary
        
    Raises:
        RuntimeError: Specified data type is not supported.
    
    Example:
        from level1.lev1_meta_nc import get_data_attributes
        att = get_data_attributes('data','data_type')    
    """
    
    if data_type in ('1B01', '1B11', '1B21'):        
        attributes = dict(ATTRIBUTES_COM, **eval('ATTRIBUTES_'+ data_type))
        
    elif data_type == '1C01':
        attributes = dict(ATTRIBUTES_COM, **ATTRIBUTES_1B01, **ATTRIBUTES_1B11, **ATTRIBUTES_1B21)
        
    else:
        raise RuntimeError(['Data type '+ data_type +' not supported for file writing.'])
        
    for key in list(rpg_variables):
        if key in attributes:
            rpg_variables[key].set_attributes(attributes[key])
        else:
            del rpg_variables[key]

    index_map = {v: i for i, v in enumerate(attributes)}
    rpg_variables = dict(sorted(rpg_variables.items(), key=lambda pair: index_map[pair[0]]))

    return rpg_variables


FIELDS = (
    'long_name',
    'standard_name',
    'units',
    'definition',
    'comment')

MetaData = namedtuple('MetaData', FIELDS)
MetaData.__new__.__defaults__ = (None,) * len(MetaData._fields)


ATTRIBUTES_COM = {
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
}


DEFINITIONS_1B01 = {
    'quality_flag':
        ('\n'
         'Bit 1: Missing TB-value\n'
         'Bit 2: TB threshold (lower range)\n'
         'Bit 3: TB threshold (upper range)\n'
         'Bit 4: Spectral consistency threshold\n'
         'Bit 5: Receiver sanity\n'
         'Bit 6: Rain flag\n'
         'Bit 7: Solar flag\n'
         'Bit 8: TB offset threshold')
}

ATTRIBUTES_1B01 = {
    'frequency': MetaData(
        long_name='Frequency of microwave channels',
        standard_name='sensor_band_central_radiation_frequency',
        units='GHz',
    ),
    'bandwidth': MetaData(
        long_name='Bandwidth of the central frequency',
        standard_name='sensor_band_spectral_width',
        units='GHz',
        comment='center frequency of single of upper side-band',
    ),
    'sideband_count': MetaData(
        long_name='Single, double, or double-double sideband',
        comment='1, 2, 4 are possible values',
    ),
    'sideband_IF_separation': MetaData(
        long_name='56.xx +/- X +/- Y',
    ),
    'beamwidth': MetaData(
        long_name='Full width at half maximum',
        units='degrees',
    ),
    'freq_shift': MetaData(
        long_name='frequency shift applied to correct measured brightness temperature for frequency offset of microwave radiometer channel',
        units='GHz',
    ),
    'tb': MetaData(
        long_name='Brightness temperatures',
        standard_name='brightness_temperature',
        units='K',
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
    'tb_accuracy': MetaData(
        long_name='Total absolute calibration uncertainty of brightness temperature, one standard deviation',
        units='K',
        comment='specify here source of this variable, e.g. literature value, specified by manufacturer, result of validation effort (updated irregularily) For RDX systems, derived from analysis performed by Tim Hewsion (Tim J. Hewison, 2006: Profiling Temperature and Humidity by Ground-based Microwave Radiometers, PhD Thesis, University of Reading.) Derived from sensitivity analysis of LN2 calibration plus instrument noise levels (ACTRIS work), currently literature values (Maschwitz et al. for HATPRO, ? for radiometrics)',
    ),
    'tb_cov': MetaData(
        long_name='Error covariance matrix of brightness temperature channels',
        units='K*K',
        comment='specify here standardized method used to estimate the covariance matrix',
    ),
    'quality_flag': MetaData(
        long_name='Quality_flag',
        standard_name='quality_flag',
        units='1 (bit variable)',
        definition=DEFINITIONS_1B01['quality_flag'],
        comment='Bit 4 and 8 are calculated centrally at the E-PROFILE processing hub',
    ),
    'pointing_flag': MetaData(
        long_name='Pointing Flag',
        standard_name='pointing_flag',
        units='1 (bit variable)',
        comment='Flag indicating observation mode - single vs. multiple pointing',
    ),
    't_amb': MetaData(
        long_name='Ambient target temperature',
        units='K',
    ),
    't_rec': MetaData(
        long_name='Receiver temperature',
        units='K',
    ),
}


ATTRIBUTES_1B11 = {
    'ir_wavelength': MetaData(
        long_name='Wavelength of infrared channels',
        standard_name='sensor_band_central_radiation_wavelength',
        units='µm',
    ),
    'ir_bandwidth': MetaData(
        long_name='Bandwidth of the infrared channels central frequency',
        standard_name='sensor_band_spectral_width',
        units='µm',
        comment='channel centre frequency',
    ),
    'ir_beamwidth': MetaData(
        long_name='Full width at half maximum of the infrared channels',
        units='degrees',
    ),
    'irt': MetaData(
        long_name='Infrared brightness temperatures',
        units='K'
    ),
    'irt_accuracy': MetaData(
        long_name='Total absolute calibration uncertainty of infrared brightness temperature, one standard deviation',
        units='K',
        comment='pecify here source of this variable, e.g. literature value, specified by manufacturer, result of validation effort (updated irregularly).',
    ),
}


DEFINITIONS_1B21 = {
    'met_valid_flag':
        ('\n'
         'Bit 1: quality of air temperature data\n'
         'Bit 2: quality of relative humidity data\n'
         'Bit 3: quality of air pressure data\n'
         'Bit 4: quality of rain rate data\n'
         'Bit 5: quality of wind direction data\n'
         'Bit 6: quality of wind speed data\n'
         'Bit 7: not used\n'
         'Bit 8: not used')
}

ATTRIBUTES_1B21 = {
    'air_temperature': MetaData(
        long_name='Air temperature',
        standard_name='air_temperature',
        units='K',
    ),
    'relative_humidity': MetaData(
        long_name='Relative humidity',
        standard_name='relative_humidity',
        units='1',
    ),
    'air_pressure': MetaData(
        long_name='Air pressure',
        standard_name='air_pressure',
        units='Pa',
    ),
    'rain_rate': MetaData(
        long_name='Precipitation amount',
        standard_name='rainfall_rate',
        units='mm/h',
    ),    
    'wind_direction': MetaData(
        long_name='Wind direction',
        standard_name='wind_from_direction',
        units='degrees',
    ),    
    'wind_speed': MetaData(
        long_name='Wind speed',
        standard_name='wind_speed',
        units='m/s',
    ),     
    'met_valid_flag': MetaData(
        long_name='Meterological data validity flag',
        standard_name='met_valid_flag',
        units='1 (bit variable)',
        definition=DEFINITIONS_1B21['met_valid_flag'],
        comment='bit variable: 0=invalid, 1=valid data',
    ),    
}
