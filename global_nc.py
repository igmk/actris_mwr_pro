"""Initialize Metadata of RPG MWR Level 1 + 2 variables for NetCDF file writing."""

from collections import namedtuple


FIELDS = ('description')

SiteConfig = namedtuple('SiteConfig', FIELDS)
SiteConfig.__new__.__defaults__ = (None,) * len(SiteConfig._fields)


GLOBAL_ALL = {
    'conventions': SiteConfig(
        description='Name of the conventions followed by the dataset',
    ),
    'title': SiteConfig(
        description='A succinct description of what is in the dataset, composed of instrument type and site name',
    ),
    'history': SiteConfig(
        description='Versioning of the datasets (containing date and software version)',
    ),
    'institution': SiteConfig(
        description='Where the original data was produced',
    ),
    'source': SiteConfig(
        description='The method of production of the original data',
    ),
    'comment': SiteConfig(
        description='Miscellaneous Information about the dataset or methods used to produce it',
    ),
    'references': SiteConfig(
        description='References that describe the data or methods used to produce it',
    ),
    'site_location': SiteConfig(
        description='Name of measurement station',
    ),
    'instrument_id': SiteConfig(
        description='E-PROFILE instrument identifier',
    ),
    'wigos_station_id': SiteConfig(
        description='WIGOS Station identifier acording to WIGOS convention',
    ),
    'principal_investigator': SiteConfig(
        description='Department responsible for the instrument',
    ),
    'instrument_manufacturer': SiteConfig(
        description='Manufacturer of the instrument',
    ),
    'instrument_model': SiteConfig(
        description='Instrument model',
    ),
    'instrument_generation': SiteConfig(
        description='Instrument generation',
    ),
    'instrument_hw_id': SiteConfig(
        description='Specific to mainboard',
    ),
    'network_name': SiteConfig(
        description='Name of network(s) that instrument may be part of',
    ),
    'campaign_name': SiteConfig(
        description='Name of campaign instrument may collect data for',
    ),
    'dependencies': SiteConfig(
        description='List of files the data set is depending on',
    ),
    'license': SiteConfig(
        description='Data license',
    ),
    'instrument_calibration_status': SiteConfig(
        description='Status of instrument absolute calibration',
    ),
    'date_of_last_absolute_calibration': SiteConfig(
        description='Time of last (automatic or manual) absolute calibration; LN2 or sky tipping as YYYYMMDD',
    ),
    'date_of_last_covariance_matrix': SiteConfig(
        description='Time of last covariance update as YYYYMMDD',
    ),
    'type_of_automatic_calibrations': SiteConfig(
        description='Type of automatic calibrations including information on calibration interval and respective integration time',
    ),
    'instrument_history': SiteConfig(
        description='Logbook repair/replacement work performed',
    ),
    'ir_instrument_manufacturer': SiteConfig(
        description='Manufacturer of the infrared radiometer',
    ),
    'ir_instrument_model': SiteConfig(
        description='Infrared radiometer model',
    ),
    'ir_instrument_fabrication_year': SiteConfig(
        description='Fabrication year of the infrared radiometer',
    ),
    'ir_accuracy': SiteConfig(
        description='Total absolute calibration uncertainty of infrared brightness temperature, one standard deviation. Unit: K',    
    ),
    'ir_instrument_history': SiteConfig(
        description='Logbook repair/replacement work performed',
    ),
    'met_instrument_manufacturer': SiteConfig(
        description='Manufacturer of the weather station',
    ),
    'met_instrument_model': SiteConfig(
        description='Weather station model',
    ),
    'met_instrument_fabrication_year': SiteConfig(
        description='Fabrication year of the weather station',
    ),
    'met_instrument_history': SiteConfig(
        description='Logbook repair/replacement work performed',
    ),    
    'air_temperature_accuracy': SiteConfig(
        description='Air temperature accuracy. Unit: K.',
    ),
    'relative_humidity_accuracy': SiteConfig(
        description='Relative humidity accuracy. Unit: 1.',
    ),
    'air_pressure_accuracy': SiteConfig(
        description='Air pressure accuracy. Unit: hPa.',
    ),
    'rain_rate_accuracy': SiteConfig(
        description='Rain rate accuracy. Unit: mm/h.',
    ),    
    'wind_direction_accuracy': SiteConfig(
        description='Wind direction accuracy. Unit: degrees.',
    ),
    'wind_speed_accuracy': SiteConfig(
        description='Wind speed accuracy. Unit: m/s.',
    ),      
}