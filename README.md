# MWRpy

[![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black)

MWRpy is a Python based software to process RPG Microwave Radiometer data and is developed at the University of Cologne, Germany as part of the [Aerosol, Clouds and Trace Gases Research Infrastructure (ACTRIS)](https://actris.eu/). 
The software features reading raw data, Level 1 quality control, generation of Level 2 data products and visualization.

The data format including metadata information, variable names and file naming is designed to be compliant with the data structure and naming convention developed in the [EUMETNET Profiling Programme E-PROFILE](https://www.eumetnet.eu/). 

![MWRpy example output](https://atmos.meteo.uni-koeln.de/~hatpro/quicklooks/obs/site/jue/tophat/actris/level2/2022/10/29/20221029_juelich_temperature.png)

## MWRpy Structure

`mwrpy/rpg_mwr.py` contains the base class <b>RpgArray</b> for storing variables as netCDF4.

### `mwrpy/site_config/`

This folder contains subfolders for each site, where retrieval coeffiecients are stored in `coefficients/` and `config.py` defines site and instrument specific information (including input and output data paths), which is used for processing purposes and metadata generation.

### `mwrpy/level1/`

<b>*lev1_to_nc*</b> in `write_lev1.py` reads the raw binary files (.BRT, .BLB, .IRT, .MET, .HKD) stored in the same folder containing data of one day, applies quality control (`quality_control.py`) and writes it into a netCDF4 file using metadata defined in `lev1_meta_nc.py`.

#### Quality flags (bit variable)
    # Bit 1: missing_tb
    # Bit 2: tb_below_threshold
    # Bit 3: tb_above_threshold
    # Bit 4: spectral_consistency_above_threshold
    # Bit 5: receiver_sanity_failed
    # Bit 6: rain_detected
    # Bit 7: sun_in_beam
    # Bit 8: tb_offset_above_threshold
    
#### Level 1 Data Types
* 1B01: MWR brightnesss temperatures from .BRT and .BLB files
* 1B11: IR brightnesss temperatures from .IRT files
* 1B21: Weather station data from .MET files
* 1C01: Combined data type with time corresponding to 1B01

### `mwrpy/level2/`

<b>*lev2_to_nc*</b> in `write_lev2.py` reads Level 1 files, applies retrieval coefficients read in by `get_ret_coeff.py` for Level 2 products and writes it into a netCDF4 file using metadata defined in `lev2_meta_nc.py`. For the LWP product an offset correction is applied (`lwp_offset.py`).

#### Level 2 Data Types
* 2I01: Liquid water path (LWP)
* 2I02: Integrated water vapor (IWV)
* 2P01: Temperature profiles from single-pointing observations
* 2P02: Temperature profiles from multiple-pointing observations
* 2P03: Absolute humidity profiles
* 2P04: Relative humidity profiles (derived from 2P01/2P02 + 2P03)
* 2P07: Potential temperature (derived from 2P01/2P02 + 2P03)
* 2P08: Equivalent potential temperature (derived from 2P01/2P02 + 2P03)
* 2S02: Brightness temperature spectrum

### `mwrpy/plots/`

<b>*generate_figure*</b> in `generate_plots.py` creates .png figures using plot specific metadata defined in `plot_meta.py`.

## How to run the software

Running the software is based on a shell script (`mwrpy/actris_process_mwr.sh`) and calls `mwrpy/process_mwr_pro.py`, where functions for generating and visualizing Level 1 and Level 2 products are called (<b>*lev1_to_nc*, *lev2_to_nc*, *generate_figure*</b>).

```
mwrpy/actris_process_mwr.sh --site site_name
```
This call runs the software for the current day and reprocesses the previous day as default. The site name must be same as the subfolder name in `mwrpy/site_config/`.

Adding the option `--date YYYYMMDD` allows to process a single day and the option `--date_e YYYYMMDD` defines a date range to be processed between *date* and *date_e*.

