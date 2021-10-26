## actris_mwr_pro

This Python/Jupyter-notebook code to process ground based Microwave Radiometer 
data is developed within the Aerosol Clouds and Trace Gases Researge 
Infrastructure [ACTRIS project](https://actris.eu/).
The program processes the measurement data ,starting from the binary data 
comming from the instrument, and produces quality controlled Level 1 output
and Level 2 products (LWP, IWV and temperature profiles). 

The code is based on the scientific knowledge builed at the Institute for 
Geophysic and Meteorologie (IGMK) at the University of Cologne, Cologne, 
Germany. The program was devdeloped and is mentained by Tobias Marke and 
Lukas Pfitzenmaier in cooperation with Bernhard Pospichal and Ulrich Loehnert. 


## Content of this documentation: ##

0. General stucture

1. How to call the program

2. Setting data processing options (config file)

3. Input and output files
    3.1. Supported input files    
    3.2. Output files: naming convention and types
    
4. Short description of the data processing
    4.1 Level 1 - processing and data filter
    4.2 Level 2 - processing and data filter
    
 X. References


####################### DESCRIPTION OF THE PROCESSING #####################


# 0. General stucture #

All code for this program is contained in direcotry scripts, and its 
subdirectories. It is assuemd that the "folders" containing the functions
are located in the same directory where you execute the program. 
For Level 1 the "site_config.py" the Site Specific Information have to be 
set that the output contains all relefant golbal atrigutes (Level1 > 
site_config.py >get_site_specs).
The data "input path" has to be set in the XYZ in "Level1".
The output format of the Level 1 processing can be selected - for a 
description see Section 3.1.  


# 1. How to call the program #

Call the scripts with:
- Run Level 1 processing - MWR_PRO()
    to process today's and yesterdays data
- Run Level 2 processing - MWR_PRO()
    to process the given date


# 2. Setting data processing options (config file) #

The settings required for running the program have to be set in the 
configuration file (config_radarname, see above for calling the program). 
This file includes:
- input and output directories of the data
- Level 1 output options (Level 1 data files comparable to the sigle 
  input files, daily file contain all quality controlled data)
- Level 2 processing options:
    - which 
    - overwriting of output file(s)
    - debugging
    - definition of output file type
    - whether moments are calculated from spectra or copied from RPG lv1 
      file (this option is not available at the moment 10.7.2019)

See config_example.m for details.


# 3. Input and output files #

## 3.1 Supported input files ##

The program can handles the binary-files from RPG-HatPro:
    i) binary files of type 1 created with RPG radar software version 1.
    ii) netcdf files created from i) where no additional processing 
        has been applied. 
    iii) binary files of type 2 created with RPG radar software version 2. 
    iv) binary files of type 3 created with RPG radar software version 3-5.
This program automatically identifies the file type and adjusts reprocessing
and creates a unified netcdf file for any of the file types i) to iv). 
Input data directory defined in config file (config.datapath). Data is 
assumed to be located in config.datapath/yyyy/mm/dd/
Note that you need to have read access to the directory.
The option of reading RPG lv1 binary files and converting them to netcdf is
not working at the moment (10.7.2019).


## 3.2 Output files: Naming convention and types ##

Two types of putput files are available. 
    i) general file: includes all metadata information, all flags, 
       all spectra, all moments
    ii) compact file: includes only moments (no spectra), some metadata 
Set in the config file if i), ii), or both should be created.

Naming convention for output files:
i) radarname_station_yyyymmddHHMMSS_program_scan.nc,
    example: joyrad94_nya_20190710000000_P01_ZEN.nc
ii) radarname_station_yyyymmddHHMMSS_program_scan_compact.nc,
    example: joyrad94_nya_20190710000000_P01_ZEN_compact.nc
where radarname and station are specified in config file, program refers to 
the program number (chirp table) defined in the RPG radar software, and 
scan is the scanning strategy also set in the RPW radar software. 'program' 
and 'scan' are copied from the input file name.

Output data directory defined in config file (config.outputpath). Data will
be written in config.outputpath/yyyy/mm/dd/*.nc If directories are not 
existing they will be created by the program. Note that you need to have 
write access to the target directory.


# 4. Short description of the whole data processing (input -> Level 1 -> Level 2) #
