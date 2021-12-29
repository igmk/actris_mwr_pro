#!/bin/bash

#Run this script if you are installing processing for a new mwr. 
#Look for "USER MODIFY begin" and "USER MODIFY end" to adjust to your specifications.
#Last change: U. Loehnert 16.06.2015

#USER MODIFY begin 

#1.) radiometer measurement (user-defined three letter code)
mwr_meas=xxx

#USER MODIFY end

#get path parameters from path_info.txt
declare -a ARRAY

for LINE in `cat path_info.txt`; do 
 ARRAY[$count]=$LINE
 ((count++))
done
  
mwr_data_path=${ARRAY[0]}
  
#create directory for new instrument
mwr_dir=${mwr_data_path}/${mwr_meas}
mkdir ${mwr_dir}

#create data directory for final data in netcdf format
mkdir ${mwr_dir}/data
#sub-directory for voltages
mkdir ${mwr_dir}/data/level0
#physical measurements (e.g. TB)
mkdir ${mwr_dir}/data/level1
#products (IWV, LWP; WDL, T, q, rh, ...)
mkdir ${mwr_dir}/data/level2
#raw data (TBs)
mkdir ${mwr_dir}/data/raw
#calibration data
mkdir ${mwr_dir}/data/calibration
#temporary data
mkdir ${mwr_dir}/data/temp

#create data directory for quicklooks
mkdir ${mwr_dir}/plots
#voltagess
mkdir ${mwr_dir}/plots/level0
#TBs
mkdir ${mwr_dir}/plots/level1
#products
mkdir ${mwr_dir}/plots/level2

#move error and uncertainty specifications to data directory
mv uncertainty ${mwr_dir}/data/

#copy manual filter file template to mwr directory
if [ -f ${mwr_dir}/filter_${mwr_meas}.dat ]; then
 echo 'Be careful, existing filter file will be overwritten if you decide to proceed!'
 echo 'Backup file will be created (filter.backup)'
 cp ${mwr_dir}/filter_${mwr_meas}.dat ${mwr_dir}/filter.backup 
 cp -i filter_template.dat ${mwr_dir}/filter_${mwr_meas}.dat
 cp radar_filter_${mwr_meas}.dat radar_filter.backup 
 cp -i radar_filter_template.dat radar_filter_${mwr_meas}.dat
else
 cp -i filter_template.dat ${mwr_dir}/filter_${mwr_meas}.dat
 cp -i radar_filter_template.dat radar_filter_${mwr_meas}.dat
fi

#comment out the next lines if you wish NOT to copy some example data from Juelich to your MWR data directory.
#If you do execute the follwing lines, you may run pl_mk.sh and the netcdf files and quicklooks
#are produced for September 18, 2012

mv ../example_data/mwr/* ${mwr_dir}/data/raw                                                         

echo 'Now modify "data_path_radar", "filter_file_radar" and "ceilo_path" manually in par_mwr_pro.pro.'
echo 'Example radar and ceilometer data have been copied to your mwr_pro/example_data directory.'
echo 'Then execute pl_mk.sh to see if everything works properly.'
echo 'Afterwards: check data and quicklook directories if everything went smoothely...'
