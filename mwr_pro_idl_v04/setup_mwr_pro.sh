#!/bin/bash

#This script sets up the directory structure needed for MWR_PRO.
#Run this script ONLY ONCE when installing MWR_PRO. 
#Look for "USER MODIFY BEGIN" and "USER MODIFY END" to adjust to your specifications.
#Last change: U. Loehnert 16.06.2015

#USER MODIFY BEGIN

#1.) path structure - these paths must exist!
#this directory will contain your MWR data
mwr_data_path=/net/ecir/testv04
#this directory will contain executables, source code and retrieval files
home_path=/home/hatpro/testv04
#this directory will contain symbolic links of your quicklooks to be displayed in the internet
www_path=/home/hatpro/public_html

#2.) path to your IDL executeable
idl_path=/usr/local/bin/idl

#USER MODIFY END

#create necessary directorires
pro_dir=${home_path}/mwr_pro
mkdir ${pro_dir}

#create directory where processing scripts are edited & run
scripts_dir=${pro_dir}/scripts
mkdir ${scripts_dir}

#create file containing essential path information
echo ${mwr_data_path} > ${scripts_dir}/path_info.txt
echo ${idl_path} >> ${scripts_dir}/path_info.txt
echo ${home_path} >> ${scripts_dir}/path_info.txt
echo ${www_path} >> ${scripts_dir}/path_info.txt

#create directory where retrieval files (*.RET) are located
ret_dir=${pro_dir}/retrievals
mkdir ${ret_dir} 
#copy retrieval file examples to retrieval directory
cp -r retrieval_examples/*.* ${ret_dir}

#copy script templates to scripts directory
cp pl_mk_template.sh ${scripts_dir}/pl_mk.sh 

#copy par_mwr_pro.pro to scripts directory
cp par_mwr_pro.pro ${scripts_dir}

#copy tb_bias_example.sav to scripts directory
cp tb_offset_xxx.sav ${scripts_dir}

#copy sort_hatpro_raw_data.sh to scripts directory
cp sort_hatpro_raw_data.sh ${scripts_dir}

#copy uncertainty directory to scripts directory
cp -r uncertainty ${scripts_dir}

#copy source directory to mwr_pro directory
if [ -d ${pro_dir}/source ]; then
 echo 'Caution - this is NOT your first time running setup.sh!'
 echo 'All of your self-modified *.pro routines in your source directory will be overwritten if you decide to procede!'
 echo 'Backup of source directory will be created: source_backup'
 echo ' Press any key to continue.'
 cp -r ${pro_dir}/source ${pro_dir}/source_backup 
 cp -ir source ${pro_dir} 
else 
 cp -ir source ${pro_dir}
fi 

#copy color table to scripts directory
cp ${pro_dir}/source/col1 ${pro_dir}/scripts 
    
#copy setup_new_meas.sh to scripts directory
cp setup_new_meas.sh ${scripts_dir}

#copy filter_template.dat to scripts directory
cp filter_template.dat ${scripts_dir}

#copy filter_template.dat to scripts directory
cp radar_filter_template.dat ${scripts_dir}

#move example data to mwr_pro home directory
cp -r example_data ${pro_dir}

echo 'Done setting up mwr_pro!'
echo 'Now continue in your mwr_pro/scripts directory.'
