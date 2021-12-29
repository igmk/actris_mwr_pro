#!/bin/bash
#Please the header of this file carefully!
#Running this MWR_PRO script will produce daily netcdf files and plots of
#1.) level1 (measurement products (TBs))
#2.) level2a (integrated retrieval products)
#(needs lwp, iwv, zwd or att retrievals)
#3.) level2b&c (T&q profile retrieval products)
#(needs hze, tze or tel retrievals)
#from raw radiometer data on a monthly basis. You can process (test-wise) only one day by setting the variable "day_files" 
#to one specific day (see below).
#--> Look for ####SPECIFY begin and ####SPECIFY end in this file to specify the dates you would like to process.
#--> You must specify the retrieval products and the rest of processing options
#    in the file par_mwr_pro.pro in this directory. Getting aquainted with par_mwr_pro.pro is imperative for succesful MWR data processing!

####SPECIFY begin
#3 letter code for identification of measurement setup
mwr_meas=xxx
#times to process:
#set the following 1 line if you want to reprocess certain months (mmdd)
#months=(1301 1302 1303 1304 1305 1306 1307 1308 1309 1310 1311 1312)
months=(1209)
#set the following 4 lines if you want to process only the current day (online-processing mode)
#tm=$(date +%y%m)
#thisday=$(date +%y%m%d)
#yesterday=$(date --date 'yesterday' +%y%m%d)
##months=$tm
####SPECIFY end

#copy parameter file to source directory where it will be executed
cp 'par_mwr_pro.pro' '../source/par_mwr_pro.pro'

#specify add mwr_pro directory to IDL path
IDL_PATH="<IDL_DEFAULT>"
path_idl="!path +': ../source'"

#get path parameters from path_info.txt
declare -a ARRAY

for LINE in `cat path_info.txt`; do 
 ARRAY[$count]=$LINE
 ((count++))
done

mwr_data_path=${ARRAY[0]}
idl_exec=${ARRAY[1]}      
home_path=${ARRAY[2]}
www_path=${ARRAY[3]}/${mwr_meas}
mwr_dir=${mwr_data_path}/${mwr_meas}

# check whether month directories exist, if not create
for month in ${months[*]}; do

 if [ ! -d ${mwr_dir}/data/level0/${month} ]; then
  mkdir -p ${mwr_dir}/data/level0/${month}
  mkdir -p ${mwr_dir}/data/level1/${month}
  mkdir -p ${mwr_dir}/data/level2/${month}
 fi
 
 if [ ! -d ${mwr_dir}/plots/level0/${month}  ]; then
  mkdir -p ${mwr_dir}/plots/level0/${month}
  mkdir -p ${mwr_dir}/plots/level1/${month}
  mkdir -p ${mwr_dir}/plots/level2/${month}
  mkdir -p ${www_path}/TB/${month}
  mkdir -p ${www_path}/IWV_LWP/${month}
  mkdir -p ${www_path}/VAPs/radar_mwr/${month}
  mkdir -p ${www_path}/VAPs/scans/${month} 
  mkdir -p ${www_path}/VAPs/Tq_profiles/${month}
  mkdir -p ${www_path}/VAPs/cili/${month}
 fi

 data_path=${mwr_dir}/data/raw/${month}
 plot_path0=${mwr_dir}/plots/level0/${month}
 plot_path1=${mwr_dir}/plots/level1/${month}
 plot_path2=${mwr_dir}/plots/level2/${month}

####SPECIFY begin
#set the following line if you would like to process only current day
#day_files=`echo $yesterday $thisday`
#set day_files to data path if you are reprocessing a whole month
day_files=$(ls ${data_path})
#set day_files "by hand" to process only this day (e.g. testing purpose)
#day_files=(130901)
#day_files=(130729)
####SPECIFY end
 for day_file in $day_files; do
  day_file=`echo $day_file | cut -c -6`
  echo $day_file

#copy parameter file to source directory where it will be executed
  cp par_mwr_pro.pro source/par_mwr_pro.pro

#call pl_mk.pro
  echo 'calling pl_mk for: '$day_file
  pl_mk_call='pl_mk, date='"'"${day_file}"'"',mwr_meas='"'"${mwr_meas}"'"',mwr_dir='"'"${mwr_dir}"'"',home_path='"'"${home_path}"'"
  echo "!path="${path_idl} > @l0_${mwr_meas}
  echo ${pl_mk_call} >> @l0_${mwr_meas} 
  echo 'exit' >> @l0_${mwr_meas}
  ${idl_exec} @l0_${mwr_meas}
  rm @l0_${mwr_meas}

#convert ps-files file to png-files:
#Please comment out lines, where no plots are produced
#1.) daily time series tbs
  ps_file=${plot_path1}/${day_file}_${mwr_meas}_l1.ps
  if [ -f ${ps_file} ]; then
   png_file=${plot_path1}/${day_file}_${mwr_meas}_l1.png
   echo 'creating: ' ${png_file} 
   pstoimg -flip r90 -density 200 ${ps_file}
   rm ${ps_file}
   ln -sf ${png_file} ${www_path}/TB/${month}
  fi
#2.) mean daily tb spectrum
  ps_file=${plot_path1}/${day_file}_${mwr_meas}_l1_sp.ps
  if [ -f ${ps_file} ]; then
   png_file=${plot_path1}/${day_file}_${mwr_meas}_l1_sp.png
   echo 'creating: ' ${png_file} 
   pstoimg -flip r90 -density 200 ${ps_file}
   rm ${ps_file}
   ln -sf ${png_file} ${www_path}/TB/${month}
 fi
#3.) daily flag time series
  ps_file=${plot_path1}/${day_file}_${mwr_meas}_flag.ps
  if [ -f ${ps_file} ]; then
   png_file=${plot_path1}/${day_file}_${mwr_meas}_flag.png
   echo 'creating: ' ${png_file} 
   pstoimg -flip r90 -density 200 ${ps_file}
   rm ${ps_file}
   ln -sf ${png_file} ${www_path}/TB/${month}
 fi
#4.) daily time series lwp/iwv/wdl
  ps_file=${plot_path2}/${day_file}_${mwr_meas}_l2a.ps
  if [ -f ${ps_file} ]; then  
   png_file=${plot_path2}/${day_file}_${mwr_meas}_l2a.png 
   echo 'creating: ' ${png_file}       
   pstoimg -flip r90 -density 200 ${ps_file}                          
   rm ${ps_file}                      
   ln -sf ${png_file} ${www_path}/IWV_LWP/${month}
  fi                 
#5.) daily time series of iwv/lwp azimuth-time contours
  ps_file=${plot_path2}/${day_file}_${mwr_meas}_l2a_iwv_aztp.ps
  if [ -f ${ps_file} ]; then
   png_file=${plot_path2}/${day_file}_${mwr_meas}_l2a_iwv_aztp.png
   echo 'creating: ' ${png_file}
   pstoimg -flip r90 -density 200 ${ps_file}
   rm ${ps_file}
   ln -sf ${png_file} ${www_path}/VAPs/scans/${month}
  fi
  ps_file=${plot_path2}/${day_file}_${mwr_meas}_l2a_lwp_aztp.ps
  if [ -f ${ps_file} ]; then
   png_file=${plot_path2}/${day_file}_${mwr_meas}_l2a_lwp_aztp.png
   echo 'creating: ' ${png_file}
   pstoimg -flip r90 -density 200 ${ps_file}
   rm ${ps_file}
   ln -sf ${png_file} ${www_path}/VAPs/scans/${month}
  fi
#6.) daily time series of IRT, ceilometer, LWP/IWV
  ps_file=${plot_path2}/${day_file}_${mwr_meas}_l2_cili.ps
  if [ -f ${ps_file} ]; then
   png_file=${plot_path2}/${day_file}_${mwr_meas}_l2_cili.png
   echo 'creating: ' ${png_file}
   pstoimg -flip r90 -density 200 ${ps_file}
   rm ${ps_file}
   ln -sf ${png_file} ${www_path}/VAPs/cili/${month}
  fi
#7.) daily time series 2b profiles (T and q along line of sight)
  ps_file=${plot_path2}/${day_file}_${mwr_meas}_l2b_tze.ps
  if [ -f ${ps_file} ]; then  
   png_file=${plot_path2}/${day_file}_${mwr_meas}_l2b_tze.png 
   echo 'creating: ' ${png_file}       
   pstoimg -flip r90 -density 200 ${ps_file}                          
   rm ${ps_file}                      
  fi                 
  ps_file=${plot_path2}/${day_file}_${mwr_meas}_l2b_hze.ps
  if [ -f ${ps_file} ]; then  
   png_file=${plot_path2}/${day_file}_${mwr_meas}_l2b_hze.png 
   echo 'creating: ' ${png_file}       
   pstoimg -flip r90 -density 200 ${ps_file}                          
   rm ${ps_file}                      
  fi                 
#8.) daily time series 2c profiles (T, Tpot, qiwv, rhlwp)
  ps_file=${plot_path2}/${day_file}_${mwr_meas}_l2c_t.ps
  if [ -f ${ps_file} ]; then  
   png_file=${plot_path2}/${day_file}_${mwr_meas}_l2c_t.png 
   echo 'creating: ' ${png_file}       
   pstoimg -flip r90 -density 200 ${ps_file}                          
   rm ${ps_file}                      
   ln -sf ${png_file} ${www_path}/VAPs/Tq_profiles/${month}
  fi                 
  ps_file=${plot_path2}/${day_file}_${mwr_meas}_l2c_tpot.ps
  if [ -f ${ps_file} ]; then  
   png_file=${plot_path2}/${day_file}_${mwr_meas}_l2c_tpot.png 
   echo 'creating: ' ${png_file}       
   pstoimg -flip r90 -density 200 ${ps_file}                          
   rm ${ps_file}                      
   ln -sf ${png_file} ${www_path}/VAPs/Tq_profiles/${month}
  fi                 
  ps_file=${plot_path2}/${day_file}_${mwr_meas}_l2c_qiwv.ps
  if [ -f ${ps_file} ]; then  
   png_file=${plot_path2}/${day_file}_${mwr_meas}_l2c_qiwv.png 
   echo 'creating: ' ${png_file}       
   pstoimg -flip r90 -density 200 ${ps_file}                          
   rm ${ps_file}                      
   ln -sf ${png_file} ${www_path}/VAPs/Tq_profiles/${month}
  fi                 
  ps_file=${plot_path2}/${day_file}_${mwr_meas}_l2c_rhlwp.ps
  if [ -f ${ps_file} ]; then  
   png_file=${plot_path2}/${day_file}_${mwr_meas}_l2c_rhlwp.png 
   echo 'creating: ' ${png_file}       
   pstoimg -flip r90 -density 200 ${ps_file}                          
   rm ${ps_file}                      
   ln -sf ${png_file} ${www_path}/VAPs/Tq_profiles/${month}
  fi
#9.) daily time series radar measurements and 2a LWP              
  ps_file=${plot_path2}/${day_file}_${mwr_meas}_l2_radar.ps
  if [ -f ${ps_file} ]; then  
   png_file=${plot_path2}/${day_file}_${mwr_meas}_l2_radar.png 
   echo 'creating: ' ${png_file}       
   pstoimg -flip r90 -density 200 ${ps_file}                          
   rm ${ps_file}                      
   ln -sf ${png_file} ${www_path}/VAPs/radar_mwr/${month}
  fi

  ps_file=${plot_path2}/${day_file}_${mwr_meas}_l2_radar_lwp.ps
  if [ -f ${ps_file} ]; then  
   png_file=${plot_path2}/${day_file}_${mwr_meas}_l2_radar_lwp.png 
   echo 'creating: ' ${png_file}       
   pstoimg -flip r90 -density 200 ${ps_file}                          
   rm ${ps_file}                      
   ln -sf ${png_file} ${www_path}/VAPs/radar_mwr/${month}
  fi
 done
done

