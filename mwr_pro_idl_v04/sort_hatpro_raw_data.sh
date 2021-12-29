#!/bin/bash
############################################################
#   This program creates day directories for the hatpro data
#   according to MWR_PRO (yymm/yymmdd) 
#   replacing RPG Structure (Yyyyy/Mmm/Ddd)
#   upper to lower case characters
#         created by Bernhard Pospichal 08-01-2012
#         changed 10-01-2012
############################################################

################## make specifications #####################

# specify year and month
year=11
months=(02 03)

# specify your paths
station=ift

# specify your path of RPG's R2CH-Software generated archive 
# with the following structure: (e.g. Y2011/M12/D24)
rpg_temp_dir=/projekt1/mikrowelle/data_hatpro/$station/raw

# specify the path for your raw data directory 
# for further processing with mwr_pro software package
raw_dir=/projekt1/mikrowelle/data_hatpro/$station/test

################# end personal specifications ###############


# change to archive
cd $rpg_temp_dir/Y20$year

if [ ! -d $rpg_temp_dir/Y20$year ]; then 
   exit
fi

# collect the data in RPG format
# loop over month directories specified above
for month in  ${months[*]}

  do
    j=M$month
    echo $j
  # check whether monthly directories exist. if not create them
    if [ ! -d $raw_dir/$year$month  ]; then
       mkdir -p $raw_dir/$year$month
    fi

  # change to month directory (Mxx)
    cd $j

    # loop over day directories
    for k in `ls`
      do 
        echo $k
        # check whether daily directories exist. if not create them
        if [ ! -d $raw_dir/$year$month/$year$month${k:1:2}  ]; then
 
          mkdir -p $raw_dir/$year$month/$year$month${k:1:2}
        fi

        # change to day directory (Dxx)
        cd $k
       
        # loop over all files in day directory
        for i in `ls`
          do
          
            # select only files with the following endings (currently BRT, HKD, IRT, MET)
            tag=`expr "$i" : '.*\(...\)'` 
            if [ "$tag" = "BRT" -o "$tag" = "HKD" -o "$tag" = "IRT" -o "$tag" = "MET" ]; then
              echo "Copying $i"
              # move or copy the data by translating upper to lower case characters 
              # mv $i  $raw_dir/$year$month/$year$month${k:1:2}/`echo "$i" | tr '[A-Z]' '[a-z]'`
              cp $i  $raw_dir/$year$month/$year$month${k:1:2}/`echo "$i" | tr '[A-Z]' '[a-z]'`
            fi     

          done
        cd .. 
     done
  cd ..
  
done
