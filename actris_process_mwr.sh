#!/bin/bash

#change to scripts directory
cd /home/hatpro/mwr_pro_actris/
source /home/tmarke/miniconda3/etc/profile.d/conda.sh
conda activate envi

do_plot=1

DATE_E=
DATE=
YYYYMMDD=
YYYYMMDD_E=


# Parse command-line arguments
if [ $1 = "--site" ]; then
    SITE=$2
else
    echo "Error: \"$1\" not understood"
    exit 1        
fi
if [ $# -gt 2 ]; then
    if [ $3 = "--date" ]; then
        if [ $# -lt 4 ]; then
            echo "No date given"
            exit 1
        fi
        DATE=$4
        if [ $# -gt 4 ]; then
           if [ $5 = "--date_e" ]; then
               if [ $# -lt 6 ]; then
                  echo "No end date given"
                  exit 1
               fi
               DATE_E=$6
           fi
        else
           DATE_E=$4
        fi
    else
        echo "Error: \"$3\" not understood"
        exit 1
    fi
fi


if [ ! "$DATE" ] ; then
    # Default option is to get the day before yesterdays data 
    YYYYMMDD=$(date --utc "+%Y%m%d" -d "today-2day")
    YYYYMMDD_E=$(date --utc "+%Y%m%d" -d "today-1day")
else
    # Check that date is YYYYMMDD (very rudimentary!)
    if [ $(echo "$DATE" | tr -d " "| wc -L) -eq 8 ] ; then
        YYYYMMDD=$DATE
    else
        echo "Error: DATE \"$DATE\" not understood"
        exit 1
    fi
    # Check that date is YYYYMMDD (very rudimentary!)
    if [ $(echo "$DATE_E" | tr -d " "| wc -L) -eq 8 ] ; then
        YYYYMMDD_E=$(date "+%Y%m%d" -d "$DATE_E + 1 day")
    else
        echo "Error: DATE \"$DATE_E\" not understood"
        exit 1
    fi

fi   


while [ "$YYYYMMDD" != "$YYYYMMDD_E" ]; do
	echo $YYYYMMDD

	python3 process_mwr_pro.py $SITE $YYYYMMDD $do_plot
	
	YYYYMMDD=$(date "+%Y%m%d" -d "$YYYYMMDD + 1 day")

done

conda deactivate

exit
