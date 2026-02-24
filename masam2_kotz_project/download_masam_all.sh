#!/bin/bash

base_url="https://noaadata.apps.nsidc.org/NOAA/G10005_V2/Data/"
file_prefix="masam2_minconc40_"
file_suffix="_v2.nc"
for year in {2012..2025}
do
  mkdir -p $year
  cd $year
  month=1
  until [[ $month -gt 12 ]]
  do 
    month_padded=$(printf '%02d' $month)
    file=$base_url${year}/${file_prefix}${year}${month_padded}${file_suffix}
    echo "Downloading $file"
    wget $file
    month=$(($month+1))
  done
  cd ../
done
