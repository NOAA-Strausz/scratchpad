#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  9 15:35:45 2019

@author: strausz
"""

import argparse
import numpy as np
import pandas as pd
import datetime as dt
from haversine import haversine




parser = argparse.ArgumentParser(description='Add ice concentration to drifter location')
parser.add_argument('infile', 
    metavar='infile', 
    type=str,
    help='full path to input file')

args=parser.parse_args()

latfile='/home/makushin/strausz/ecofoci_github/EcoFOCI_ssmi_ice/psn25lats_v3.dat'
lonfile='/home/makushin/strausz/ecofoci_github/EcoFOCI_ssmi_ice/psn25lons_v3.dat'
#locations of ice files
bootstrap = '/home/akutan/strausz/ssmi_ice/data/bootstrap/'
nrt = '/home/akutan/strausz/ssmi_ice/data/nrt/'
#latest available bootstrap year
boot_year = 2018


def decode_datafile(filename):
    #determine if it's nrt or bootstrap from filename prefix
    #note that we remove path first if it exists
    prefix = filename.split('/')[-1:][0][:2] 
    icefile = open(filename, 'rb')
    
    if prefix == 'nt':
        #remove the header
        icefile.seek(300)
        ice = np.fromfile(icefile,dtype=np.uint8)
        ice[ice >= 253] = 0
        ice = ice/2.5
    elif prefix == 'bt':
        ice = np.fromfile(icefile,dtype=np.uint16)
        ice = ice/10.
        ice[ice == 110] = 100 #110 is polar hole
        ice[ice == 120] = np.nan #120 is land
    else: 
        ice=np.nan
    
    icefile.close()
    
    return ice;

def get_date(filename):
    #gets date from filename
    #first remove path from filename if it is there
    filename = filename.split('/')[-1:][0]
    date = filename[3:11]
    date = dt.datetime.strptime(date,"%Y%m%d")
    return date;

def decode_latlon(filename):
    latlon_file = open(filename, 'rb')
    output = np.fromfile(latlon_file,dtype='<i4')
    output = output/100000.0
    #output = int(output * 1000)/1000 #sets decimal place at 3 without rounding
    latlon_file.close()
    return output;

def get_ice(data, df_ice):
    #
    df_ice['dist'] = df_ice.apply(lambda x: haversine((data.latitude, data.longitude), (x.latitude, x.longitude)), axis=1)

    #nearest_ice = df_ice.loc[df_ice.dist.idxmin()].ice_conc
    min_dist = df_ice.dist.min()
    if(min_dist > 17.68):
        print(min_dist)
    df_ice = df_ice[df_ice.dist<=17.68]
    #print(df_ice)
    
    return df_ice.ice_conc.mean()
    
def lon_360(lon):
    if lon < 0:
        return 360 + lon
    else:
        return lon
    


if args.infile:

    df = pd.read_csv(args.infile, sep='\s+', dtype={'year':str, 'doy':str, 'time':str} )
    
    df['combined'] = df.year + df.doy + df.time
    
    df['datetime'] = pd.to_datetime(df.combined, format='%Y%j%H%M')
   
    
    df.set_index(['datetime'], inplace=True)
    #df.drop(columns=['year','doy','time', 'combined'], inplace=True)
    df['longitude'] = df.longitude * -1
    #cleanup by resampling to the hour
    #df = df.resample('H').mean()
    df['lon_360'] = df.apply(lambda x: lon_360(x.longitude), axis=1)
    df['datetime'] = df.index
#    df.dropna(inplace=True)
#    df['argosid']=df.argosid.astype(int)
#    df['latitude']=df.latitude.round(decimals=3)
#    df['longitude']=df.longitude.round(decimals=3)
#    df['voltage']=df.voltage.round(decimals=2)
#    df['sst']=df.sst.round(decimals=2)
#    df['strain']=df.strain.round(decimals=2)
    #now group by doy
    #add blank ice column
    #df['ice_concentration'] = ''
    ice_conc = []
    groups = df.groupby(df.index.dayofyear)
    for name, group in df.groupby(df.index.dayofyear):
        #print(name)
        #print(group.latitude)
        date = group.iloc[0].datetime.strftime("%Y%m%d")
        if group.iloc[0].datetime.year <= boot_year:
            ice_file = bootstrap + str(group.iloc[0].datetime.year) + "/" + "bt_" + date + "_f17_v3.1_n.bin"
        else:
            ice_file = nrt + "nt_" + date + "_f17_nrt_n.bin"
        print("opening file: " + ice_file)
        wlon = group.lon_360.min() - .5
        elon = group.lon_360.max() + .5
        nlat = group.latitude.max() + .5
        slat = group.latitude.min() - .5
        data_ice={'latitude':decode_latlon(latfile), 'longitude':decode_latlon(lonfile),
          'ice_conc':decode_datafile(ice_file)}
        df_ice=pd.DataFrame(data_ice)
        #df_ice.dropna(inplace=True)
    
        df_ice['lon_360'] = df_ice.apply(lambda x: lon_360(x.longitude), axis=1)
        df_ice_chopped = df_ice[(df_ice.latitude < nlat) & (df_ice.latitude > slat) & (df_ice.lon_360 > wlon) & (df_ice.lon_360 < elon)]
        #print(df_ice_chopped)
        ice_conc = ice_conc + group.apply(lambda x: get_ice(x, df_ice_chopped), axis=1).to_list()
                
        
        #print("the ice concentration is: "+ice_concentration)
        #df_ice_chopped['dist'] = df_ice_chopped.apply(lambda x: haversine((data.latitude, data.longitude), (x.latitude, x.longitude)), axis=1)
    df['ice_concentration'] = ice_conc
    df=df.drop(['lon_360', 'datetime'], axis=1)
    outfile, ext = args.infile.split('.')
    outfile = outfile + "_12.5km_ice_boot.csv"
    df.to_csv(outfile)
    