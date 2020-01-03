#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 26 15:24:06 2019

@author: strausz
"""

import numpy as np
import pandas as pd
import datetime as dt
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.io import shapereader
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER
from scipy import interpolate 
#from matplotlib.mlab import griddata
import cmocean
import sys
import argparse
import matplotlib.ticker as mticker
import cartopy.mpl.ticker as cticker
import xarray as xr
from haversine import haversine

parser = argparse.ArgumentParser(description='Decode binary SSMI satellite data')
parser.add_argument('infile', 
    metavar='infile', 
    type=str,
    help='full path to input file')
parser.add_argument('-m', '--mooring', nargs=1, 
                    help='add mooring location, ie ck1-ck4 or bs2-bs8 or all')
parser.add_argument('-ex', '--extents', nargs=1,
                    help='chooses extents of map, options are bering, chukchi, custom, and default')
parser.add_argument('-p', '--plot', action="store_true",
                    help='plot the ice concentrations on a map')
parser.add_argument('-c', '--csv', action="store_true",
                    help='make a CSV file of the data with lat,lon, and ice percent')


args=parser.parse_args()


#data_file=sys.argv[1]
#data_file='nt_20180402_f18_nrt_n.bin'
#these are the files that contain the lats and lons. Obtained from here:
#ftp://sidads.colorado.edu/pub/DATASETS/seaice/polar-stereo/tools/
latfile='/home/makushin/strausz/ecofoci_github/EcoFOCI_ssmi_ice/psn25lats_v3.dat'
lonfile='/home/makushin/strausz/ecofoci_github/EcoFOCI_ssmi_ice/psn25lons_v3.dat'

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
        ice[ice == 110] = 100 #110 is the polar hole, assume 100% ice
        ice[ice == 120] = np.nan #120 should be land
    else: 
        ice=np.nan
    
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
    return output;

if args.infile:        
    ice = decode_datafile(args.infile).reshape((448,304))
    lat = decode_latlon(latfile).reshape((448,304))
    lon = decode_latlon(lonfile).reshape((448,304))
#    ds = xr.Dataset({'ice':(['x','y'], ice)}, coords={'lon':(['x', 'y'], lon),
#                     'lat': (['x','y'],lat)})
#    #lets try just using a pandas dataframe
#    data={'latitude':decode_latlon(latfile), 'longitude':decode_latlon(lonfile),
#          'ice_conc':decode_datafile(args.infile)}
#    df=pd.DataFrame(data)
    ice2 = ice[35:245, 15:156]
    lat2=lat[35:245, 15:156]
    lon2=lon[35:245, 15:156]
    data={'latitude':lat2.flatten(), 'longitude':lon2.flatten(),'ice_conc':ice2.flatten()}
    df=pd.DataFrame(data)
    df.to_csv('bootstrap_2018_12_31_trimmed_to_bcb.csv')