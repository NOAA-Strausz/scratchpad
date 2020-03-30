#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 16:46:44 2020

@author: strausz
"""

import argparse
import sys
import pandas as pd

parser = argparse.ArgumentParser(description='process underway data')
parser.add_argument('infile', metavar='infiles', type=str, 
                    help='input file')
args=parser.parse_args()

if args.infile:
    filename = args.infile
    col_names = ['timestamp','nitrate']
    #dtypes = {'gps_time':str,'lon':str,'lat':str,'date':str}
    df = pd.read_csv(filename, header=None, 
                     names=col_names, parse_dates=[0])
#    df['gps_datetime'] = pd.to_datetime(df.date + df.gps_time, format="%d%m%y%H%M%S")
#    
#    df['latitude'] = df.lat.str.slice(0,2).astype(int) + df.lat.str.slice(2,).astype(float) / 60
#    df['longitude'] = -1*(df.lon.str.slice(0,3).astype(int) + df.lon.str.slice(3,).astype(float) / 60)
#    
    df.set_index(['timestamp'], inplace=True)
    df_5m = df.resample('5T').mean()
    df_5m.to_csv("os1901L1_isus_5min.csv")
else:
    sys.exit("Must have input file!")

    


