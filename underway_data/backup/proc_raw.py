#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 16:46:44 2020

@author: strausz
"""

#pulls all underway system raw files together into one dataframe
#must have combined all files for each instrument into one file first
#easy and quick to do with perl one liners


import pandas as pd


instruments = {'gps':'raw_combined/os1901_gps_full_gprmc_only.csv', 
               'isus':'raw_combined/os1901_isus_full_lf_only.csv',
               'eco':'raw_combined/os1901_eco_full.csv',
               'optode':'raw_combined/os1901_optode_full.csv',
               'sbe48':'raw_combined/os1901_sbe48_full.csv',
               'tdgp':'raw_combined/os1901_tdgp_full.csv',
               'tsg':'raw_combined/os1901_tsg_full.csv'}


#
#    filename = args.infile
#    col_names = ['timestamp','nmea_code','gps_time','status','lat','lat_dir',
#             'lon','lon_dir','speed','hdg','date','mag','mag_dir','dcheksum']
#    dtypes = {'gps_time':str,'lon':str,'lat':str,'date':str}
#    df = pd.read_csv(filename, header=None, sep=',|\*', dtype=dtypes, 
#                     names=col_names, parse_dates=[0])
#    df['gps_datetime'] = pd.to_datetime(df.date + df.gps_time, format="%d%m%y%H%M%S")
#    
#    df['latitude'] = df.lat.str.slice(0,2).astype(int) + df.lat.str.slice(2,).astype(float) / 60
#    df['longitude'] = -1*(df.lon.str.slice(0,3).astype(int) + df.lon.str.slice(3,).astype(float) / 60)
#    
#    df.set_index(['timestamp'], inplace=True)
#    df_5m = df.resample('5T').mean()
#    df_5m.to_csv("os1901L1_gps_5m.csv")

def intake(instrument, filename):
    dtypes={}
    sep=','
    if instrument == 'gps':
        col_names = ['timestamp','gps_time','lat','lon','speed','hdg','date','mag']
        dtypes = {'gps_time':str,'lon':str,'lat':str,'date':str}
        usecols = [0,2,4,6,8,9,10,11]
    elif instrument == 'isus':
        col_names = ['timestamp','nitrate(ug/l)']
        usecols = [0,4]
    elif instrument == 'eco':
        col_names = ['timestamp','fl_counts']
        sep = ",|\s+"
        usecols = [0,4]
    elif instrument == 'optode':
        col_names = ['timestamp','O2_conc(uM)','O2_percent_sat','optode_temp']
        sep = ",|\s+"
        usecols = [0,3,4,5]
    elif instrument == 'sbe48':
        col_names = ['timestamp','hull_temp']
        usecols = [0,1]
    elif instrument == 'tdgp':
        col_names = ['timestamp','dgp(mb)','tdgp_temp']
        usecols = [0,7,8]
    elif instrument == 'tsg':
        col_names = ['timestamp','tsg_temp','cond(S/m)','sal']
        usecols = [0,1,2,3]
    df = pd.read_csv(filename, header=None, sep=sep,
                     names=col_names, parse_dates=[0], usecols=usecols, dtype=dtypes)
    df.set_index(['timestamp'], inplace=True)
    #df = df.resample('1T').mean()
    return df
dfs = {}
for key, value in instruments.items():
    dfs[key] = intake(key, value)

#fix the gps lats and lons
dfs['gps']['gps_datetime'] = pd.to_datetime(dfs['gps'].date + dfs['gps'].gps_time, format="%d%m%y%H%M%S")
dfs['gps']['latitude'] = dfs['gps'].lat.str.slice(0,2).astype(int) + dfs['gps'].lat.str.slice(2,).astype(float) / 60
dfs['gps']['longitude'] = -1*(dfs['gps'].lon.str.slice(0,3).astype(int) + dfs['gps'].lon.str.slice(3,).astype(float) / 60)
dfs['gps']['deltatime'] = (dfs['gps'].index - dfs['gps'].gps_datetime).dt.total_seconds()

#convert fluorometer
dfs['eco']['chl_conc(ug/l)'] = 0.0063*(dfs['eco'].fl_counts - 112)

#this resamples asll the instruments to 1min data
for key in instruments:
    if key == 'isus':
        #the isus made 11 1s samples every 3 minutes  resampling to 3 minutes
        #averages those values every 3 min, then interpolating fills in
        #some of the missing data.  Resampling again to 1m interval leaves
        #gaps so I interpolated again to get 1min data
        dfs[key] = dfs[key].resample('3T').mean()
        dfs[key] = dfs[key].interpolate(method='linear', limit=2)
        dfs[key] = dfs[key].resample('1T').mean()
        dfs[key] = dfs[key].interpolate(method='linear', limit=2)
    else:
        #the other data all has sample rates more than 1m so extra steps are
        #not needed
        dfs[key] = dfs[key].resample('1T').mean()
    

df = pd.concat(dfs.values(), axis=1, sort=False)

df_rounded = df.round({'speed':1, 'hdg':0, 'latitude':4, 'longitude':4,
                       'nitrate(ug/l)':1, 'chl_conc(ug/l)':2, 'O2_conc(uM)':2, 'O2_percent_sat':2,
                       'optode_temp':2, 'hull_temp':3, 'dgp(mb)':2, 'tdgp_temp':2,
                       'tsg_temp':3, 'cond(S/m)':3, 'sal':3})
df_rounded.drop(columns=['mag', 'deltatime', 'fl_counts'],inplace=True)



