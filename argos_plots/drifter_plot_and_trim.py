#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 16:00:10 2019

@author: strausz
"""

import pandas as pd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
#import cartopy.feature as cfeature
import cmocean
#for calculating distance
from haversine import haversine
import argparse
from datetime import datetime



parser = argparse.ArgumentParser(description='Plot drifter track on map')
parser.add_argument('infile', type=str, help='full path to input file')
parser.add_argument('-p', '--plot', type=str, 
                    help="make plot of 'sst', 'strain', or 'speed' ")
parser.add_argument('-f', '--file', action="store_true",
                    help="output csv file of data")
parser.add_argument('-ph', '--phyllis', action="store_true",
                    help="output format for phyllis friendly processing")
parser.add_argument('-d', '--date', action="store_true",
                    help="add occasional date to track")
parser.add_argument('-e', '--erddap', action="store_true",
                    help="input file has erddap csv format")
parser.add_argument('-c', '--cut', nargs=2, 
                    type=lambda x: datetime.strptime(x, "%Y-%m-%dT%H:%M:%S"),
                    help="date span in format '2019-01-01T00:00:00 2019-01-01T00:00:01' for example")
args=parser.parse_args()




filename=args.infile

def get_extents(df):
    nlat = df.latitude.max() + 2
    slat = df.latitude.min() - 2
    wlon = df.longitude.max() - 20
    elon = df.longitude.min() + 20

    extents = [wlon, elon, slat, nlat]
    return extents

def plot_variable(dfin, var):
    proj = ccrs.LambertConformal(central_longitude=-165, central_latitude=60)
    ax = plt.axes(projection=proj)
    ax.natural_earth_shp(name='land', resolution='50m' )
    ax.coastlines(resolution='50m')
    if var == 'speed' :
        vmin, vmax, cmap = 0, 120, cmocean.cm.speed
    elif var == 'sst':
        vmin, vmax, cmap = -2, 20, cmocean.cm.thermal
    elif var == 'strain':
        vmin, vmax, cmap = 0, 20, cmocean.cm.haline
    plotted = ax.scatter(dfin['longitude'], dfin['latitude'], s=10, c=dfin[var], transform=ccrs.PlateCarree(), 
               cmap=cmap, vmin=vmin, vmax=vmax )
    plt.colorbar(plotted)
    #ax.plot(dfin['longitude'], dfin['latitude'], transform=ccrs.PlateCarree())
    ax.set_extent(get_extents(dfin))
    title = str(dfin.argosid[0]) + " " + var
    ax.set_title(title)
    
    return ax

def trim_data(df, delta_t):
    start = delta_t[0].strftime('%Y-%m-%d %H:%M:%S')
    end   = delta_t[1].strftime('%Y-%m-%d %H:%M:%S')
    return df[start:end]


if args.erddap:
    names = ['argosid','strain','voltage','datetime','latitude','sst','longitude']
    df=pd.read_csv(filename, header=0, names=names, parse_dates=[3])
    df['longitude'] = df.longitude - 360
    df['datetime'] = df.datetime.dt.tz_localize(None) #to remove timezone info
    df.set_index(['datetime'], inplace=True)
else:
    df = pd.read_csv(filename)
    
    df['datetime'] = pd.to_datetime(df['year_doy_hhmm'], format='%Y-%m-%d %H:%M:%S')
    df['datetime'] = df.datetime.dt.tz_localize(None) #to remove timezone info
    
    df.set_index(['datetime'], inplace=True)
    df.drop(columns=['year_doy_hhmm','year_doy_hhmm.1'], inplace=True)
    df['longitude'] = df.longitude * -1

if args.cut:
    df = trim_data(df, args.cut)

#resample data to on an even hour
df_hour=df.resample('H').mean()

#use linear interpolation to fill in gaps
df_hour.interpolate(inplace=True)

#now calculate distance for drifter speed calculation
df_hour['time'] = df_hour.index
df_hour['next_lat'] = df_hour.latitude.shift(-1)
df_hour['next_lon'] = df_hour.longitude.shift(-1)
#then calculate distance between points with the haversine function
df_hour['dist'] = df_hour.apply(lambda x: haversine((x.latitude, x.longitude), (x.next_lat, x.next_lon)), axis=1)
#next shift up the 'dist' column
df_hour.dist.shift()
#now calculate the time difference

df_hour['time2'] = df_hour.time.shift(-1)
#make the time_delta
df_hour['time_delta'] = df_hour.time2 - df_hour.time
#now make new column of seconds
df_hour['seconds'] = df_hour.time_delta.dt.total_seconds()
#now calculate speed in m/s
df_hour['speed'] = df_hour.dist * 100000 / df_hour.seconds
df_hour['argosid'] = df_hour.argosid.astype(int)

#now can do plotting stuff if selected

if args.plot:
    ax = plot_variable(df_hour, args.plot)
#    ax.scatter(df_hour['longitude'], df_hour['latitude'], s=10, c=df_hour.speed, transform=ccrs.PlateCarree(), 
#               cmap=cmocean.cm.speed, vmin=0, vmax=120 )


if args.file:
    df_out = df_hour[['argosid','latitude','longitude','sst','strain','voltage','speed']]
    df_out = df_out.round({'latitude':3, 'longitude':3,'sst':2,'strain':1,'voltage':1,'speed':1})
    outfile = str(df_out.argosid[0]) + '_trimmed.csv'
    df_out.to_csv(outfile)
    
    
if args.phyllis:
    df_hour['doy'] = df_hour.index.strftime('%j')
    df_hour['hour'] = df_hour.index.strftime('%H')
    df_hour['minute'] = df_hour.index.strftime('%M')
    df_phy = df_hour[['doy','hour','minute','latitude','longitude']]
    df_phy['longitude'] = df_phy.longitude * -1
    df_phy = df_phy.round({'latitude':3,'longitude':3})
    outfile = str(df_hour.argosid[0]) + '_for_phyllis.csv'
    df_phy.to_csv(outfile, sep=" ", index=False)
    