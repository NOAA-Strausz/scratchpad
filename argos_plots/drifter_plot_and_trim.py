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
import numpy as np



parser = argparse.ArgumentParser(description='Plot drifter track on map')
parser.add_argument('infile', type=str, help='full path to input file')
parser.add_argument('-p', '--plot', type=str, 
                    help="make plot of 'sst', 'strain', or 'speed' ")
parser.add_argument('-f', '--file', action="store_true",
                    help="output csv file of data")
parser.add_argument('-i', '--ice', action="store_true",
                    help="add ice concentration as last field and output file")
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


#the following is info needed for adding the ice concentration 
latfile='/home/makushin/strausz/ecofoci_github/EcoFOCI_ssmi_ice/psn25lats_v3.dat'
lonfile='/home/makushin/strausz/ecofoci_github/EcoFOCI_ssmi_ice/psn25lons_v3.dat'
#locations of ice files
bootstrap = '/home/akutan/strausz/ssmi_ice/data/bootstrap/'
nrt = '/home/akutan/strausz/ssmi_ice/data/nrt/'
#latest available bootstrap year will need to be changed as new data comes in
boot_year = 2018

filename=args.infile

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
        ice[ice == 110] = 0 #110 is land
        ice[ice == 120] = 100 #120 is polar hole
    else: 
        ice=np.nan
    
    icefile.close()
    
    return ice;

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
    nearest_ice = df_ice.loc[df_ice.dist.idxmin()].ice_conc
    return nearest_ice
    
def lon_360(lon):
    if lon < 0:
        return 360 + lon
    else:
        return lon

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
    elif var == 'ice_concentration':
        vmin, vmax, cmap = 0, 100, cmocean.cm.ice
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
    df=pd.read_csv(filename, skiprows=1, header=0, names=names, parse_dates=[3])
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
if args.ice:
    df_hour['lon_360'] = df_hour.apply(lambda x: lon_360(x.longitude), axis=1)
    df_hour['datetime'] = df_hour.index
    df_hour.dropna(inplace=True)
    df_hour['argosid']=df_hour.argosid.astype(int)
    df_hour['latitude']=df_hour.latitude.round(decimals=3)
    df_hour['longitude']=df_hour.longitude.round(decimals=3)
    df_hour['voltage']=df_hour.voltage.round(decimals=2)
    df_hour['sst']=df_hour.sst.round(decimals=2)
    df_hour['strain']=df_hour.strain.round(decimals=2)
    #now group by doy
    #add blank ice column
    #df_hour['ice_concentration'] = ''
    ice_conc = []
    groups = df_hour.groupby(df_hour.index.dayofyear)
    for name, group in df_hour.groupby(df_hour.index.dayofyear):
        #print(name)
        #print(group.latitude)
        date = group.iloc[0].datetime.strftime("%Y%m%d")
        if group.iloc[0].datetime.year <= boot_year:
            ice_file = bootstrap + str(group.iloc[0].datetime.year) + "/" + "bt_" + date + "_f17_v3.1_n.bin"
        else:
            ice_file = nrt + "nt_" + date + "_f18_nrt_n.bin"
        print("opening file: " + ice_file)
        wlon = group.lon_360.min() - .25
        elon = group.lon_360.max() + .25
        nlat = group.latitude.max() + .25
        slat = group.latitude.min() - .25
        data_ice={'latitude':decode_latlon(latfile), 'longitude':decode_latlon(lonfile),
          'ice_conc':decode_datafile(ice_file)}
        df_hour_ice=pd.DataFrame(data_ice)
        df_hour_ice.dropna(inplace=True)
    
        df_hour_ice['lon_360'] = df_hour_ice.apply(lambda x: lon_360(x.longitude), axis=1)
        df_hour_ice_chopped = df_hour_ice[(df_hour_ice.latitude < nlat) & (df_hour_ice.latitude > slat) & (df_hour_ice.lon_360 > wlon) & (df_hour_ice.lon_360 < elon)]
        ice_conc = ice_conc + group.apply(lambda x: get_ice(x, df_hour_ice_chopped), axis=1).to_list()
                
        
        #print("the ice concentration is: "+ice_concentration)
        #df_hour_ice_chopped['dist'] = df_hour_ice_chopped.apply(lambda x: haversine((data.latitude, data.longitude), (x.latitude, x.longitude)), axis=1)
    df_hour['ice_concentration'] = ice_conc
    df_hour=df_hour.drop(['lon_360', 'datetime'], axis=1)
    df_out = df_hour[['argosid','latitude','longitude','sst','strain','voltage','speed','ice_concentration']]
    df_out = df_out.round({'latitude':3, 'longitude':3,'sst':2,'strain':1,'voltage':1,'speed':1, 'ice_concentration':1})
    outfile = str(df_hour.argosid[0]) + "_with_ice.csv"
    df_out.to_csv(outfile)

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
    
