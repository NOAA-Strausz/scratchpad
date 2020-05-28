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
parser.add_argument('-t', '--temp', action="store_true", 
                    help="make plot of surface temperature")
parser.add_argument('-s', '--strain', action="store_true",
                    help="make plot of drogue strain")
parser.add_argument('-d', '--date', action="store_true",
                    help="add occasional date to track")
parser.add_argument('-c', '--cut', nargs=2, 
                    type=lambda x: datetime.strptime(x, "%Y-%m-%dT%H:%M:%S"),
                    help="date span in format '2019-01-01T00:00:00 2019-01-01T00:00:01' for example")
args=parser.parse_args()




filename=args.infile


df = pd.read_csv(filename)

df['datetime'] = pd.to_datetime(df['year_doy_hhmm'], format='%Y-%m-%d %H:%M:%S')

df.set_index(['datetime'], inplace=True)
df.drop(columns=['year_doy_hhmm','year_doy_hhmm.1'], inplace=True)
df['longitude'] = df.longitude * -1
#let's add columns for calculating speed
#first add next_lat and next_lon
df['next_lat'] = df.latitude.shift(-1)
df['next_lon'] = df.longitude.shift(-1)
#then calculate distance between points with the haversine function
df['dist'] = df.apply(lambda x: haversine((x.latitude, x.longitude), (x.next_lat, x.next_lon)), axis=1)
#next shift up the 'dist' column
df.dist.shift()
#now calculate the time difference
df['time'] = df.index
df['time2'] = df.time.shift(-1)
#make the time_delta
df['time_delta'] = df.time2 - df.time
#now make new column of seconds
df['seconds'] = df.time_delta.dt.total_seconds()
#now calculate speed in m/s
df['speed'] = df.dist * 1000 / df.seconds

proj = ccrs.LambertConformal(central_longitude=-165, central_latitude=60)
ax = plt.axes(projection=proj)
ax.natural_earth_shp(name='land', resolution='50m' )
ax.coastlines(resolution='50m')

#auto find extents
nlat = df.latitude.max() + 5
slat = df.latitude.min() - 5
wlon = df.longitude.max() - 20
elon = df.longitude.min() + 20

extents = [wlon, elon, slat, nlat]
ax.set_extent(extents)

#t = ax.scatter(df['longitude'], df['latitude'], s=5, c=df.sst, transform=ccrs.PlateCarree(), 
#               cmap=cmocean.cm.thermal, vmin=-2, vmax=10 )
#plt.colorbar(t)
plt.title(filename)
plt.show()

#now make new data frame that only has one point per day

df_day = df[df.index.hour == 0]
df_week = df_day[df_day.index.dayofweek == 0]
#adds in first and last points 
df_week = pd.concat([df.head(1),df_week,df.tail(1)])

df_week['date'] = df_week.index.strftime('%m/%d')
df_week = df_week[['longitude','latitude','date']]
ax.scatter(df_week['longitude'], df_week['latitude'], marker='+', s=100, transform=ccrs.PlateCarree())

#let's add columns for calculating speed
#only on a daily basis, so use df_day
#reduce amount of times, only one 00 hour per day
df_day['time'] = df_day.index
df_day['prev_time'] = df_day.time.shift()
df_day.drop(df_day[df_day.time.dt.day == df_day.prev_time.dt.day].index, inplace=True)

#first add next_lat and next_lon
df_day['next_lat'] = df_day.latitude.shift(-1)
df_day['next_lon'] = df_day.longitude.shift(-1)
#then calculate distance between points with the haversine function
df_day['dist'] = df_day.apply(lambda x: haversine((x.latitude, x.longitude), (x.next_lat, x.next_lon)), axis=1)
#next shift up the 'dist' column
df_day.dist.shift()
#now calculate the time difference

df_day['time2'] = df_day.time.shift(-1)
#make the time_delta
df_day['time_delta'] = df_day.time2 - df_day.time
#now make new column of seconds
df_day['seconds'] = df_day.time_delta.dt.total_seconds()
#now calculate speed in m/s
df_day['speed'] = df_day.dist * 100000 / df_day.seconds

speed = ax.scatter(df_day['longitude'], df_day['latitude'], s=100, c=df_day.speed, transform=ccrs.PlateCarree(), 
               cmap=cmocean.cm.speed, vmin=0, vmax=120 )
plt.colorbar(speed)
#this will keep only points where the hour is 00

#now let's keep only every 7th day  

#note that annotate requires the _as_mpl_transfomr(ax) for some unknown reason
#df_week.apply(lambda x: ax.annotate(x.date, xy=(x.longitude,x.latitude), xycoords=ccrs.PlateCarree()._as_mpl_transform(ax)), axis=1)
#here's another way to do it using plt.text
df_week.apply(lambda x: ax.text(x.longitude,x.latitude,'  '+x.date, transform=ccrs.PlateCarree()), axis=1)

df_hour=df.resample('H').mean()


#ax.scatter(df_hour['longitude'], df_hour['latitude'], s=5, c=df_hour.sst, transform=ccrs.PlateCarree(), 
#               cmap=cmocean.cm.thermal, vmin=-2, vmax=10 )

df_day2=df.resample('D').mean()

#ax.scatter(df_day2['longitude'], df_day2['latitude'], s=100, c=df_day2.sst, transform=ccrs.PlateCarree(), 
#               cmap=cmocean.cm.thermal, vmin=-2, vmax=10 )

#first add next_lat and next_lon
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

ax.scatter(df_hour['longitude'], df_hour['latitude'], s=10, c=df_hour.speed, transform=ccrs.PlateCarree(), 
               cmap=cmocean.cm.speed, vmin=0, vmax=120 )


#ax.text(-165, 65, 'boogers', transform=ccrs.PlateCarree())