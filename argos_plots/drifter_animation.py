#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 25 15:22:45 2019

@author: strausz
"""

import pandas as pd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import matplotlib.animation
#import cartopy.feature as cfeature
import cmocean

filename='139914_2018.csv'



df = pd.read_csv(filename)

df['datetime'] = pd.to_datetime(df['year_doy_hhmm'], format='%Y-%m-%d %H:%M:%S')

df.set_index(['datetime'], inplace=True)
df.drop(columns=['year_doy_hhmm','year_doy_hhmm.1'], inplace=True)
df['longitude'] = df.longitude * -1

df = df.resample('W').mean()

proj = ccrs.LambertConformal(central_longitude=-165, central_latitude=60)
fig = plt.figure()
ax = plt.axes(projection=proj)
ax.natural_earth_shp(name='land', resolution='50m' )
ax.coastlines(resolution='50m')
ax.set_extent([-150, -170, 70, 75])

def make_plot(x):
    
    ax.scatter(x['longitude'], x['latitude'], s=5, c=x.sst, transform=ccrs.PlateCarree(), 
                  cmap=cmocean.cm.thermal, vmin=-2, vmax=10 )
    
    
#for index, row in df.iterrows():
#    make_plot(row)