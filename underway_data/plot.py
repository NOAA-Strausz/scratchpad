#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 24 14:45:54 2020

@author: strausz
"""

import pandas as pd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cmocean

filename = '1min_time_fixed_gps.csv'

df = pd.read_csv(filename, parse_dates=[0])

proj = ccrs.LambertConformal(central_longitude=-165, central_latitude=60)
fig = plt.figure(figsize=(10,10))
ax = plt.axes(projection=proj)
ax.natural_earth_shp(name='land', resolution='50m' )
ax.coastlines(resolution='50m')
ax.set_title("OS1901 SST")
wlon = 180
elon = -150
slat = 51
nlat = 73
extents = [wlon, elon, slat, nlat]
ax.set_extent(extents)
#
#chl = ax.scatter(df['longitude'], df['latitude'], s=10, c=df.chl_conc, transform=ccrs.PlateCarree(), 
#               cmap=cmocean.cm.algae, vmin=0, vmax=15 )
#plt.colorbar(chl)

hull_temp = ax.scatter(df['longitude'], df['latitude'], s=10, c=df.hull_temp, transform=ccrs.PlateCarree(), 
               cmap=cmocean.cm.thermal, vmin=0, vmax=15 )

plt.colorbar(hull_temp)

fig.savefig('sst_map.png')