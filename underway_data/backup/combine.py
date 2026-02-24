#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 23 14:37:24 2020

@author: strausz
"""

import pandas as pd
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import cmocean


gps_file = 'os1901L1_gps_5min.csv'
eco_file = 'os1901L1_eco_5min.csv'
isus_file = 'os1901L1_isus_5min.csv'
sbe48_file = 'os1901L1_sbe48_5min.csv'
optode_file = 'os1901L1_optode_5min.csv'
tsg_file = 'os1901L1_tsg_5min.csv'
tdgp_file = 'os1901L1_tdgp_5min.csv'

df_gps = pd.read_csv(gps_file, parse_dates=[0])
df_eco = pd.read_csv(eco_file, parse_dates=[0])
df_isus = pd.read_csv(isus_file, parse_dates=[0])
df_sbe48 = pd.read_csv(sbe48_file, parse_dates=[0])
df_optode = pd.read_csv(optode_file, parse_dates=[0])
df_tsg = pd.read_csv(tsg_file, parse_dates=[0])
df_tdgp = pd.read_csv(tdgp_file, parse_dates=[0])

df_gps.set_index(['timestamp'], inplace=True)
df_eco.set_index(['timestamp'], inplace=True)
df_isus.set_index(['timestamp'], inplace=True)
df_sbe48.set_index(['timestamp'], inplace=True)
df_optode.set_index(['timestamp'], inplace=True)
df_tsg.set_index(['timestamp'], inplace=True)
df_tdgp.set_index(['timestamp'], inplace=True)

df_sbe48.rename(columns={'temp':'hull_temp'},inplace=True)
df_optode.rename(columns={'temp':'O2_temp'},inplace=True)
df_tsg.rename(columns={'temp':'tsg_temp'},inplace=True)
df_tdgp.rename(columns={'temp':'tdgp_temp'},inplace=True)

df=pd.concat([df_gps,df_eco,df_isus,df_sbe48,df_optode,df_tsg,df_tdgp], axis=1, sort=False)

proj = ccrs.LambertConformal(central_longitude=-165, central_latitude=60)
ax = plt.axes(projection=proj)
ax.natural_earth_shp(name='land', resolution='50m' )
ax.coastlines(resolution='50m')
wlon = 180
elon = -150
slat = 50
nlat = 75
extents = [wlon, elon, slat, nlat]
ax.set_extent(extents)

chl = ax.scatter(df['longitude'], df['latitude'], s=10, c=df.hull_temp, transform=ccrs.PlateCarree(), 
               cmap=cmocean.cm.thermal, vmin=0, vmax=15 )
plt.colorbar(chl)