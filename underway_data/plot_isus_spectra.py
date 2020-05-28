#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 14:40:04 2020

@author: strausz
"""

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates





spectra_file='isus_204_spectra.csv'
isus_file='raw_combined/os1901_isus_full_df_only.csv'


#pull spectra into list from file 
spectra = [line.rstrip() for line in open(spectra_file)]

df = pd.read_csv(isus_file, parse_dates=[0], header=None)

df.set_index(df[0], inplace=True)

dfs=df[df.columns[21:277]]

dfs.columns=spectra

dfs_3m=dfs.resample('3T').mean()

fig, ax = plt.subplots(figsize=(9,4))

image = ax.pcolormesh(dfs_3m.index, dfs_3m.columns, dfs_3m.T, cmap=plt.cm.plasma, vmin=1000, vmax=22000)

plt.ylabel('Y(nm)')
myfmt = mdates.DateFormatter('%m-%Y')
ax.xaxis.set_major_formatter(myfmt)
ax.xaxis.set_major_locator(plt.MaxNLocator(6))
ax.yaxis.set_major_locator(plt.MaxNLocator(10))
plt.title("OS1901 Isus 204 Spectra")

fig.colorbar(image)
