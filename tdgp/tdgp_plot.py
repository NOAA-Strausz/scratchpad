#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 11:29:02 2020

@author: strausz
"""

import pandas as pd
import matplotlib.pyplot as plt


filename = 'mini_tdgp_seacat_test.csv'

dtypes={1:str, 2:str, 3:str, 4:str, 5:str, 6:str}

df = pd.read_csv(filename, sep=' |,',dtype=dtypes,header=None)

df['date']=df[1]+df[2]+df[3]+df[4]+df[5]+df[6]
df['datetime'] = pd.to_datetime(df['date'], format='%Y%m%d%H%M%S')

df.set_index(['datetime'], inplace=True)

df.drop(range(7),axis=1,inplace=True)
df.drop(range(10,14), axis=1,inplace=True)
df.drop(['date'],axis=1,inplace=True)


df_hour = df.resample('H').mean()
df_hour.dropna(inplace=True)
df_hour.columns=['DGP', 'Temp', 'Volts']

fig = plt.figure()
ax = plt.axes()
ax.plot(df_hour.index, df_hour.DGP)