#!/usr/bin/env python3

#masam2 nearest neighbor
#7/30/2025 this is an attempt to use this dataset:
#https://nsidc.org/data/g10005/versions/2
#program will find nearest data point to given lat and lon

import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math
import glob

dist=20 

#points of desired data

kotz = {'KS08':[66.935175,-163.824225], 'KS04':[66.61012917,-163.6313292],
                'KS01':[67.0714625,-163.7723917], 'KS02':[66.78301667,-163.7495125],
                'KS03':[66.683175,-164.4552458], 'KS07':[66.23546667,-162.2614]}

def find_box(lat1, lon1, dist, nm):
        
        #formula pulled from this website:
        #http://www.movable-type.co.uk/scripts/latlong.html
        
        if nm:
            r=3440
        else:
            r=6371
        
        dist = dist/2
        
        wlon = math.radians(lon1) + math.atan2(math.sin(math.radians(270)) *
                            math.sin(dist/r) * math.cos(math.radians(lat1)),
                            math.cos(dist/r) - math.sin(math.radians(lat1)) *
                            math.sin(math.radians(lat1)))
        elon = math.radians(lon1) + math.atan2(math.sin(math.radians(90)) *
                            math.sin(dist/r) * math.cos(math.radians(lat1)),
                            math.cos(dist/r) - math.sin(math.radians(lat1)) *
                            math.sin(math.radians(lat1)))
        nlat = math.asin(math.sin(math.radians(lat1)) * math.cos(dist/r) + 
                         math.cos(math.radians(lat1)) * math.sin(dist/r) * 
                         math.cos(math.radians(0)))
        slat = math.asin(math.sin(math.radians(lat1)) * math.cos(dist/r) + 
                         math.cos(math.radians(lat1)) * math.sin(dist/r) * 
                         math.cos(math.radians(180)))
        
        wlon = round(math.degrees(wlon), 4)
        elon = round(math.degrees(elon), 4)
        nlat = round(math.degrees(nlat), 4)
        slat = round(math.degrees(slat), 4)
        
        return([nlat,slat,wlon,elon])

def make_plot(lats, lons, ice, tlat, tlon):

    lat = lats
    lon = lons
    ice = ice
    
    # Target location in decimal degrees
    target_lat = tlat
    target_lon = tlon
    
    # Find nearest grid point
    lat_diff = abs(lat - target_lat)
    lon_diff = abs(lon - target_lon)
    total_diff = lat_diff + lon_diff
    min_index = total_diff.argmin(dim=["y", "x"])
    y_index = int(min_index["y"].values)
    x_index = int(min_index["x"].values)
    
    # Extract value and location
    matched_lat = float(lat[y_index, x_index].values)
    matched_lon = float(lon[y_index, x_index].values)
    ice_val = float(ice[14, y_index, x_index].values)  # Jan 15 is index 14
    
    print(f"Nearest value: {ice_val:.2f} at lat={matched_lat:.3f}, lon={matched_lon:.3f}")
    
    # Extract a 5x5 neighborhood around the point
    window = 2
    ys = slice(max(0, y_index - window), y_index + window + 1)
    xs = slice(max(0, x_index - window), x_index + window + 1)
    
    ice_subset = ice[14, ys, xs]
    lat_subset = lat[ys, xs]
    lon_subset = lon[ys, xs]
    
    # Flatten for scatter plotting
    lat_vals = lat_subset.values.flatten()
    lon_vals = lon_subset.values.flatten()
    ice_vals = ice_subset.values.flatten()
    
    # Plot
    plt.figure(figsize=(8, 6))
    sc = plt.scatter(lon_vals, lat_vals, c=ice_vals, cmap="Blues", edgecolors="black", s=100)
    plt.scatter(matched_lon, matched_lat, color="red", marker="x", s=120, label="Target Point")
    plt.colorbar(sc, label="Sea Ice Concentration")
    plt.title("Sea Ice Concentration Near Target Location\n(January 15, 2025)")
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    

# Load dataset
years = list(range(2012,2026))

files = []

for i in years:
    year = str(i)
    path = year + '/'
    files = files + glob.glob(path + '*.nc')
    files = sorted(files)

data=[]

for key, value in kotz.items():
    name = key
    target_lat = value[0]
    target_lon = value[1]
    filename = "masam2_"+name+"_"+str(target_lat)+"_"+str(target_lon)+"_"+str(dist)+"km.csv"
    print("Processing for "+name)

    for i in files:
    
        
    
        #file_path = "masam2_minconc40_202501_v2.nc" 
        ds = xr.open_dataset(i)
        ice = ds["sea_ice_concentration"]
        time = ds["time"]
        lats = ds["latitude"]
        lons = ds["longitude"]
        #make_plot(lats, lons, ice, target_lat, target_lon)
              
        #make pandas dataframe of just lats lons and ice %
        
        flat_lats = ds.latitude.values.flatten()
        flat_lons = ds.longitude.values.flatten()
        
        
        for t in range(len(ds.time)):
            date = pd.to_datetime(ds.time[t].values)
        
            flat_ice = ds.sea_ice_concentration.isel(time=t).values.flatten()
            
            build_dataframe = {'latitude': flat_lats, 'longitude': flat_lons,
                'sea_ice_concentration': flat_ice}
                
            df_ice = pd.DataFrame(build_dataframe)
            
                
            nlat, slat, wlon, elon = find_box(target_lat, target_lon, dist, nm='')
            
            df_ice_chopped = df_ice[(df_ice.latitude <= nlat) & (df_ice.latitude >= slat) & 
                                            (df_ice.longitude >= wlon) & (df_ice.longitude <= elon)]
            
            print("Ice concentration for "+name+" on "+date.strftime("%Y-%m-%d")+" is "+str(df_ice_chopped.sea_ice_concentration.mean()))
           
            data.append({"date": date, "sea_ice_concentration": df_ice_chopped.sea_ice_concentration.mean()})
        
    df_seaice = pd.DataFrame(data)
    df_seaice.set_index('date',inplace=True)
        
    df_seaice.to_csv(filename)
    data=[]
    
