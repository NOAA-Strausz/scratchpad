#!/usr/bin/env python3

import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import sys

# The name of the NetCDF file
file_name = '2025/masam2_minconc40_202507_v2.nc'

# Open the NetCDF file using xarray
try:
    ds = xr.open_dataset(file_name)
    print("Dataset loaded successfully.")
except FileNotFoundError:
    print(f"Error: The file '{file_name}' was not found. Please make sure it is in the same directory as this script.")
    sys.exit()

#subset the data
min_lat, max_lat = 66, 67.5
min_lon, max_lon = -165, -160

mask = (ds['latitude'] >= min_lat) & \
       (ds['latitude'] <= max_lat) & \
       (ds['longitude'] >= min_lon) & \
       (ds['longitude'] <= max_lon)

subset = ds.where(mask, drop=True)

subset_ice_only_first_day = subset.sea_ice_concentration.isel(time=0)



# Create the plot figure and axes
plt.figure(figsize=(10, 8))

# --- NEW: Set the map projection for the axes using Cartopy's Plate Carree projection ---
ax = plt.axes(projection=ccrs.PlateCarree())

# Store the plot object in a variable called 'mappable_plot'
# --- NEW: Tell xarray to use the Plate Carree projection for the data ---
mappable_plot = subset_ice_only_first_day.plot(cmap='Blues_r', transform=ccrs.PlateCarree())

# Set the plot limits to zoom in on Kotzebue Sound
# Use the axes object 'ax' instead of plt
ax.set_extent([-165, -160, 66, 67.5], crs=ccrs.PlateCarree())

# Add coastline features to the map
# --- NEW: Add coastline features using Cartopy's built-in features ---
ax.add_feature(cfeature.COASTLINE, linewidth=0.5)
ax.add_feature(cfeature.BORDERS, linestyle=':')
ax.add_feature(cfeature.LAND, facecolor='lightgray')

# Add titles and labels for clarity
ax.set_title('MASAM2 Sea Ice Concentration (Kotzebue Sound)')
ax.set_xlabel('Longitude')
ax.set_ylabel('Latitude')

# Pass the mappable_plot to the colorbar function
plt.colorbar(mappable_plot, label='Sea Ice Concentration (0-1)')

# Save the plot to a file
output_filename = 'kotzebue_sound_sea_ice_with_coast.png'
plt.savefig(output_filename, dpi=300)
print(f"Plot saved to '{output_filename}'")

# Close the dataset
ds.close()
