import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from pyproj import Proj, Transformer

# Load the dataset
file_path = "NSIDC0081_SEAICE_PS_N25km_20250506_v2.0.nc"
ds = xr.open_dataset(file_path)

# Load ice data and x/y coordinates
ice_data = ds["F17_ICECON"].squeeze()
x = ds["x"].values
y = ds["y"].values
xx, yy = np.meshgrid(x, y)

# Convert projection coordinates to lat/lon
proj_stere = Proj(proj='stere', lat_ts=70, lat_0=90, lon_0=-45, x_0=0, y_0=0, ellps='WGS84')
transformer = Transformer.from_proj(proj_stere, "epsg:4326", always_xy=True)
lon, lat = transformer.transform(xx, yy)

# Mask land/missing data (values >= 251 are invalid, per NSIDC convention)
ice_data_masked = np.ma.masked_where((ice_data >= 251) | (ice_data > 100), ice_data)

# Normalize values to 0–1 range for better color interpretation (optional)
ice_data_norm = ice_data_masked / 100.0  # values between 0.0 and 1.0

# Create plot
fig, ax = plt.subplots(figsize=(10, 8), subplot_kw={'projection': ccrs.PlateCarree()})
ax.set_title("Sea Ice Concentration – Bering Sea Region (F17 Sensor)\n2025-05-06", fontsize=14)

# Add land and coastlines
ax.add_feature(cfeature.LAND, color='lightgray')
ax.add_feature(cfeature.COASTLINE)

# Set lat/lon bounds for Bering Sea
ax.set_extent([-180, -140, 50, 72], crs=ccrs.PlateCarree())

# Plot the data using a reversed colormap (white = ice, blue = open water)
ice_plot = ax.pcolormesh(lon, lat, ice_data_norm,
                         transform=ccrs.PlateCarree(),
                         cmap='Blues_r', shading='auto')

# Add colorbar with reversed label scale
cbar = plt.colorbar(ice_plot, ax=ax, orientation='vertical', shrink=0.6)
cbar.set_label('Sea Ice Concentration (fraction)', fontsize=11)

plt.tight_layout()
plt.show()

