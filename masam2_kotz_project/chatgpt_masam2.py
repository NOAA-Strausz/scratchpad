import xarray as xr
import numpy as np
import matplotlib.pyplot as plt
#7/30/2025 this is an attempt to use this dataset:
#https://nsidc.org/data/g10005/versions/2
#program will find nearest data point to given lat and lon


# Load your NetCDF file
file_path = "masam2_minconc40_202501_v2.nc"  # Adjust path if needed
ds = xr.open_dataset(file_path)

# Extract coordinates and concentration
lat = ds["latitude"]
lon = ds["longitude"]
ice = ds["sea_ice_concentration"]

# Target location in decimal degrees
target_lat = 66 + 56.1250 / 60     # 66° 56.1250′ N
target_lon = -(163 + 49.3320 / 60) # 163° 49.3320′ W

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

import matplotlib.patches as patches

# Determine corner of the box (based on center point and grid spacing)
# Assumes uniform grid, 4 km per grid cell
# Convert to approximate degrees:
km_per_deg_lat = 111.0
km_per_deg_lon = 111.0 * np.cos(np.radians(matched_lat))

delta_lat = 4 / km_per_deg_lat  # ~0.036°
delta_lon = 4 / km_per_deg_lon  # depends on latitude

# Create the rectangle: lower-left corner, width, height
lower_left_lon = matched_lon - delta_lon / 2
lower_left_lat = matched_lat - delta_lat / 2

rect = patches.Rectangle(
    (lower_left_lon, lower_left_lat),
    delta_lon,
    delta_lat,
    linewidth=2,
    edgecolor='red',
    facecolor='none',
    label="4km x 4km Cell"
)

plt.gca().add_patch(rect)
plt.show()
