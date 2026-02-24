import xarray as xr
import matplotlib.pyplot as plt

# The name of the NetCDF file
file_name = 'masam2_minconc40_202507_v2.nc'

# Open the NetCDF file using xarray
try:
    ds = xr.open_dataset(file_name)
    print("Dataset loaded successfully.")
except FileNotFoundError:
    print(f"Error: The file '{file_name}' was not found. Please make sure it is in the same directory as this script.")
    exit()

# Select the sea ice concentration variable
sea_ice_concentration = ds['sea_ice_concentration']

# Select a single time step to plot
first_day_concentration = sea_ice_concentration.isel(time=0)

# Create the plot
plt.figure(figsize=(10, 8))

# Store the plot object in a variable called 'mappable_plot'
mappable_plot = first_day_concentration.plot(cmap='Blues_r')

# Set the plot limits to zoom in on Kotzebue Sound
plt.xlim(-165, -160)
plt.ylim(66, 67.5)

# Add titles and labels for clarity
plt.title('MASAM2 Sea Ice Concentration (Kotzebue Sound)')
plt.xlabel('Longitude')
plt.ylabel('Latitude')

# Pass the mappable_plot to the colorbar function
plt.colorbar(mappable_plot, label='Sea Ice Concentration (0-1)')

# Save the plot to a file
output_filename = 'kotzebue_sound_sea_ice.png'
plt.savefig(output_filename, dpi=300)
print(f"Plot saved to '{output_filename}'")

# Close the dataset
ds.close()
