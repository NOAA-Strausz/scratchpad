#!/usr/bin/env python3
import xarray as xr
import numpy as np
import pandas as pd
import os
import sys

# --- CONFIGURATION ---
SEARCH_DIR = 'data/all_000' 
OUTPUT_FILE = 'kotzebue_hs_timeseries.csv'
FILE_PATTERN = '.nc'

#target variables

TARGET_VARIABLES = ['hs', 'hi', 'strength', 'sigP']

# Locations
KOTZ_LOCATIONS = {
    'KS08': [66.935175, -163.824225],
    'KS04': [66.61012917, -163.6313292],
    'KS01': [67.0714625, -163.7723917],
    'KS02': [66.78301667, -163.7495125],
    'KS03': [66.683175, -164.4552458],
    'KS07': [66.23546667, -162.2614]
}

# --- HELPER FUNCTIONS ---
def get_files_recursively(directory, pattern):
    """Finds all files matching pattern in directory and subdirectories."""
    matches = []
    for root, dirnames, filenames in os.walk(directory):
        for filename in filenames:
            if filename.endswith(pattern):
                matches.append(os.path.join(root, filename))
    return sorted(matches)

def find_nearest_indices(lat_grid, lon_grid, target_lat, target_lon):
    """Finds (j, i) indices for nearest grid point."""
    # Normalize to [-180, 180]
    lon_grid_norm = (lon_grid + 180) % 360 - 180
    target_lon_norm = (target_lon + 180) % 360 - 180
    
    # Cosine correction
    cos_factor = np.cos(np.deg2rad(target_lat))
    sq_diff = (lat_grid - target_lat)**2 + ((lon_grid_norm - target_lon_norm) * cos_factor)**2
    
    return np.unravel_index(sq_diff.argmin(), lat_grid.shape)

# --- MAIN EXECUTION ---
print(f"Searching for {FILE_PATTERN} files in '{SEARCH_DIR}'...")
files = get_files_recursively(SEARCH_DIR, FILE_PATTERN)

if not files:
    print("No files found!")
else:
    print(f"Found {len(files)} files. Starting processing...")
    
    results_list = []
    cached_indices = {}
    
    for filepath in files:
        filename = os.path.basename(filepath)
        
        try:
            # We use 'with' to ensure the file closes, but for interactive work
            # you might sometimes want to keep 'ds' open. 
            # In this loop, we close it to save memory.
            with xr.open_dataset(filepath) as ds:
                
                # --- EXTRACT DATE FROM TIME VARIABLE ---
                if 'time' in ds and ds.time.size > 0:
                    raw_time = ds.time.values[0]
                    try:
                        formatted_date = pd.to_datetime(raw_time).strftime('%Y-%m-%d')
                    except:
                        formatted_date = str(raw_time).split()[0]
                else:
                    formatted_date = "Unknown"
                    
                print(f"Processing: {filename} -> Date: {formatted_date}")

                # --- Grid Setup (Run only once) ---
                if not cached_indices:
                    print("  > Calculating grid indices (first run)...")
                    
                    possible_lats = ['lat', 'latitude', 'nav_lat', 'TLAT']
                    possible_lons = ['lon', 'longitude', 'nav_lon', 'TLON']
                    
                    lat_name = next((n for n in possible_lats if n in ds), None)
                    lon_name = next((n for n in possible_lons if n in ds), None)
                    
                    if not lat_name or not lon_name:
                        print(f"  Skipping {filename}: Coordinates not found.")
                        continue
                        
                    lat_grid = ds[lat_name].values
                    lon_grid = ds[lon_name].values
                    
                    for station, coords in KOTZ_LOCATIONS.items():
                        j, i = find_nearest_indices(lat_grid, lon_grid, coords[0], coords[1])
                        cached_indices[station] = (j, i)

                # --- Extract Data ---
                for station, (j, i) in cached_indices.items():
                    try:
                        # Handle dimensions
                        #work through target variabls
                        for var_name in TARGET_VARIABLES:


                        if 'time' in ds['hs'].dims:
                            val = ds['hs'].isel({ds['hs'].dims[0]: 0, ds['hs'].dims[-2]: j, ds['hs'].dims[-1]: i}).item()
                        else:
                            val = ds['hs'].isel({ds['hs'].dims[-2]: j, ds['hs'].dims[-1]: i}).item()
                        
                        results_list.append({
                            'Date': formatted_date,
                            'Location': station,
                            'hs_value': val,
                            'Source_File': filename
                        })
                    except KeyError:
                        print(f"  Warning: 'hs' variable missing in {filename}")
                        
        except Exception as e:
            print(f"  Error processing {filename}: {e}")

    # --- Save Results ---
    if results_list:
        df = pd.DataFrame(results_list)
        df = df.sort_values(by=['Date', 'Location'])
        
        df.to_csv(OUTPUT_FILE, index=False)
        print(f"\nProcessing complete!")
        print(f"Data saved to: {OUTPUT_FILE}")
        print(df.head())
    else:
        print("No data extracted.")
