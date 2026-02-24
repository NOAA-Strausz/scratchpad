#!/usr/bin/env python3
import xarray as xr
import numpy as np
import pandas as pd
import os

# --- CONFIGURATION ---
SEARCH_DIR = 'data/all_000' 
FILE_PATTERN = '.nc'

# Add or remove variables you want to extract right here!
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
    lon_grid_norm = (lon_grid + 180) % 360 - 180
    target_lon_norm = (target_lon + 180) % 360 - 180
    cos_factor = np.cos(np.deg2rad(target_lat))
    sq_diff = (lat_grid - target_lat)**2 + ((lon_grid_norm - target_lon_norm) * cos_factor)**2
    return np.unravel_index(sq_diff.argmin(), lat_grid.shape)

def main():
    print(f"Searching for {FILE_PATTERN} files in '{SEARCH_DIR}'...")
    files = get_files_recursively(SEARCH_DIR, FILE_PATTERN)
    
    if not files:
        print("No files found!")
        return

    print(f"Found {len(files)} files. Starting processing...\n")
    
    results_list = []
    cached_indices = {}
    
    for filepath in files:
        filename = os.path.basename(filepath)
        
        try:
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

                # --- Extract Data for Every Location ---
                for station, (j, i) in cached_indices.items():
                    row_data = {
                        'Date': formatted_date,
                        'Location': station,
                        'Source_File': filename
                    }
                    
                    for var_name in TARGET_VARIABLES:
                        if var_name in ds:
                            try:
                                if 'time' in ds[var_name].dims:
                                    val = ds[var_name].isel({ds[var_name].dims[0]: 0, ds[var_name].dims[-2]: j, ds[var_name].dims[-1]: i}).item()
                                else:
                                    val = ds[var_name].isel({ds[var_name].dims[-2]: j, ds[var_name].dims[-1]: i}).item()
                                row_data[var_name] = val
                            except Exception:
                                row_data[var_name] = np.nan
                        else:
                            row_data[var_name] = np.nan
                    
                    results_list.append(row_data)

        except Exception as e:
            print(f"  Error processing {filename}: {e}")

    # --- Save Results to Separate Files ---
    if results_list:
        print("\nSaving data to files...")
        df = pd.DataFrame(results_list)
        
        # Reorder columns to put Date, Location, File first, then variables
        cols = ['Date', 'Location'] + TARGET_VARIABLES + ['Source_File']
        cols = [c for c in cols if c in df.columns] 
        df = df[cols]
        
        # Sort globally by Date first
        df = df.sort_values(by=['Date'])
        
        # Loop through each unique location and save its own file
        for station in df['Location'].unique():
            # Filter dataframe for just this station
            station_df = df[df['Location'] == station]
            
            # Create filename (e.g., KS08_riops.csv)
            output_filename = f"{station}_riops.csv"
            
            # Save to CSV
            station_df.to_csv(output_filename, index=False)
            print(f"  -> Saved {len(station_df)} rows to {output_filename}")
            
        print("\nProcessing complete!")
    else:
        print("No data extracted.")

if __name__ == "__main__":
    main()
