#!/usr/bin/env python3

import xarray as xr
import numpy as np
import pandas as pd
import sys

def find_nearest_indices(lat_grid, lon_grid, target_lat, target_lon):
    """
    Finds the indices (j, i) of the grid point nearest to the target lat/lon.
    Handles 0-360 vs -180/180 longitude mismatches.
    """
    # 1. Normalize both grids to [-180, 180] range to ensure they match
    # This handles cases where model is 0-360 and target is negative (West)
    lon_grid_norm = (lon_grid + 180) % 360 - 180
    target_lon_norm = (target_lon + 180) % 360 - 180
    
    # 2. Calculate squared Euclidean distance with cosine correction
    cos_factor = np.cos(np.deg2rad(target_lat))
    sq_diff = (lat_grid - target_lat)**2 + ((lon_grid_norm - target_lon_norm) * cos_factor)**2
    
    # 3. Find index of minimum distance
    min_idx = sq_diff.argmin()
    
    return np.unravel_index(min_idx, lat_grid.shape)

def main():
    # --- CONFIGURATION ---
    # Update this to your full filename if needed
    filename = 'data/2025/2025010100_000_iau.nc' 
    output_file = 'kotzebue_hs_data.csv'
    
    # --- LOAD DATA ---
    try:
        # Load dataset
        # decode_times=False is sometimes safer for non-standard calendars, 
        # but usually True is fine.
        ds = xr.open_dataset(filename)
    except FileNotFoundError:
        print(f"Error: Could not find file '{filename}'.")
        return
    except Exception as e:
        print(f"Error opening file: {e}")
        return

    # --- DETECT COORDINATES ---
    # Added TLAT/TLON to the list
    possible_lats = ['lat', 'latitude', 'nav_lat', 'TLAT']
    possible_lons = ['lon', 'longitude', 'nav_lon', 'TLON']
    
    # Check both coordinates and data variables
    lat_name = next((name for name in possible_lats if name in ds), None)
    lon_name = next((name for name in possible_lons if name in ds), None)

    if not lat_name or not lon_name:
        print(f"Error: Could not find latitude/longitude coordinates.")
        print(f"Variables found: {list(ds.keys())}")
        return

    print(f"Using coordinates: {lat_name}, {lon_name}")
    lat_grid = ds[lat_name].values
    lon_grid = ds[lon_name].values

    # --- DEFINE LOCATIONS ---
    kotz = {
        'KS08': [66.935175, -163.824225],
        'KS04': [66.61012917, -163.6313292],
        'KS01': [67.0714625, -163.7723917],
        'KS02': [66.78301667, -163.7495125],
        'KS03': [66.683175, -164.4552458],
        'KS07': [66.23546667, -162.2614]
    }

    results_list = []
    print("\nExtracting data...")
    print(f"{'Station':<8} | {'Model Lat':<10} | {'Model Lon':<10} | {'hs Value':<10}")
    print("-" * 50)

    # --- PROCESS LOOP ---
    for station, coords in kotz.items():
        target_lat, target_lon = coords
        
        # Find nearest grid point
        j_idx, i_idx = find_nearest_indices(lat_grid, lon_grid, target_lat, target_lon)
        
        try:
            # Extract 'hs' value
            # Handle potential dimensions: (time, y, x) or (y, x)
            if 'time' in ds['hs'].dims:
                # Assuming time is the first dimension, taking the first time step
                 data_slice = ds['hs'].isel({ds['hs'].dims[0]: 0, 
                                             ds['hs'].dims[-2]: j_idx, 
                                             ds['hs'].dims[-1]: i_idx})
            else:
                 data_slice = ds['hs'].isel({ds['hs'].dims[-2]: j_idx, 
                                             ds['hs'].dims[-1]: i_idx})
            
            val = data_slice.values.item()
            
            # Get the actual model coordinates for verification
            model_lat_pt = lat_grid[j_idx, i_idx]
            model_lon_pt = lon_grid[j_idx, i_idx]
            
            # Print to screen
            print(f"{station:<8} | {model_lat_pt:<10.4f} | {model_lon_pt:<10.4f} | {val:.6f}")

            # Append to results
            results_list.append({
                'Station_ID': station,
                'Target_Lat': target_lat,
                'Target_Lon': target_lon,
                'Model_Lat': model_lat_pt,
                'Model_Lon': model_lon_pt,
                'hs_value': val
            })
            
        except Exception as e:
            print(f"Error extracting {station}: {e}")

    # --- SAVE TO CSV ---
    if results_list:
        df = pd.DataFrame(results_list)
        df.to_csv(output_file, index=False)
        print(f"\nSuccess! Data saved to {output_file}")
    else:
        print("\nNo data extracted.")

if __name__ == "__main__":
    main()
