#!/usr/bin/env python3

import geopandas as gpd
from shapely.geometry import Point
import pandas as pd
from pyproj import CRS, Transformer
import glob
import re
from datetime import datetime
import os

# 1. Define File and CRS Information
# Assuming the shapefile component ARCTIC20250103.shp is in the same directory
SHAPEFILE_PATH = "ARCTIC20250103/ARCTIC20250103.shp" 

#get list of all files

years = list(range(2025,2026)) #adjust to decrease range

files = []
base_path = "/home/makushin/strausz/ecofoci_github/scratchpad/sigrid3/data/"
for i in years:
    year = str(i)
    path = base_path + year + '/'
    files = files + glob.glob(path + '*.zip')
    files = sorted(files)

# CRS from ARCTIC20250103.prj 
SHAPEFILE_WKT = 'PROJCS["WGS_1984_Stereographic_North_Pole",GEOGCS["GCS_WGS_1984",DATUM["D_WGS_1984",SPHEROID["WGS_1984",6378137.0,298.257223563]],PRIMEM["Greenwich",0.0],UNIT["Degree",0.0174532925199433]],PROJECTION["Stereographic_North_Pole"],PARAMETER["False_Easting",0.0],PARAMETER["False_Northing",0.0],PARAMETER["Central_Meridian",180.0],PARAMETER["Standard_Parallel_1",60.0],UNIT["Meter",1.0]]'


# Fields of interest (from ARCTIC20250103.shp.xml)
FIELDS_OF_INTEREST = ['ICECODE', 'CT', 'CA', 'CB', 'CC', 'SA', 'SB', 'SC', 'POLY_TYPE']

#station locations

kotz = {'KS08':[66.935175,-163.824225], 'KS04':[66.61012917,-163.6313292],
                'KS01':[67.0714625,-163.7723917], 'KS02':[66.78301667,-163.7495125],
                'KS03':[66.683175,-164.4552458], 'KS07':[66.23546667,-162.2614]}

mml_moorings = {
    'CH01': [75.1003383, -168.0005900],
    'BF02': [71.7514667, -154.4712500],
    'PB01': [71.2066800, -158.0140700],
    'WT01': [71.0458667, -160.5089000],
    'IC03': [71.8339500, -165.9033000],
    'IC01': [70.7983833, -163.0811167],
    'PH01': [67.9074500, -168.2026500],
    'NM01': [64.8580167, -168.4480333]
}

#one item for testing purposes, comment out above if used
#kotz = {'KS08':[66.935175,-163.824225], 'KS04':[66.61012917,-163.6313292]}

# 3. Load the Shapefile

def load_and_standardize(shp_path):
    """
    Loads a SIGRID-3 shapefile and standardizes its columns to match the 
    modern (post-2008/2026) format based on specific era variations.
    """
    print(f"DEBUG: Attempting to load {shp_path}...")
    
    if not os.path.exists(shp_path):
        print(f"ERROR: File not found at {shp_path}")
        return None

    try:
        # Load the file
        gdf = gpd.read_file(shp_path)
        print(f"DEBUG: Loaded {len(gdf)} rows. Original Columns: {list(gdf.columns)}")

        # --- 1. PROJECTION FIX ---
        if gdf.crs is None:
            print("DEBUG: CRS is missing. Forcing Arctic Stereographic.")
            custom_crs = "+proj=stere +lat_0=90 +lat_ts=60 +lon_0=180 +datum=WGS84 +units=m +no_defs"
            gdf.set_crs(custom_crs, inplace=True)
        else:
            print(f"DEBUG: CRS found: {gdf.crs.name}")

        # --- 2. COLUMN MAPPING ---
        rename_map = {}
        # Map old Ice Stage columns to new standard
        if 'CN' in gdf.columns: rename_map['CN'] = 'SO'
        if 'CD' in gdf.columns: rename_map['CD'] = 'SD'
        
        # Map old Form columns (CF) to new standard (FP)
        if 'CF' in gdf.columns and 'FP' not in gdf.columns:
            rename_map['CF'] = 'FP'
            
        # Map old Geometry columns
        if 'AREA' in gdf.columns: rename_map['AREA'] = 'Shape_Area'
        if 'PERIMETER' in gdf.columns: rename_map['PERIMETER'] = 'Shape_Leng'

        if rename_map:
            print(f"DEBUG: Renaming columns: {rename_map}")
            gdf.rename(columns=rename_map, inplace=True)

        # --- 3. FILL MISSING DATA ---
        if 'ICECODE' not in gdf.columns:
            gdf['ICECODE'] = "N/A (Old Format)"
        
        # Fill NaN values with empty strings to prevent errors
        string_cols = ['CT', 'CA', 'CB', 'CC', 'SA', 'SB', 'SC', 'FP', 'SO', 'SD', 'ICECODE']
        for col in string_cols:
            if col in gdf.columns:
                gdf[col] = gdf[col].fillna("")

        print("DEBUG: Standardization complete.")
        return gdf

    except Exception as e:
        print(f"CRITICAL ERROR in load_and_standardize: {e}")
        return None



# 4. Prepare for Reprojection
# Source CRS is WGS 84 (standard for Lat/Lon)
source_crs = CRS.from_epsg(4326) # EPSG:4326 is WGS 84


# 5. Process and Query Points
results = []
for name, coords in kotz.items():
    # Reproject (Lon, Lat) to (X, Y) in the target CRS
    filename_prefix = name
    lat, lon = coords
    # 6. Perform Point-in-Polygon Query
    # Use spatial join to find the polygon that contains the point
    # Create a temporary GeoSeries for the query point
    
    

    for file in files:
        #get date from filename
        pattern = r'(\d{8})\.'
        match = re.search(pattern, file)
        date_string = None
        if match:
        # 4. Extract the captured group (the raw date string: "20250103")
            date_string = match.group(1)
    
            # 5. Parse the date string and reformat it
            try:
            # datetime.strptime() parses the string from the format 'YYYYMMDD'
                date_object = datetime.strptime(date_string, '%Y%m%d')
        
                # strftime() formats the date object into the desired 'YYYY-MM-DD' output format
                file_date = date_object.strftime('%Y-%m-%d')
        
            except ValueError as e:
                print(f"Error: The extracted string '{date_string}' is not a valid date. {e}")
        else:
            print("Error: 8-digit date string not found in the filename.")

        print("Working on location " + name + " loading file for " + file_date)
        gdf_ice = load_and_standardize(file)
        if gdf_ice is not None and not gdf_ice.empty:
    
            # --- OPTIMIZATION START ---
            # Force the creation of a spatial index (R-tree) to speed up searching
            # This acts like an index in a book, letting pandas skip unrelated data.
            _ = gdf_ice.sindex 
            # --- OPTIMIZATION END ---
        
        
            target_crs = gdf_ice.crs
            transformer = Transformer.from_crs(source_crs, target_crs, always_xy=True)
            x, y = transformer.transform(lon, lat) 
            query_point_projected = Point(x, y)
            
            # 4. FAST SEARCH (Filter first, then check geometry)
            # Instead of checking all 10,000 polygons, we use the index to find 
            # only the polygons whose bounding box overlaps our point.
            possible_matches_index = list(gdf_ice.sindex.query(query_point_projected, predicate='contains'))
            possible_matches = gdf_ice.iloc[possible_matches_index]
        
            # 5. Precise Check (on the small subset only)
            # Now we do the heavy 'contains' math only on the candidates
            precise_matches = possible_matches[possible_matches.contains(query_point_projected)]
            
            if not precise_matches.empty:
                row = precise_matches.iloc[0]
                # ... process your results ...
                print(f"Found ice info: {row.get('ICECODE', 'N/A')}")
            
        
        # 6. Perform Point-in-Polygon Query
        # Use spatial join to find the polygon that contains the point
        # Create a temporary GeoSeries for the query point
        query_gdf = gpd.GeoDataFrame(geometry=[query_point_projected], crs=target_crs)
    
        # Spatial join: 'sjoin_within' finds features in gdf_ice that contain the points in query_gdf
        # Use 'how=right' to keep all original query points, even if they fall outside all polygons
        join_result = gpd.sjoin(query_gdf, gdf_ice[FIELDS_OF_INTEREST + ['geometry']], how="left", predicate='within')
        
        # 7. Extract Results
        if not join_result.empty and join_result.iloc[0]['ICECODE'] is not None:
            # Get the first (and should be only) matching row's attributes
            data = join_result.iloc[0].to_dict()
            
            # Format the output
            result_entry = {
                'Date': file_date,
                'Latitude': lat,
                'Longitude': lon,
                'ICECODE': data.get('ICECODE'),
                'Total Conc. (CT)': data.get('CT'),
                'Conc. 1st Thickest (CA)': data.get('CA'),
                'Conc. 2nd Thickest (CB)': data.get('CB'),
                'Conc. 3rd Thickest (CC)': data.get('CC'),
                'Stage 1st Thickest (SA)': data.get('SA'),
                'Stage 2nd Thickest (SB)': data.get('SB'),
                'Stage 3rd Thickest (SC)': data.get('SC'),
                'Polygon Type (POLY_TYPE)': data.get('POLY_TYPE')
            }
        else:
            # No polygon found at this location
            result_entry = {   
                'Date': file_date,             
                'Latitude': lat,
                'Longitude': lon,
                'ICECODE': 'N/A',
                'Total Conc. (CT)': 'N/A',
                'Conc. 1st Thickest (CA)': 'N/A',
                'Conc. 2nd Thickest (CB)': 'N/A',
                'Conc. 3rd Thickest (CC)': 'N/A',
                'Stage 1st Thickest (SA)': 'N/A',
                'Stage 2nd Thickest (SB)': 'N/A',
                'Stage 3rd Thickest (SC)': 'N/A',
                'Polygon Type (POLY_TYPE)': 'N/A'
            }
            
        results.append(result_entry)
    

    # 8. Report the Final Results
    #print("\n--- Sea Ice Chart Query Results (2025/01/03) ---")
    results_df = pd.DataFrame(results)
    #print(results_df.to_markdown(index=False))
    filename = filename_prefix + "_" + str(years[0]) + "_" + str(years[-1]) + "_sigrid3.csv"
    results_df.to_csv(filename)
    results = []
