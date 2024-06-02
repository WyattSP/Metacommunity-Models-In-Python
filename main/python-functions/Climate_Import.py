#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 27 17:59:42 2024

@author: wyattpetryshen
"""

# Import climate data from .rds files 
import pyreadr # this library allows for .rds files to be imported
import rasterio # this allows for raster files to be imported
import numpy as np
#import xarray
#import rioxarray # Need to pip install this...

def import_climate_data(path_to_file, file_type, engine):
    if file_type == "RDS":
        print("Import RDS file.")
        # This may be very very unstable and will not work with S4 R objects. I do not know about S3 objects.
        # https://github.com/ofajardo/pyreadr
        file = pyreadr.read_r(path_to_file) 
        print("Import Complete")
        return(file)
    elif file_type == "Geotif" and engine == "rasterio":
        print("Importing Raster file.")
        # This is likely the most stable and reliable way to import data as a Geotif 
        file = rasterio.open(path_to_file)
        # Create empty array to store climate data
        climate_array = np.empty((file.shape[0], file.shape[1], file.count))
        # For loop to fill in empty array
        for i in np.arange(1, file.count + 1):
            in_layer = file.read(int(i))
            climate_array[:,:,i - 1] = in_layer
        print("Import Complete")
        return(climate_array)
    #elif file_type == "Geotif" and engine == "xarray":
    #    x_array_obj = rioxarray.open_rasterio(path_to_file)
    #    print("Import Complete")
    #    return(x_array_obj)
    else:
        return(print("Incorrect file type provided."))
    

def return_climate_array(coord_list, start_index, end_index, global_HadCM3_climate):
    # convert the lat lon into grid positions 
    # this assumes the 3.75 by 2.5 grid dimensions
    lat_centers = [90.0000,87.5000,85.0000,82.5000,80.0000,77.5000,75.0000,72.5000,
                 70.0000,67.5000,65.0000,62.5000,60.0000,57.5000,55.0000,52.5000,
                 50.0000,47.5000,45.0000,42.5000,40.0000,37.5000,35.0000,32.5000,
                 30.0000,27.5000,25.0000,22.5000,20.0000,17.5000,15.0000,12.5000,
                 10.0000,7.5000,5.0000,2.5000,0.0000,-2.5000,-5.0000,-7.5000,
                 -10.0000,-12.5000,-15.0000,-17.5000,-20.0000,-22.5000,-25.0000,-27.5000,
                 -30.0000,-32.5000,-35.0000,-37.5000,-40.0000,-42.5000,-45.0000,-47.5000,
                 -50.0000,-52.5000,-55.0000,-57.5000,-60.0000,-62.5000,-65.0000,-67.5000,
                 -70.0000,-72.5000,-75.0000,-77.5000,-80.0000,-82.5000,-85.0000,-87.5000,
                 -90.0000]
    lon_centers = [0.0000, 3.7500,7.5000,11.2500,15.0000,18.7500,22.5000,26.2500,
                 30.0000,33.7500,37.5000,41.2500,45.0000,48.7500,52.5000,56.2500,
                 60.0000,63.7500,67.5000,71.2500,75.0000,78.7500,82.5000,86.2500,
                 90.0000,93.7500,97.5000,101.2500,105.0000,108.7500,112.5000,116.2500,
                 120.0000,123.7500,127.5000,131.2500,135.0000,138.7500,142.5000,146.2500,
                 150.0000,153.7500,157.5000,161.2500,165.0000,168.7500,172.5000,176.2500,
                 180.0000,183.7500,187.5000,191.2500,195.0000,198.7500,202.5000,206.2500,
                 210.0000,213.7500,217.5000,221.2500,225.0000,228.7500,232.5000,236.2500,
                 240.0000,243.7500,247.5000,251.2500,255.0000,258.7500,262.5000,266.2500,
                 270.0000,273.7500,277.5000,281.2500,285.0000,288.7500,292.5000,296.2500,
                 300.0000,303.7500,307.5000,311.2500,315.0000,318.7500,322.5000,326.2500,
                 330.0000,333.7500,337.5000,341.2500,345.0000,348.7500,352.5000,356.2500]
    
    #lat_centers = np.trim_zeros(lat_centers)  
    #lon_centers = np.trim_zeros(lon_centers) 
    
    coord_list = np.trim_zeros(coord_list) 
    
    climate_out = np.zeros((end_index-start_index, len(coord_list)))
    count = 0
    for i in coord_list:
        index_lon = lon_centers.index(i[0])
        index_lat = lat_centers.index(i[1])
        
        climate_out[start_index:end_index, count] = global_HadCM3_climate[index_lat, index_lon, start_index:end_index]
        count += 1
    return(climate_out)

