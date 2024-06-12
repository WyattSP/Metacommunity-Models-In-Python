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
import matplotlib.pyplot as plt
from scipy import signal
import scipy
#import 
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
        
        climate_out[:, count] = global_HadCM3_climate[index_lat, index_lon, start_index:end_index]
        count += 1
    return(climate_out)

###################################
#   Climate Alteration Functions  #
###################################


# Two potential methods:
    # 1) skipping every n-number of points and connecting
    # 2) using a smoothing spline from scipy 

# Skipping n-points to remove high frequency variation

def linear_climate_interpolation(climate, sampling_interval, supress_warnings = False):
    climate_in = climate

    # Original step vector
    kya_2_step = np.arange(len(climate_in))
    
    # Sample every second point
    step_interval = kya_2_step[0::sampling_interval]
    
    # Index original time-series
    new_clim = climate_in[step_interval]
    
    interp_climate = np.array([])
    # New climate interpolated
    for i in range(len(step_interval) - 1):
        # Set point values X
        s_index = step_interval[i]
        e_index = step_interval[i + 1]
        # Set climate values Y
        cs_index = new_clim[i]
        ce_index = new_clim[i + 1]
        
        # New point indexes
        new_i =  np.zeros(sampling_interval-1)
        next_i = s_index + 1
        for i in range(len(new_i)):
            new_i[i] = next_i
            next_i += 1
        
        # Get new points using a linear interpolation
        new_point = np.interp(new_i, [s_index,e_index], [cs_index,ce_index])
        
        # Combine point into array
        in_stack = np.hstack([cs_index, new_point])
        
        # Add to master array
        interp_climate = np.hstack([interp_climate,in_stack])
    
    # Add end point
    end_c = new_clim[-1]
    interp_climate = np.hstack([interp_climate,end_c])
    
    # Check lengths and throw warning if not equal to original
    if len(kya_2_step) != len(interp_climate):
        if supress_warnings == True:
            print("Warning: Climte input and function output lengths are different %d" % len(kya_2_step), "versus %d" % len(interp_climate))
    
    return(interp_climate)

def interpolate_climate_array(climate_array, sampling_interval):
    # Define new empty array
    new_array = np.zeros(climate_array.shape)
    
    # Sample each column 
    for cols in range(climate_array.shape[1]):
        new_vec = linear_climate_interpolation(climate_array[:,cols], sampling_interval)
        new_array[:,cols] = new_vec
        
    return(new_array)

def fourier_smooth_climate(climate_array, lfreq, hfreq, plot_spec_matplot = True, plot_scipy = False, verbose = False):
    # Retrieve current t-series length
    spectra = climate_array
    n = len(spectra)
    # Check that time series is even length
    if n % 2 != 0:
        # drop first element to make even
        spectra = spectra[1:]
        n = len(spectra)
    else:
        if verbose == True:
            print("Spectra is even")
        
    # Set frequency interval
    lfreq = lfreq
    hfreq = hfreq
    
    #### Scipy implementation with rfft, rfftfreq, irfft
    # fourier coefficients
    rfft = scipy.fft.rfft(spectra)
    # frequencies
    rfreq = scipy.fft.rfftfreq(n)
    # get mean
    ###mean = abs(rfft[0])
    # get Nyquist
    ####nyq = np.sqrt(np.imag(rfft[-1])**2 + np.real(rfft[-1])**2)
    # get phases
    phase = np.arctan2(np.imag(rfft), np.real(rfft))
    # get amplitudes
    amp = np.sqrt(np.imag(rfft)**2 + np.real(rfft)**2)
    # power
    powr = amp**2
    #### fourier manipulation
    # find index positions to reset
    zero_freq_index = np.where((rfreq >= lfreq) & (rfreq <= hfreq))[0]
    # Change frequencies to zero; phases don't matter
    amp[zero_freq_index] = 0
    #phase[zero_freq_index] = 0 
    # recompile time series 
    new_series = np.array([complex(amp[i]*np.cos(phase[i]), amp[i]*np.sin(phase[i])) for i in range(len(rfreq))])
    # inverse rfft
    X = scipy.fft.irfft(new_series)
    
    # FFT via numpy
    #s_fft = np.fft.fft(spectra)
    #mean = abs(s_fft[0])
    #nyq = np.sqrt(np.imag(s_fft[n//2])**2 + np.real(s_fft[n//2])**2)
    #phase = np.arctan2(np.imag(s_fft), np.real(s_fft))
    #s_fft_l = s_fft[1:n//2]
    #amp = np.sqrt(np.imag(s_fft_l)**2 + np.real(s_fft_l)**2)
    #powr = np.real(s_fft_l)**2
    #phsr = np.imag(s_fft_l)
    #freq = np.array([i/n for i in range(1, n//2)])
    #zero_freq_index = np.where((freq >= lfreq) & (freq <= hfreq))[0]
    #amp[zero_freq_index] = 0
    #new_aseries = np.concatenate(([mean], amp, [nyq], amp[::-1]))
    #new_pseries = np.concatenate(([0], phase[1:n//2-1], [0], phase[n//2:]))
    #new_series = np.array([complex(new_aseries[i]*np.cos(new_pseries[i]), new_aseries[i]*np.sin(new_pseries[i])) for i in range(len(new_aseries))])
    #X = np.real(np.fft.ifft(new_series)) / n
    
    # Check Lengths
    L = len(X)
    if L != n:
        print("Old and new series length do not match")

    # Plot if desired
    if plot_spec_matplot:
        fig, axs = plt.subplots(2, 2, figsize=(12, 8))
        axs[0, 0].plot(range(len(climate_array)), climate_array, color='black', linestyle='-', label='Original Time Series')
        axs[0, 0].set_xlabel('Observations')
        axs[0, 0].set_ylabel('Unit')
        axs[0, 0].set_title('Original Time Series')

        axs[0, 1].plot(range(len(X)), X, color='red', linestyle='-', label='Modified Time Series')
        axs[0, 1].set_xlabel('Observations')
        axs[0, 1].set_ylabel('Unit')
        axs[0, 1].set_title('Modified Time Series')

        axs[1, 0].plot(rfreq, np.log(powr), color='blue', linestyle='-')
        axs[1, 0].set_xlim([-0.02, 0.52])
        axs[1, 0].set_xlabel('frequency [Hz]')
        axs[1, 0].set_ylabel('Log PSD [V**2/Hz]')

        axs[1, 1].plot(rfreq, np.log(amp**2), color='green', linestyle='-')
        axs[1, 1].set_xlim([-0.02, 0.52])
        axs[1, 1].set_xlabel('frequency [Hz]')
        axs[1, 1].set_ylabel('Log PSD [V**2/Hz]')

        plt.tight_layout()
        plt.show()
    if plot_scipy:
        f, Pxx_den = signal.periodogram(X, 1)
        plt.semilogy(f, Pxx_den)
        #plt.ylim([1e-7, 1e2])
        plt.xlabel('frequency [Hz]')
        plt.ylabel('PSD [V**2/Hz]')
        plt.show()  

    return X
    
def get_timeseries_frequency_range(climate_array):
    spectra = climate_array
    n = len(spectra)
    # Check that time series is even length
    if n % 2 != 0:
        # drop first element to make even
        spectra = spectra[1:]
        n = len(spectra)
        print("Spectra is odd: dropped first observation at index = 0")
    else:
        print("Spectra is even")
        
    rfreq = scipy.fft.rfftfreq(n)
    return(rfreq)
        
def get_fourier_smooth_climate(climate_array, lfreq, hfreq,):
    # Define new empty array
    new_array = np.zeros(climate_array.shape)
    # Sample each column 
    for cols in range(climate_array.shape[1]):
        new_vec = fourier_smooth_climate(climate_array[:,cols], lfreq, hfreq, False, False)
        new_array[:,cols] = new_vec
        
    return(new_array)
    
def year_to_frequency(year, measurement_unit):
    return(1/(year/measurement_unit))
    

