#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 09:49:00 2024

@author: wyattpetryshen
"""
# Set working directory
import os
os.chdir('/Users/wyattpetryshen/Documents/GitHub/Metacommunity-Models-In-Python/main/python-functions')

# Import functions
#from MCME_Speciation import MCME
from MCME_Allopatric import MCME
import MCME_Functions
import MC_Properies_Functions
import numpy as np
from MCME_Plot_Functions import mcme_plot, species_richness_plot, plot_origination_extinction, plot_richness_map, get_richness
import MCME_Plot_Functions
import butterworth_filter
import matplotlib.pyplot as plt
import sys
import Climate_Import
import pandas as pd
import random
import math

### General improvements for speed
# Convert distance matrix to off-diagonal (upper or lower triangle) only, or sparse-matrix
# Having a static-global distance matrix with a resistence matrix may be better
# Resistence matrix defined by habitat niche
### You need to work through an example with 1 starting species; can't just copy old interaction matrix's

#########################
#  Set Function Inputs  #
#########################

Tmax_4kyr = Climate_Import.import_climate_data("/Users/wyattpetryshen/Documents/Nature Paper 2024/Climate/Proc_Bricks/Tmax.mon.4kyr.tif", "Geotif", "rasterio")


#########################
#   Global Parameters   #
#########################

S = 25 # set sarting species number
max_growth_r = 5 # set max growth rate
M = 165 # set patch number

###################################
#  Import Continental Coordinates #
###################################

# Metacommunity on North America
Cont_Coords = pd.read_csv("/Users/wyattpetryshen/Documents/GitHub/Metacommunity-Models-In-Python/main/World_Continent_Coords.csv")

# Extract only NA coords
# Keep these below lines for full dataset
mask = Cont_Coords.Cont == "South America"
Clim_array = Cont_Coords[mask].reset_index()

# Convert Data Frame into array of zip (lon,lat)
Clim_Coords = list()
for i in range(len(Clim_array)):
    Clim_Coords.append((Clim_array.Lon_d[i], Clim_array.Lat_d[i]))
    
# Randomly sample points in North America
coord_index = np.random.choice(range(len(Clim_array)), M, replace = False)
coords = [Clim_Coords[i] for i in coord_index]

# Define distance matrix
distance_matrix = MCME_Functions.get_graph_distance_matrix_HadCM3(coords)

###################################
#       Interpolate Climate       #
###################################
# Intervals: 4 kya, 50 kya, 100 kya

# Native resolution climate for desired coordinates
climate_4kya = Climate_Import.return_climate_array(coords, 400, 901, Tmax_4kyr) - 273.5

# Interpolated Climate
# Resample at new intervals
climate_16kya = Climate_Import.interpolate_climate_array(climate_4kya, 4) # skip every second points
climate_40kya = Climate_Import.interpolate_climate_array(climate_4kya, 10) # skip every fifteen points
climate_100kya = Climate_Import.interpolate_climate_array(climate_4kya, 25) # skip every fifty points

######################################
#       Fourier Smooth Climate       #
######################################
# Smooth out variance at: 4 kya, 40 kya, 100 kya
# lfreq: 
#tseries_freq = Climate_Import.get_timeseries_frequency_range(climate_4kya[:,1])
sample_rate = 4 # t series sample rate
# double sampling 4 * 3 = 12 kya
lfreq_1 = round(1 / (16 / sample_rate), 3)
# 10 times sampling 4 * 10 = 40 kya
lfreq_2 = round(1 / (40 / sample_rate), 3)
# 25 times sampling 4 * 25 = 100 kya
lfreq_3 = round(1 / (100 / sample_rate), 3)
# Get new climates
fClimate_16kya = Climate_Import.get_fourier_smooth_climate(climate_4kya[1:], lfreq_1, 0.5)
fClimate_40kya = Climate_Import.get_fourier_smooth_climate(climate_4kya[1:], lfreq_2, 0.5)
fClimate_100kya = Climate_Import.get_fourier_smooth_climate(climate_4kya[1:], lfreq_3, 0.5)

##########################################
#      Butterworth Low-Pass Filter       #
##########################################
sample_rate = 4 # t series sample rate
# double sampling 4 * 3 = 16 kya
lfreq_1 = round(1 / (16 / sample_rate), 3)
# 10 times sampling 4 * 10 = 40 kya
lfreq_2 = round(1 / (40 / sample_rate), 3)
# 25 times sampling 4 * 25 = 100 kya
lfreq_3 = round(1 / (100 / sample_rate), 3)

# Note that for evenly spaced data you will use analog = False in the butter function 
lowpass_freq = list([lfreq_1, lfreq_2, lfreq_3])

BFc16 = butterworth_filter.get_the_butter(climate_array = climate_4kya[1:], lpass = lfreq_1, order = 9)
BFc40 = butterworth_filter.get_the_butter(climate_array = climate_4kya[1:], lpass = lfreq_2, order = 9)
BFc100 = butterworth_filter.get_the_butter(climate_array = climate_4kya[1:],lpass = lfreq_3, order = 9)

###################################
#    Set up Model Parameters      #
###################################
# Max Simulation time
sim_time = climate_4kya.shape[0]

# Simulation times
time_in = [200, 10, 100, sim_time] # seed time, seed inteval, burn-in time, simulation time

# Set intial species niche optimums; this will stay static between all simulations
niche_optimum_ini_4 = MCME_Functions.initial_species_niche_optimum(S, M, climate_4kya[0,:])
niche_optimum_ini_100 = MCME_Functions.initial_species_niche_optimum(S, M, fClimate_100kya[0,:])

###################################
#     Iteratable Parameters       #
###################################

# Dispersal ability
#dispersal_rate_ini = MCME_Functions.random_dispersal_vector(16, 1e-5, 1)
dispersal_rate_ini = [0.01]

# Niche Breadth
niche_breadth = [0.5]

# Specify the speciation rate
speciation_threshold = [0]

# Species interaction end-members
end_member = ["equal", "stabalizing", "mixed", "neutral"]
# Equal species interactions
alpha_matrix_equal = MCME_Functions.initalize_aij("equal", 0.0, 1.0, S, 0.3)
# Stabalizing species interactions
alpha_matrix_stabalizing = MCME_Functions.initalize_aij("stabalizing", 0.0, 0.5, S, 0.3)
# Mixed species interactions
alpha_matrix_mixed = MCME_Functions.initalize_aij("mixed", 0.0, 1.5, S, 0.3)
# Neutral species interactions
alpha_matrix_neutral = MCME_Functions.initalize_aij("neutral", 0.0, 1.0, S, 0.3)

alpha_matrix_list = list([alpha_matrix_equal, alpha_matrix_stabalizing, alpha_matrix_mixed, alpha_matrix_neutral])


#########################
#   Run Simulation Exp  #
#########################

# 4 kya forcings
results, niche_opt, alpha_matrix, phylogeny, divergence_time, patch_origin = MCME(time_in =  [100, 10, 100, 500],
                                                                   species_number_ini = 25,
                                                                   max_growth_r = 5,
                                                                   speciation_threshold = 0.25,
                                                                   alpha_matrix_ini = alpha_matrix_stabalizing,
                                                                   distance_matrix = distance_matrix,
                                                                   climate_input = climate_4kya,
                                                                   niche_optimum_ini = niche_optimum_ini_4,
                                                                   dispersal_rate_ini = 0.1,
                                                                   niche_breadth_ini = 0.5,
                                                                   end_member = "stabalizing")


# 100 kya forcings
results_100, niche_opt_100, alpha_matrix_100, phylogeny_100, divergence_time_100, patch_origin_100 = MCME(time_in =  [100, 10, 100, 500],
                                                                   species_number_ini = 25,
                                                                   max_growth_r = 5,
                                                                   speciation_threshold = 0.25,
                                                                   alpha_matrix_ini = alpha_matrix_stabalizing,
                                                                   distance_matrix = distance_matrix,
                                                                   climate_input = fClimate_100kya,
                                                                   niche_optimum_ini = niche_optimum_ini_100,
                                                                   dispersal_rate_ini = 0.1,
                                                                   niche_breadth_ini = 0.5,
                                                                   end_member = "stabalizing")

#########################
#      Profile          #
#########################
from cProfile import Profile
from pstats import SortKey, Stats

# When profiling set the return function to nothing so you don't overflow your console

with Profile() as profile:
    print(f"{ MCME(time_in =  [100, 10, 100, 500],species_number_ini = 25,max_growth_r = 5,speciation_threshold = 0.25,alpha_matrix_ini = alpha_matrix_stabalizing,distance_matrix = distance_matrix, climate_input = fClimate_100kya,niche_optimum_ini = niche_optimum_ini_100,dispersal_rate_ini = 0.1,niche_breadth_ini = 0.5,end_member = 'stabalizing') = }")
    (
     Stats(profile)
     .strip_dirs()
     .sort_stats(SortKey.CALLS)
     .print_stats()
     )

sys.getsizeof(results)

######################################
# Notes for Performance Improvements #
######################################
import time
start = time.time()
np.random.poisson(10000,int(1e7))
end = time.time()
print(end-start)


##########################
# Calculate Comm Metrics #
##########################

gamma = MC_Properies_Functions.gamma_diversity(results_100)
alpha = MC_Properies_Functions.alpha_richness(results_100)

# t_gamma = MC_Properies_Functions.temporal_gamma_diversity(N_save)

#########################
#          Plot         # 
#########################

# Plot of gamma and alpha diversity
x = np.arange(0,len(alpha))
plt.plot(x, gamma, label = "gamma")
plt.plot(x, alpha, c = "green", label = "alpha")
plt.legend(loc="upper left")


# Plot of species richness change
orgination, extinction = plot_origination_extinction(results, True)
plt.plot(range(len(np.cumsum(orgination[:,1]))), np.cumsum(orgination[:,1]), c = "blue")
plt.plot(range(len(np.cumsum(extinction[:,1]))), np.cumsum(extinction[:,1]), c = "red")

# Plot map of geographic points with colours illustrating community richness
plot_richness_map(results[-1], Clim_array, coord_index, False)

#np.savetxt("/Users/wyattpetryshen/Documents/GitHub/Metacommunity-Models-In-Python/test_out_raster.csv", t, delimiter=",")

##########################################
#        Plots of Climate Inputs         #
##########################################
# Mean Climate Untransformed
climate_mean = np.mean(climate_4kya[1:],1)
# Mean Butterworths
Bfc16_mean = np.mean(BFc16,1)
Bfc40_mean = np.mean(BFc40,1)
Bfc100_mean = np.mean(BFc100,1)
# Mean Fourier
fC16_mean = np.mean(fClimate_16kya, 1)
fC40_mean = np.mean(fClimate_40kya, 1)
fC100_mean = np.mean(fClimate_100kya, 1)
# Mean Interpolations
Sc16_mean = np.mean(climate_16kya, 1)
Sc40_mean = np.mean(climate_40kya, 1)
Sc100_mean = np.mean(climate_100kya, 1)

##########################
#### Plot Time Domain ####
##########################
fig, axs = plt.subplots(4, 1, figsize=(12, 8))
plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None)
# Plot of Original Data
axs[0].plot(range(len(climate_mean)), climate_mean, c = "black")
axs[0].set_ylim([23, 28])
axs[0].text(0, 23.25,'Mean Untransformed Climate')
axs[0].xaxis.set_ticks([0, 100, 200, 300, 400, 500]) 
axs[0].xaxis.set_ticklabels(["2000", "1600", "1200", "800", "400", "0"])
axs[0].set_ylabel('Temperature (\u00B0C)')

# Butterworth Plot
axs[1].plot(range(len(Bfc16_mean)), Bfc16_mean)
axs[1].plot(range(len(Bfc40_mean)), Bfc40_mean)
axs[1].plot(range(len(Bfc100_mean)), Bfc100_mean)
axs[1].set_ylim([23, 28])
axs[1].text(0, 23.25,'Mean Butterworth Climate')
axs[1].xaxis.set_ticks([0, 100, 200, 300, 400, 500]) 
axs[1].xaxis.set_ticklabels(["2000", "1600", "1200", "800", "400", "0"])
axs[1].set_ylabel('Temperature (\u00B0C)')
# Fourier Plot
axs[2].plot(range(len(fC16_mean)), fC16_mean)
axs[2].plot(range(len(fC40_mean)), fC40_mean)
axs[2].plot(range(len(fC100_mean)), fC100_mean)
axs[2].set_ylim([23, 28])
axs[2].text(0, 23.25,'Mean Fourier Climate')
axs[2].xaxis.set_ticks([0, 100, 200, 300, 400, 500]) 
axs[2].xaxis.set_ticklabels(["2000", "1600", "1200", "800", "400", "0"])
axs[2].set_ylabel('Temperature (\u00B0C)')
# Smooth Interpolated Plot
axs[3].plot(range(len(Sc16_mean)), Sc16_mean)
axs[3].plot(range(len(Sc40_mean)), Sc40_mean)
axs[3].plot(range(len(Sc100_mean)), Sc100_mean)
axs[3].set_ylim([23, 28])
axs[3].text(0, 23.25,'Mean Smooth Interpolated Climate')
axs[3].xaxis.set_ticks([0, 100, 200, 300, 400, 500]) 
axs[3].xaxis.set_ticklabels(["2000", "1600", "1200", "800", "400", "0"])
axs[3].set_xlabel("Thousand Years")
axs[3].set_ylabel('Temperature (\u00B0C)')
# Legend
fig.legend(labels=['4 kya', '16 kya', '40 kya', '100 kya'], loc='lower left', 
           shadow=False, ncol=4, bbox_to_anchor=(0.075, 0.05))

fig.savefig('/Users/wyattpetryshen/Library/CloudStorage/GoogleDrive-wyatt.petryshen@yale.edu/My Drive/Conferences/WBF 2024/figures/mean_climates.png', format='png', dpi=300, transparent=False)

###############################
#### Plot Frequency Domain ####
###############################





