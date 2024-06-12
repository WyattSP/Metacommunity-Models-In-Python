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
from MCME_Plot_Functions import mcme_plot, species_richness_plot, plot_origination_extinction, plot_richness_map
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

S = 5  # set sarting species number
max_growth_r = 1 # set max growth rate
M = 5 # set patch number

###################################
#  Import Continental Coordinates #
###################################

# Metacommunity on North America
Cont_Coords = pd.read_csv("/Users/wyattpetryshen/Documents/GitHub/Metacommunity-Models-In-Python/main/World_Continent_Coords.csv")

# Extract only NA coords
# Keep these below lines for full dataset
mask = Cont_Coords.Cont == "North America"
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
climate_2kya = Climate_Import.return_climate_array(coords, 400, 901, Tmax_4kyr) - 273.5

# If you need to check output climate lengts use below function
Climate_Import.linear_climate_interpolation(climate_2kya[:,0], 25, True)

# Resample at new intervals
climate_4kya = Climate_Import.interpolate_climate_array(climate_2kya, 2) # skip every second points
climate_50kya = Climate_Import.interpolate_climate_array(climate_2kya, 25) # skip every fifteen points
climate_100kya = Climate_Import.interpolate_climate_array(climate_2kya, 50) # skip every fifty points

# Put climates into a giant array
climate_input_list = list([climate_2kya, climate_4kya, climate_50kya, climate_100kya])
climate_input_names = ["2kya_step", "4kya_step", "50kya_step", "100kya_step"]

###################################
#    Set up Model Parameters      #
###################################
# Max Simulation time
sim_time = climate_2kya.shape[0]

# Simulation times
time_in = [200, 10, 200, sim_time] # seed time, seed inteval, burn-in time, simulation time

# Set intial species niche optimums; this will stay static between all simulations
niche_optimum_ini = MCME_Functions.initial_species_niche_optimum(S, M, climate_2kya[0,:])

###################################
#     Iteratable Parameters       #
###################################

# Dispersal ability
#dispersal_rate_ini = MCME_Functions.random_dispersal_vector(16, 1e-5, 1)
dispersal_rate_ini = [0.1]

# Niche Breadth
niche_breadth = [5]

# Specify the speciation rate
speciation_rate = [0.1, 0.01]

# Species interaction end-members
end_member = ["equal", "stabalizing", "mixed", "neutral"]
# Equal species interactions
alpha_matrix_equal = MCME_Functions.initalize_aij("equal", 0.0, 1.0, S, 0.3)
# Stabalizing species interactions
alpha_matrix_stabalizing = MCME_Functions.initalize_aij("stabalizing", 0.0, 1.0, S, 0.3)
# Mixed species interactions
alpha_matrix_mixed = MCME_Functions.initalize_aij("mixed", 0.0, 1.5, S, 0.3)
# Neutral species interactions
alpha_matrix_neutral = MCME_Functions.initalize_aij("neutral", 0.0, 1.0, S, 0.3)

alpha_matrix_list = list([alpha_matrix_equal, alpha_matrix_stabalizing, alpha_matrix_mixed, alpha_matrix_neutral])


#########################
#   Run Simulation Exp  #
#########################
import time

start = time.time()
results, niche_opt, alpha_matrix, phylogeny, divergence_time = MCME(time_in =  [100, 10, 100, 200],
                                                                   species_number_ini = 5,
                                                                   max_growth_r = 2,
                                                                   max_population = 50000,
                                                                   speciation_rate = 0.001,
                                                                   alpha_matrix_ini = alpha_matrix_stabalizing,
                                                                   distance_matrix = distance_matrix,
                                                                   climate_input = climate_4kya,
                                                                   niche_optimum_ini = niche_optimum_ini,
                                                                   dispersal_rate_ini = 0.01,
                                                                   niche_breadth_ini = 10,
                                                                   end_member = "stabalizing")

end = time.time()
print(end - start)

# Test Save

np.savez("/Users/wyattpetryshen/Documents/sTime/test_results.npz", *results)
reload = np.load("/Users/wyattpetryshen/Documents/sTime/test_results.npz")
list2 = [reload[k] for k in reload]

#########################
#      Profile          #
#########################
from cProfile import Profile
from pstats import SortKey, Stats

# When profiling set the return function to nothing so you don't overflow your console

with Profile() as profile:
    print(f"{ MCME(time_in =[200, 10, 200, 200], species_number_ini = S, max_growth_r = 5, speciation_rate = 0.1, alpha_matrix_ini = alpha_matrix_stabalizing, distance_matrix = distance_matrix, climate_input = climate_input_list[1], niche_optimum_ini = niche_optimum_ini, dispersal_rate_ini = 0.001, niche_breadth_ini = 2.5, end_member = end_member[1] ) = }")
    (
     Stats(profile)
     .strip_dirs()
     .sort_stats(SortKey.CALLS)
     .print_stats()
     )


sys.getsizeof(results)
##########################
# Calculate Comm Metrics #
##########################

gamma = MC_Properies_Functions.gamma_diversity(results)
alpha = MC_Properies_Functions.alpha_richness(results)

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
plot_richness_map(results[0], Clim_array, coord_index, False)
plot_richness_map(results[-1], Clim_array, coord_index, False)

#np.savetxt("/Users/wyattpetryshen/Documents/GitHub/Metacommunity-Models-In-Python/test_out_raster.csv", t, delimiter=",")


####### 
#PARAMS for individual run
time_in =  [100, 10, 200, 200]
species_number_ini = 5
max_growth_r = 2
max_population = 50000
speciation_rate = 0.001
alpha_matrix_ini = alpha_matrix_stabalizing
distance_matrix = distance_matrix
climate_input = climate_4kya
niche_optimum_ini = niche_optimum_ini
dispersal_rate_ini = 0.5
niche_breadth_ini = 5
end_member = "stabalizing"


