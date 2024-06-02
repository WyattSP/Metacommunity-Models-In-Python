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
from MCME_Speciation import MCME
import MCME_Functions
import MC_Properies_Functions
import numpy as np
from MCME_Plot_Functions import mcme_plot, species_richness_plot, plot_origination_extinction
import matplotlib.pyplot as plt
import sys
import Climate_Import

### General improvements for speed
# Convert distance matrix to off-diagonal (upper or lower triangle) only, or sparse-matrix
# Having a static-global distance matrix with a resistence matrix may be better
# Resistence matrix defined by habitat niche
### You need to work through an example with 1 starting species; can't just copy old interaction matrix's

#########################
#  Set Function Inputs  #
#########################

# [Seed time, Seed Step, Burn-in Time, Simulation Time]
time_in = [200, 10, 200, 400]    
S = 25 # Starting number of species
max_growth_r = 5 # r max
M = 100 # Number of patchs

# Set interaction matrix
end_member = "equal"
alpha_matrix_ini = MCME_Functions.initalize_aij(end_member, 0.0, 1.0, S, 0.3)

# Distance matrix
#distance_matrix_ini = np.ones([M,M])
#np.fill_diagonal(distance_matrix_ini, 0.0)
# Select random coordinates
coords = MCME_Functions.random_coordinates(100)

distance_matrix = MCME_Functions.get_graph_distance_matrix_HadCM3(coords)

# General Note: Some form of area normalization will be required to normalize all coordinates across the globe. Potential normalized carrying capacity?
# Another idea is to regrid cliamte data and resample from new grid centriods

# Provide Climate
# Import Climate at a particular geographic location 
# This can either be an entire continent or a single geographic grid cell
# We are going to randomly pick 20 local communities from within North America
# Import from RDS files used in your environmental_variance work
# Tmax 4kyr
Tmax_4kyr = Climate_Import.import_climate_data("/Users/wyattpetryshen/Documents/Nature Paper 2024/Climate/Proc_Bricks/Tmax.mon.4kyr.tif", "Geotif", "rasterio")
# Climate input array
climate_input = Climate_Import.return_climate_array(coords, 0, 800, Tmax_4kyr) - 273.5

# Expand Climate input to run for desired simulation runtime

# Provide initial niche optimum
#niche_optimum_ini = np.repeat(np.random.rand(species_number_ini), M).reshape(species_number_ini, M).T 
# initial niche optimum is the intial climate +- some value
#niche_optimum_ini = np.random.rand(S)
niche_optimum_ini = MCME_Functions.initial_species_niche_optimum(S, M, climate_input[0,:])

# Provide Initial dispersal ability
# dispersal_rate_ini = np.random.choice(MCME_Functions.random_dispersal_vector(16))
dispersal_rate_ini = MCME_Functions.random_dispersal_vector(16)[14]

# Provide Initial niche breadth
# niche_breadth_ini = np.random.choice(MCME_Functions.random_sigma_vector(13))
# niche_breadth_ini = MCME_Functions.random_sigma_vector(13)[10]
niche_breadth_ini = 10

# Specify the speciation rate
speciation_rate = 0.001

#########################
#   Run Simulation Exp  #
#########################
import time

start = time.time()
N_save, lambda_save, env_save, den_save, niche_opt, alpha_matrix, phylogeny, divergence_time = MCME(time_in, 
                                                                                                    S, 
                                                                                                    max_growth_r,
                                                                                                    speciation_rate,
                                                                                                    alpha_matrix_ini,
                                                                                                    distance_matrix, 
                                                                                                    climate_input, 
                                                                                                    niche_optimum_ini, 
                                                                                                    dispersal_rate_ini, 
                                                                                                    niche_breadth_ini,
                                                                                                    end_member)
end = time.time()
print(end - start)
#########################
#      Profile          #
#########################
from cProfile import Profile
from pstats import SortKey, Stats

# When profiling set the return function to nothing so you don't overflow your console

with Profile() as profile:
    print(f"{MCME(time_in,S,max_growth_r,speciation_rate,alpha_matrix_ini,distance_matrix,climate_input, niche_optimum_ini, dispersal_rate_ini, niche_breadth_ini, end_member) = }")
    (
     Stats(profile)
     .strip_dirs()
     .sort_stats(SortKey.CALLS)
     .print_stats()
     )


sys.getsizeof(N_save)
##########################
# Calculate Comm Metrics #
##########################

gamma = MC_Properies_Functions.gamma_diversity(N_save)
alpha = MC_Properies_Functions.alpha_richness(N_save)

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
orgination, extinction = plot_origination_extinction(N_save, True)

# Plot map of geographic points with colours illustrating community richness

