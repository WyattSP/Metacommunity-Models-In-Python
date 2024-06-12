#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jun  1 14:44:48 2024

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
from MCME_Plot_Functions import mcme_plot, species_richness_plot, plot_origination_extinction, plot_richness_map
import matplotlib.pyplot as plt
import sys
import Climate_Import
import pandas as pd
import random
import math

# Large imports 
# Climate rasters
Tmax_4kyr = Climate_Import.import_climate_data("/Users/wyattpetryshen/Documents/Nature Paper 2024/Climate/Proc_Bricks/Tmax.mon.4kyr.tif", "Geotif", "rasterio")

# Metacommunity on North America
Cont_Coords = pd.read_csv("/Users/wyattpetryshen/Documents/GitHub/Metacommunity-Models-In-Python/main/World_Continent_Coords.csv")

# Extract only NA coords
# Keep these below lines for full dataset
mask = Cont_Coords.Cont == "North America"
NA_array = Cont_Coords[mask]

# Convert Data Frame into array of zip (lon,lat)
Na_Coords = list()
for i in range(len(NA_array)):
    Na_Coords.append((NA_array.Lon_d[i], NA_array.Lat_d[i]))
    
# Randomly sample points in North America
coord_index = [random.choice(range(len(NA_array))) for i in range(10)]
coords = [Na_Coords[i] for i in coord_index]

# coords = MCME_Functions.random_coordinates(100) # Actual random global coordinates
# Calculate distance matrix
distance_matrix = MCME_Functions.get_graph_distance_matrix_HadCM3(coords)

climate_input = Climate_Import.return_climate_array(coords, 0, 800, Tmax_4kyr) - 273.5

# Try to recreate figure 3 from Thompson et al. 2020.
# Set intial paramters 
time_in = [200, 10, 200, 800]    
S = 25; max_growth_r = 5; M = 10 

end_member = "equal"; alpha_matrix_ini = MCME_Functions.initalize_aij(end_member, 0.0, 1.0, S, 0.3)

niche_optimum_ini = MCME_Functions.initial_species_niche_optimum(S, M, climate_input[0,:])

# Provide Initial dispersal ability
# dispersal_rate_ini = np.random.choice(MCME_Functions.random_dispersal_vector(16))
dispersal_rate_ini = MCME_Functions.random_dispersal_vector(16)
niche_breadth = [0.5, 10]

# Specify the speciation rate
speciation_rate = 0.001

# Save list for each experiment
MCME_niche0_5_equal = list()
total = len(dispersal_rate_ini)

counter = 0
for dis in dispersal_rate_ini:
    N_save, lambda_save, env_save, den_save, niche_opt, alpha_matrix, phylogeny, divergence_time = MCME(time_in, 
                                                                                                        S, 
                                                                                                        max_growth_r,
                                                                                                        speciation_rate,
                                                                                                        alpha_matrix_ini,
                                                                                                        distance_matrix, 
                                                                                                        climate_input, 
                                                                                                        niche_optimum_ini, 
                                                                                                        dis, 
                                                                                                        niche_breadth[0],
                                                                                                        end_member)
    counter += 1
    print("Completed dispersal: %d" % counter, "of %d" % total) 
    MCME_niche0_5_equal.append(N_save)
    
    
MCME_niche10_equal = list()
counter = 0
for dis in dispersal_rate_ini:
    N_save, lambda_save, env_save, den_save, niche_opt, alpha_matrix, phylogeny, divergence_time = MCME(time_in, 
                                                                                                        S, 
                                                                                                        max_growth_r,
                                                                                                        speciation_rate,
                                                                                                        alpha_matrix_ini,
                                                                                                        distance_matrix, 
                                                                                                        climate_input, 
                                                                                                        niche_optimum_ini, 
                                                                                                        dis, 
                                                                                                        niche_breadth[1],
                                                                                                        end_member)
    counter += 1
    print("Completed dispersal: %d" % counter, "of %d" % total) 
    MCME_niche10_equal.append(N_save)
    
    
# Calculate gamma and alpha

sigma10_gamma = list()
sigma10_alpha = list()

end_slice = len(MCME_niche10_equal[0]) - 1

for i, output in enumerate(MCME_niche10_equal):
    t_gamma = MC_Properies_Functions.gamma_diversity(output)[-1]
    t_alpha = MC_Properies_Functions.alpha_richness(output)[-1]
    sigma10_gamma.append(t_gamma)
    sigma10_alpha.append(t_alpha)
    

# Sigma 0.5
sigma05_gamma = list()
sigma05_alpha = list()

end_slice = len(MCME_niche0_5_equal[0]) - 1

for i, output in enumerate(MCME_niche0_5_equal):
    t_gamma = MC_Properies_Functions.gamma_diversity(output)[-1]
    t_alpha = MC_Properies_Functions.alpha_richness(output)[-1]
    sigma05_gamma.append(t_gamma)
    sigma05_alpha.append(t_alpha)
    


#########################
#          Plot         # 
#########################

# Not Log
# Plot of gamma and alpha diversity
plt.plot(dispersal_rate_ini, sigma10_gamma, label = "gamma")
plt.plot(dispersal_rate_ini, sigma10_alpha, c = "green", label = "alpha")
plt.legend(loc="upper left")

# Plot of gamma and alpha diversity
plt.plot(dispersal_rate_ini, sigma05_gamma, label = "gamma")
plt.plot(dispersal_rate_ini, sigma05_alpha, c = "green", label = "alpha")
plt.legend(loc="upper left")


# Log
plt.plot([math.log(i) for i in dispersal_rate_ini], sigma10_gamma, label = "gamma")
plt.plot([math.log(i) for i in dispersal_rate_ini], sigma10_alpha, c = "green", label = "alpha")
plt.legend(loc="upper left")


plt.plot([math.log(i) for i in dispersal_rate_ini], sigma05_gamma, label = "gamma")
plt.plot([math.log(i) for i in dispersal_rate_ini], sigma05_alpha, c = "green", label = "alpha")
plt.legend(loc="upper left")

#########################
#     Spatial Plot      #
#########################

plot_richness_map(MCME_niche10_equal[14][-1], NA_LatLon, coord_index)


