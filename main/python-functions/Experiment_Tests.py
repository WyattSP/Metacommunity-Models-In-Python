#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 22 09:49:00 2024

@author: wyattpetryshen
"""

# Import functions
from MCME_Speciation import MCME
import MCME_Functions
import MC_Properies_Functions
import numpy as np
from MCME_Plot_Functions import mcme_plot
import matplotlib.pyplot as plt
import sys

#########################
#  Set Function Inputs  #
#########################
# [Seed time, Seed Step, Burn-in Time, Simulation Time]
time_in = [100, 10, 200, 500]    

species_number_ini = 50 # Starting number of species
max_growth_r = 5 # r max
M = 20 # Number of patchs

# Set interaction matrix
alpha_matrix_ini = MCME_Functions.initalize_aij("stabalizing", 0.0, 1.0, species_number_ini, 0.3)

# Distance matrix
distance_matrix_ini = np.ones([M,M])
np.fill_diagonal(distance_matrix_ini, 0.0)

# Provide Climate
# Import Climate at a particular geographic location 
# This can either be an entire continent or a single geographic grid cell
# Import from RDS files used in your environmental_variance work
climate_input = (np.random.rand(M*time_in[3]) * 1.0).reshape(M,time_in[3]) 

# Expand Climate input to run for desired simulation runtime

# Provide initial niche optimum
#niche_optimum_ini = np.repeat(np.random.rand(species_number_ini), M).reshape(species_number_ini, M).T 
niche_optimum_ini = np.random.rand(species_number_ini)

# Provide Initial dispersal ability
# dispersal_rate_ini = np.random.choice(MCME_Functions.random_dispersal_vector(16))
dispersal_rate_ini = MCME_Functions.random_dispersal_vector(16)[10]

# Provide Initial niche breadth
# niche_breadth_ini = np.random.choice(MCME_Functions.random_sigma_vector(13))
# niche_breadth_ini = MCME_Functions.random_sigma_vector(13)[10]
niche_breadth_ini = 10

# Specify the speciation rate
speciation_rate = 0.1

#########################
#   Run Simulation Exp  #
#########################

N_save, lambda_save, env_save, den_save, niche_opt, alpha_matrix, phylogeny, divergence_time = MCME(time_in, 
                                                                                                    species_number_ini, 
                                                                                                    max_growth_r,
                                                                                                    speciation_rate,
                                                                                                    alpha_matrix_ini,
                                                                                                    distance_matrix_ini, 
                                                                                                    climate_input, 
                                                                                                    niche_optimum_ini, 
                                                                                                    dispersal_rate_ini, 
                                                                                                    niche_breadth_ini)

sys.getsizeof(N_save)
##########################
# Calculate Comm Metrics #
##########################

gamma = MC_Properies_Functions.gamma_diversity(N_save)
alpha = MC_Properies_Functions.alpha_richness(N_save)

t_gamma = MC_Properies_Functions.temporal_gamma_diversity(N_save)

#########################
#          Plot         # 
#########################

x = np.arange(0,len(alpha))
plt.plot(x, gamma, label = "gamma")
plt.plot(x, alpha, c = "green", label = "alpha")
plt.legend(loc="upper right")


plt.plot(range(100), t_gamma, c = "red")


# Troubleshouting variables
S = 50
M = 100
r = 5
climate = (np.random.rand(M*time_in[3]) * 1.0).reshape(M,time_in[3]) 
z = np.repeat(np.random.rand(species_number_ini), M).reshape(species_number_ini, M).T 
sigma_niche = MCME_Functions.random_sigma_vector(13)[10]
alpha_matrix =  MCME_Functions.initalize_aij("stabalizing", 0.0, 1.0, species_number_ini, 0.3)
dispersal_rate = MCME_Functions.random_dispersal_vector(16)[10]
distance_matrix = distance_matrix_ini
